
#include "cvode_config.hpp"
/* RHS 7-bucket + forcing I/O wall-clock timer dumps.
 * Header is empty under default (SHUD_ENABLE_DIAGNOSTICS undefined). */
#include "../Model/MD_diagnostics.hpp"

#include <errno.h>   /* errno (SHUD_SPGMR_MAXL strtol parse) */
#include <stdlib.h>  /* getenv, exit, EXIT_FAILURE (SHUD_LINSOL) */
#include <string.h>  /* strcmp (SHUD_LINSOL) */
#include <dlfcn.h>   /* dlopen / dlclose (Hypre runtime probe) */

/* Experimental BoomerAMG (algebraic multigrid) linear-solver wrapper,
 * selected at run time with SHUD_LINSOL=amg; off by default. */
#include "sunlinsol_hypre.h"

/* Runtime linear-solver selector enum + factory dispatch. The
 * `SHUD_LINSOL` env var picks the linear solver passed to CVODE.
 * Unset or `=spgmr` (the default) calls
 * SUNLinSol_SPGMR(udata, PREC_NONE, get_spgmr_maxl_from_env(), sunctx),
 * so default results do not depend on the selector. `=amg` opts in to
 * the experimental SUNLinSol_Hypre wrapper instead. Any other value
 * fatal-exits BEFORE `CVodeCreate`, so no CVODE state is allocated
 * on the error path. */
typedef enum {
    LINSOL_SPGMR = 0,
    LINSOL_AMG = 1,
    LINSOL_UNKNOWN = -1
} linsol_t;

/* `source` reports whether the selection came from the env var
 * ("env") or fell through to the default ("default"). Used to tag
 * the `[shud] linsol=...` stdout marker for log post-mortem
 * traceability (matches existing `[CVODE] SPGMR maxl=...` pattern). */
struct linsol_selection {
    linsol_t sel;
    const char *sel_name;
    const char *source;
};

static linsol_selection parse_linsol_env(void)
{
    const char *env = getenv("SHUD_LINSOL");

    /* Unset / empty / whitespace-only -> default SPGMR (no env tag). */
    if (env == NULL || env[0] == '\0') {
        return linsol_selection{LINSOL_SPGMR, "spgmr", "default"};
    }
    /* Strip leading/trailing whitespace via local pointer math.
     * Case-sensitive compare ("amG" / "AMG" / "Amg" all rejected). */
    const char *start = env;
    while (*start == ' ' || *start == '\t' || *start == '\n') ++start;
    if (*start == '\0') {
        return linsol_selection{LINSOL_SPGMR, "spgmr", "default"};
    }
    /* Compute end without modifying the env-var string. */
    const char *end = start + strlen(start);
    while (end > start && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n')) --end;
    size_t len = (size_t)(end - start);

    if (len == 5 && strncmp(start, "spgmr", 5) == 0) {
        return linsol_selection{LINSOL_SPGMR, "spgmr", "env"};
    }
    if (len == 3 && strncmp(start, "amg", 3) == 0) {
        return linsol_selection{LINSOL_AMG, "amg", "env"};
    }

    /* Unrecognized value — fatal-exit BEFORE CVodeCreate: the stderr
     * message names the offender + the accepted set, the process
     * exits non-zero, and no CVODE state is left allocated. */
    fprintf(stderr,
            "[shud] FATAL: SHUD_LINSOL=%s unrecognized; accepted: spgmr, amg\n",
            env);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

/* Factory: SPGMR backend (the default path). Wraps the
 * SUNLinSol_SPGMR call with PREC_NONE + get_spgmr_maxl_from_env() +
 * sunctx; `y` is the state N_Vector that SetCVODE receives as
 * `udata`. */
static SUNLinearSolver create_spgmr_ls(N_Vector y, SUNContext sunctx);

/* Factory: AMG backend via the SUNLinSol_Hypre wrapper. The wrapper
 * constructor returns NULL + a stderr message on invalid arguments
 * or allocation failure; a NULL return is fatal (the factory exits
 * rather than falling back to SPGMR). */
static SUNLinearSolver create_amg_ls(N_Vector y, Model_Data *MD,
                                     SUNContext sunctx);

int check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;
    
    /* Check if SUNDIALS function returned NULL pointer - no memory
     * allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        myexit(ERRCVODE);
        
    }
    /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *)flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            myexit(ERRCVODE);
        }
    }
    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        myexit(ERRCVODE);
    }
    return (0);
}
void PrintFinalStats(void *cvode_mem, FILE *fout)
{
    long int lenrw, leniw;
    long int lenrwLS, leniwLS;
    long int nst, nfe, nsetups, nni, ncfn, netf;
    long int nli, npe, nps, ncfl, nfeLS;
    int flag;

    flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
    check_flag(&flag, "CVodeGetWorkSpace", 1);
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(&flag, "CVodeGetNumSteps", 1);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_flag(&flag, "CVodeGetNumRhsEvals", 1);
    flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_flag(&flag, "CVodeGetNumErrTestFails", 1);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

//    flag = CVSpilsGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    flag = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    check_flag(&flag, "CVSpilsGetWorkSpace", 1);
//    flag = CVSpilsGetNumLinIters(cvode_mem, &nli);
    flag = CVodeGetNumLinIters(cvode_mem, &nli);
    check_flag(&flag, "CVSpilsGetNumLinIters", 1);
//    flag = CVSpilsGetNumPrecEvals(cvode_mem, &npe);
    flag = CVodeGetNumPrecEvals(cvode_mem, &npe);
    check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
//    flag = CVSpilsGetNumPrecSolves(cvode_mem, &nps);
    flag = CVodeGetNumPrecSolves(cvode_mem, &nps);
    check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
//    flag = CVSpilsGetNumConvFails(cvode_mem, &ncfl);
    flag = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
    check_flag(&flag, "CVSpilsGetNumConvFails", 1);
//    flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
    flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
    check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);

    printf("\nFinal Statistics.. \n\n");
    printf("lenrw   = %5ld     leniw   = %5ld\n", lenrw, leniw);
    printf("lenrwLS = %5ld     leniwLS = %5ld\n", lenrwLS, leniwLS);
    printf("nst     = %5ld\n", nst);
    printf("nfe     = %5ld     nfeLS   = %5ld\n", nfe, nfeLS);
    printf("nni     = %5ld     nli     = %5ld\n", nni, nli);
    printf("nsetups = %5ld     netf    = %5ld\n", nsetups, netf);
    printf("npe     = %5ld     nps     = %5ld\n", npe, nps);
    printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);

    /* Optional key=value persistence. stdout output above is the same
     * with or without `fout`, so log-scraping tools keep working.
     * Field order is deterministic + machine-parseable; the six main
     * solver counters (nfe, nfeLS, nni, nli, nsetups, netf) lead the
     * file, the remainder follow for completeness. */
    if (fout != NULL) {
        fprintf(fout, "nfe=%ld\n",     nfe);
        fprintf(fout, "nfeLS=%ld\n",   nfeLS);
        fprintf(fout, "nni=%ld\n",     nni);
        fprintf(fout, "nli=%ld\n",     nli);
        fprintf(fout, "nsetups=%ld\n", nsetups);
        fprintf(fout, "netf=%ld\n",    netf);
        fprintf(fout, "nst=%ld\n",     nst);
        fprintf(fout, "npe=%ld\n",     npe);
        fprintf(fout, "nps=%ld\n",     nps);
        fprintf(fout, "ncfn=%ld\n",    ncfn);
        fprintf(fout, "ncfl=%ld\n",    ncfl);
        fprintf(fout, "lenrw=%ld\n",   lenrw);
        fprintf(fout, "leniw=%ld\n",   leniw);
        fprintf(fout, "lenrwLS=%ld\n", lenrwLS);
        fprintf(fout, "leniwLS=%ld\n", leniwLS);
#ifdef SHUD_ENABLE_DIAGNOSTICS
        /* Diagnostics-only additions: last step size and last method
         * order, which together with nst / nfe / netf / nni / nli
         * above give the seven CVODE statistics the diagnostics
         * build reports.
         *
         * `hlast` / `qlast` are gated behind SHUD_ENABLE_DIAGNOSTICS so
         * the default build always emits exactly the 15 keys above.
         * Both are SUNDIALS 6.0.0 public API reads (post-solve, no
         * RHS path mutation), so a diagnostics build produces the
         * same model output as a default build; only this file's
         * trailing key set differs. */
        int qlast;
        realtype hlast;
        flag = CVodeGetLastStep(cvode_mem, &hlast);
        check_flag(&flag, "CVodeGetLastStep", 1);
        flag = CVodeGetLastOrder(cvode_mem, &qlast);
        check_flag(&flag, "CVodeGetLastOrder", 1);
        fprintf(fout, "hlast=%.17g\n", (double)hlast);
        fprintf(fout, "qlast=%d\n",    qlast);

        /* RHS 7-bucket + forcing I/O wall-clock dump. The seven
         * buckets cover the whole RHS, so the pct_rhs_* values sum
         * to 100% up to rounding. Buckets 0-6 are defined in
         * MD_diagnostics.hpp; their accumulators live in
         * MD_rhs_core.cpp (g_rhs_timer_ns) and TimeSeriesData.cpp
         * (g_forcing_io_ns). All values are integer nanoseconds — no
         * floating-point arithmetic in the timer path itself; the
         * percentage / seconds renders below are post-run diagnostics
         * and never feed back into RHS state. */
        long long t_total_ns = 0;
        for (int i = 0; i < shud_diag::RHS_BUCKET_COUNT; ++i) {
            t_total_ns += shud_diag::g_rhs_timer_ns[i];
        }
        fprintf(fout, "t_rhs_update=%lld\n",  shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_UPDATE]);
        fprintf(fout, "t_rhs_ET=%lld\n",      shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_ET]);
        fprintf(fout, "t_rhs_lateral=%lld\n", shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_LATERAL]);
        fprintf(fout, "t_rhs_segment=%lld\n", shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_SEGMENT]);
        fprintf(fout, "t_rhs_river=%lld\n",   shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_RIVER]);
        fprintf(fout, "t_rhs_gather=%lld\n",  shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_GATHER]);
        fprintf(fout, "t_rhs_applyDY=%lld\n", shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_APPLYDY]);
        fprintf(fout, "t_rhs_total=%lld\n",   t_total_ns);
        /* Percentages — division by 0 only possible if no RHS call was
         * ever made, which would mean the run never started; emit 0.0
         * defensively rather than NaN. */
        double total_d = (t_total_ns > 0) ? (double)t_total_ns : 1.0;
        fprintf(fout, "pct_rhs_update=%.3f\n",  100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_UPDATE]  / total_d);
        fprintf(fout, "pct_rhs_ET=%.3f\n",      100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_ET]      / total_d);
        fprintf(fout, "pct_rhs_lateral=%.3f\n", 100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_LATERAL] / total_d);
        fprintf(fout, "pct_rhs_segment=%.3f\n", 100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_SEGMENT] / total_d);
        fprintf(fout, "pct_rhs_river=%.3f\n",   100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_RIVER]   / total_d);
        fprintf(fout, "pct_rhs_gather=%.3f\n",  100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_GATHER]  / total_d);
        fprintf(fout, "pct_rhs_applyDY=%.3f\n", 100.0 * (double)shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_APPLYDY] / total_d);
        /* Forcing I/O — separate channel, NOT included in t_rhs_total. */
        fprintf(fout, "t_forcing_io_ns=%lld\n", shud_diag::g_forcing_io_ns);
        fprintf(fout, "t_forcing_io_s=%.3f\n", (double)shud_diag::g_forcing_io_ns / 1.0e9);
#endif
    }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

void CVODEstatus(void *cvode_mem, N_Vector u, realtype t){
  long int nst;
  int qu, retval;
  realtype hu, *udata;
//  udata = N_VGetArrayPointer(u);

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_flag(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_flag(&retval, "CVodeGetLastStep", 1);
  printf("t = %.2f   no. steps = %ld   order = %d   stepsize = %.4f\n",
         t, nst, qu, hu);
}


//void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS){
//
//    int flag;
//    /* allocate memory for solver */
//    /********* SUNDIALS 5.0+ ************/
//    //    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON); //v3.x
//    cvode_mem = CVodeCreate(CV_BDF);
//    check_flag((void *)cvode_mem, "CVodeCreate", 0);
//
//    flag = CVodeSetUserData(cvode_mem, MD);
//    check_flag(&flag, "CVodeSetUserData", 1);
//
//    //Model start from TIME = zero;
//    flag = CVodeInit(cvode_mem, f, MD->CS.StartTime, udata);
//    check_flag(&flag, "CVodeInit", 1);
//
//    flag = CVodeSStolerances(cvode_mem, MD->CS.reltol, MD->CS.abstol);
//    check_flag(&flag, "CVodeSStolerances", 1);
//
//    //    LS = SUNSPGMR(udata, 0, 0); //v3.x
//    LS = SUNLinSol_SPGMR(udata, 0, 0);
//    check_flag((void *)LS, "SUNLinSol_SPGMR", 0);
//
//    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
//    check_flag(&flag, "CVSpilsSetLinearSolver", 1);
//
//    flag = CVodeSetMinStep(cvode_mem, 1E-6); //Minimum time interval in cvode.dt = t(i) - t(i - 1);
//    check_flag(&flag, "CVodeSetMinStep", 1);
//
//    flag = CVodeSetMaxNumSteps(cvode_mem, 1E6); //max iterations.
//    check_flag(&flag, "CVodeSetMaxNumSteps", 1);
//
//    flag = CVodeSetInitStep(cvode_mem, MD->CS.InitStep);
//    check_flag(&flag, "CVodeSetInitStep", 1);
//
//    //force cvode run at least every x time - units.t(i) - t(i - 1) < X;
//    flag = CVodeSetMaxStep(cvode_mem, MD->CS.MaxStep);
//    check_flag(&flag, "CVodeSetMaxStep", 1);
//
////    flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
//    //flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
//}


/* Runtime env-var hook for the SPGMR Krylov subspace dimension `maxl`
 * passed to SUNLinSol_SPGMR in create_spgmr_ls(). Unset / "" / "0" / "5"
 * all give the SUNDIALS default behavior (SUNDIALS docs: maxl <= 0
 * collapses to the documented default 5; explicit 5 is the default).
 * Opt-in values {10, 15, 20, 30} flow through to SUNDIALS unchanged.
 * Every non-zero accepted value (5 included) also emits a stdout
 * provenance line so a run's log records the maxl it used. Any other
 * value (e.g. "7", "50", "foo", "-1") aborts via myexit(ERRCVODE)
 * BEFORE any SPGMR allocation, so a mistyped setting in a batch of
 * runs fails fast instead of silently running with a different maxl. */
static int get_spgmr_maxl_from_env(void)
{
    const char *env = getenv("SHUD_SPGMR_MAXL");
    /* Unset or empty -> SUNDIALS default (silent; no provenance log). */
    if (env == NULL || env[0] == '\0') {
        return 0;
    }

    /* Strict whitelist parse: reject leading whitespace ("\t 5"), leading
     * signs ("+5", "-1"), leading zeros ("05"), trailing whitespace ("5 "),
     * non-numeric suffixes ("5x", "10.0"), and any value not in the
     * allow-list. strtol(3) on its own accepts +/whitespace/leading-zeros
     * so we pre-validate the raw string character-by-character. */
    int valid_chars = 1;
    for (const char *p = env; *p != '\0'; ++p) {
        if (*p < '0' || *p > '9') { valid_chars = 0; break; }
    }
    /* Reject leading-zero forms like "05" while still permitting the
     * single character "0". env[0] != '\0' is guaranteed above. */
    int no_leading_zero = (env[0] != '0' || env[1] == '\0');
    char *endptr = NULL;
    errno = 0;
    long val = strtol(env, &endptr, 10);
    int parse_ok = (errno == 0 && endptr != NULL && *endptr == '\0' && endptr != env);
    int value_ok = valid_chars && no_leading_zero && parse_ok &&
                   (val == 0 || val == 5 || val == 10 ||
                    val == 15 || val == 20 || val == 30);
    if (!value_ok) {
        fprintf(stderr,
                "[CVODE] ERROR: SHUD_SPGMR_MAXL must be unset, 0, 5, 10, 15, 20, or 30 (got: %s)\n",
                env);
        myexit(ERRCVODE);
    }

    /* "0" is documented-equivalent to unset per SUNDIALS 6.0.0 (maxl <= 0
     * -> default 5). The default stays silent: unset / "" / "0" all
     * suppress the provenance log line. */
    if (val == 0) {
        return 0;
    }

    /* val in {5, 10, 15, 20, 30}: emit provenance line so each run's
     * log can be attributed to its maxl value. */
    fprintf(stdout, "[CVODE] SPGMR maxl=%ld pretype=PREC_NONE\n", val);
    fflush(stdout);
    return (int)val;
}

static SUNLinearSolver create_spgmr_ls(N_Vector y, SUNContext sunctx)
{
    /* Default solver: unpreconditioned SPGMR (PREC_NONE + maxl from
     * the env hook + sunctx). Changing these arguments changes the
     * default model results. */
    SUNLinearSolver LS = SUNLinSol_SPGMR(y, PREC_NONE,
                                         get_spgmr_maxl_from_env(),
                                         sunctx);
    check_flag((void *)LS, "SUNLinSol_SPGMR", 0);
    return LS;
}

static SUNLinearSolver create_amg_ls(N_Vector y, Model_Data *MD,
                                     SUNContext sunctx)
{
    /* Hardcoded (interp_type=6, coarsen_type=8). The wrapper
     * constructor rejects other pairs with stderr + NULL return. */
    SUNLinearSolver LS = SUNLinSol_Hypre(y, (void *)MD, 6, 8, sunctx);
    if (LS == NULL) {
        fprintf(stderr,
                "[shud] FATAL: AMG factory failed; check Hypre install\n");
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    return LS;
}

void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS, SUNContext &sunctx){

    int flag;

    /* Linear-solver selection is parsed BEFORE CVodeCreate so that
     * an unrecognized SHUD_LINSOL value exits with no CVODE state
     * allocated. The marker line precedes any CVODE output so logs
     * show the selection first. */
    const linsol_selection sel = parse_linsol_env();
    fprintf(stdout, "[shud] linsol=%s source=%s\n",
            sel.sel_name, sel.source);
    fflush(stdout);

    /* Hypre runtime probe (defense in depth), AMG path only.
     *
     * We attempt dlopen with multiple candidate library names; if any
     * succeeds (i.e., the dynamic loader can find Hypre via the binary
     * DT_NEEDED entries OR via LD_LIBRARY_PATH OR via the RPATH baked
     * at link time), we are confident the AMG path will work. Without
     * RTLD_NOLOAD this also forces a fresh load attempt for libraries
     * that may have been resolved lazily on macOS (where the loader's
     * versioned-soname matching is less forgiving than glibc's).
     *
     * Candidate names cover the platform matrix:
     *   - libHYPRE.so              (Linux unversioned)
     *   - libHYPRE.dylib           (Mac unversioned)
     *   - libHYPRE.so.<N>          (Linux soname; Hypre 3.x = libHYPRE.so.0)
     *   - libHYPRE.301.dylib       (Mac versioned, Hypre 3.1.0 brew)
     *   - libHYPRE-3.1.0.so        (Linux versioned-soname form,
     *                               e.g. a CMake build from source)
     *   - libHYPRE.3.1.0.dylib     (Mac versioned-soname form for
     *                               Hypre 3.1.0 brew install)
     *
     * The versioned-soname entries are required because some installs
     * (e.g. CMake-built Hypre 3.1.0 on Linux) emit a directly-
     * loadable `libHYPRE-3.1.0.so` filename WITHOUT a `libHYPRE.so`
     * symlink. Probing only the unversioned names misses such installs
     * even though the binary itself links them via DT_NEEDED.
     *
     * If a candidate dlopen succeeds, we dlclose immediately — that
     * just drops our extra reference, the loader keeps it pinned via
     * the binary's own NEEDED entry. */
    if (sel.sel == LINSOL_AMG) {
        /* Candidate order: most-common-success first, to minimize the
         * number of negative dlopen probes. A Linux CMake build of
         * Hypre 3.1.0 from source emits libHYPRE-3.1.0.so, and the Mac
         * brew install emits libHYPRE.301.dylib. Unversioned and
         * SONAME-0 fallbacks follow for other installs. */
        const char *candidates[] = {
            "libHYPRE-3.1.0.so",       /* Linux Hypre 3.1.0 built from source */
            "libHYPRE.301.dylib",      /* Mac brew Hypre 3.1.0 */
            "libHYPRE.so",             /* Linux unversioned */
            "libHYPRE.dylib",          /* Mac unversioned */
            "libHYPRE.so.0",           /* Linux SONAME-0 (Hypre 2.x) */
            "libHYPRE.3.1.0.dylib",    /* Mac alternate versioned-soname */
            NULL
        };
        void *hypre_handle = NULL;
        for (int i = 0; candidates[i] != NULL; ++i) {
            hypre_handle = dlopen(candidates[i], RTLD_LAZY);
            if (hypre_handle != NULL) break;
        }
        if (hypre_handle == NULL) {
            const char *hypre_libdir = getenv("HYPRE_LIBDIR");
            fprintf(stderr,
                    "[shud] FATAL: Hypre runtime dylib not loadable; "
                    "tried { libHYPRE-3.1.0.so, libHYPRE.301.dylib, "
                    "libHYPRE.so, libHYPRE.dylib, libHYPRE.so.0, "
                    "libHYPRE.3.1.0.dylib }. "
                    "Set LD_LIBRARY_PATH (Linux) or DYLD_LIBRARY_PATH "
                    "(macOS) to HYPRE_LIBDIR=%s or platform default "
                    "(Linux: /usr/lib, /usr/local/lib; "
                    "macOS: $(brew --prefix hypre)/lib).\n",
                    hypre_libdir ? hypre_libdir : "<unset>");
            fflush(stderr);
            exit(EXIT_FAILURE);
        }
        dlclose(hypre_handle);
    }

    /********* SUNDIALS 6.0+ ************/
    /* Allocate memory, and set problem data, initial values, tolerances */
//    u = N_VNew_Serial(NY, sunctx);
//    check_flag(void *)u, "N_VNew_Serial", 0));
//    data = AllocUserData();
//    check_flag(void *)data, "AllocUserData", 2);
//    InitUserData(data);
//    SetInitialProfiles(u, data->dx, data->dy);

    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    check_flag((void *)cvode_mem, "CVodeCreate", 0);

    flag = CVodeSetUserData(cvode_mem, MD);
    check_flag(&flag, "CVodeSetUserData", 1);

    //Model start from TIME = zero;
    flag = CVodeInit(cvode_mem, f, MD->CS.StartTime, udata);
    check_flag(&flag, "CVodeInit", 1);

    /* Optional CVODE relative tolerance override for tolerance
     * sensitivity experiments. When SHUD_CVODE_RELTOL is set to a
     * value in (0, 1), it overrides MD->CS.reltol (parsed from
     * cfg.para) at the CVodeSStolerances call site below. Same pattern
     * as the SHUD_CVODE_EPSLIN hook below — strtod parse with strict
     * range gate + fatal exit on malformed input + stderr provenance
     * line. env unset → reltol_effective = MD->CS.reltol, i.e. the
     * cfg.para value is used unchanged. */
    double reltol_effective = MD->CS.reltol;
    {
        const char *env_rel = getenv("SHUD_CVODE_RELTOL");
        if (env_rel != NULL && env_rel[0] != '\0') {
            char *endp = NULL;
            double rel = strtod(env_rel, &endp);
            if (endp == env_rel || rel <= 0.0 || rel >= 1.0) {
                fprintf(stderr,
                    "[shud] FATAL: SHUD_CVODE_RELTOL=%s is not a valid float in (0, 1)\n",
                    env_rel);
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
            reltol_effective = rel;
            fprintf(stderr,
                "[shud] PR-Z1 hook: CVODE reltol overridden %.4e -> %.4e\n",
                MD->CS.reltol, rel);
        }
        /* env unset → uses MD->CS.reltol from cfg.para. */
    }

    flag = CVodeSStolerances(cvode_mem, reltol_effective, MD->CS.abstol);
    check_flag(&flag, "CVodeSStolerances", 1);

    /* Factory dispatch. The default path (LINSOL_SPGMR) invokes
     * create_spgmr_ls; the AMG path is taken only when
     * SHUD_LINSOL=amg. */
    LS = (sel.sel == LINSOL_AMG)
       ? create_amg_ls(udata, MD, sunctx)
       : create_spgmr_ls(udata, sunctx);

    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    check_flag(&flag, "CVSpilsSetLinearSolver", 1);

    /* Optional CVODE linear convergence safety factor override, for
     * investigating nonlinear convergence failures (e.g. on the
     * experimental AMG path). When SHUD_CVODE_EPSLIN is set to a
     * value in (0, 1), it overrides the SUNDIALS default (0.05) via
     * CVodeSetEpsLin. Same style as the SHUD_SPGMR_MAXL hook above —
     * strtod parse with strict range gate + fatal exit on malformed
     * input + stderr provenance line. env unset → CVodeSetEpsLin is
     * not called, so CVODE's default 0.05 applies. */
    {
        const char *env_eps = getenv("SHUD_CVODE_EPSLIN");
        if (env_eps != NULL && env_eps[0] != '\0') {
            char *endp = NULL;
            double eps_lin = strtod(env_eps, &endp);
            if (endp == env_eps || eps_lin <= 0.0 || eps_lin >= 1.0) {
                fprintf(stderr,
                    "[shud] FATAL: SHUD_CVODE_EPSLIN=%s is not a valid float in (0, 1)\n",
                    env_eps);
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
            int rc = CVodeSetEpsLin(cvode_mem, eps_lin);
            if (rc != CV_SUCCESS) {
                fprintf(stderr,
                    "[shud] FATAL: CVodeSetEpsLin(%.4f) returned %d\n",
                    eps_lin, rc);
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,
                "[shud] PR-X1 hook: CVodeSetEpsLin(%.4f) applied\n",
                eps_lin);
        }
        /* env unset → CVODE's default (0.05) stays in effect. */
    }

    flag = CVodeSetMinStep(cvode_mem, 1E-6); //Minimum time interval in cvode.dt = t(i) - t(i - 1);
    check_flag(&flag, "CVodeSetMinStep", 1);
    
    flag = CVodeSetMaxNumSteps(cvode_mem, 1E6); //max iterations.
    check_flag(&flag, "CVodeSetMaxNumSteps", 1);
    
    flag = CVodeSetInitStep(cvode_mem, MD->CS.InitStep);
    check_flag(&flag, "CVodeSetInitStep", 1);
    
    //force cvode run at least every x time - units.t(i) - t(i - 1) < X;
    flag = CVodeSetMaxStep(cvode_mem, MD->CS.MaxStep);
    check_flag(&flag, "CVodeSetMaxStep", 1);
    
//    flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
    //flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
}
