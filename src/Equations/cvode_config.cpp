
#include "cvode_config.hpp"
/* S5c-B (#174): RHS 7-bucket + forcing I/O wall-clock timer dumps.
 * Header is empty under default (SHUD_ENABLE_DIAGNOSTICS undefined). */
#include "../Model/MD_diagnostics.hpp"

#include <errno.h>   /* errno (SHUD_SPGMR_MAXL strtol parse, p8tune-spgmr-maxl PR-C) */
#include <stdlib.h>  /* getenv, exit, EXIT_FAILURE (p8tune-g0 PR-0 SHUD_LINSOL) */
#include <string.h>  /* strcmp (p8tune-g0 PR-0 SHUD_LINSOL) */

/* P8-tune.G0 PR-0 (openspec change p8tune-g0-instrumented-amg-smoke).
 * BoomerAMG wrapper at SHUD/src/Equations/sunlinsol_hypre.{h,cpp}. */
#include "sunlinsol_hypre.h"

/* P8-tune.G0 PR-0 — runtime linear-solver selector enum + factory
 * dispatch. `SHUD_LINSOL` env var picks the inner solver passed to
 * CVODE; default-compat (`unset` OR `=spgmr`) replicates the pre-G0
 * SUNLinSol_SPGMR(udata, PREC_NONE, get_spgmr_maxl_from_env(), sunctx)
 * call site EXACTLY (G0-1 bit-identical anchor). `=amg` opt-in
 * dispatches the SUNLinSol_Hypre wrapper instead. Any other value
 * fatal-exits BEFORE `CVodeCreate`, so no CVODE state is allocated
 * on the error path (auditable: no leak). See design.md §D3. */
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

    /* Unrecognized value — fatal-exit BEFORE CVodeCreate. Spec
     * REQ-G0 fatal-exit-on-unknown invariant: stderr message names
     * the offender + the accepted set, exits non-zero, leaves no
     * CVODE state allocated. */
    fprintf(stderr,
            "[shud] FATAL: SHUD_LINSOL=%s unrecognized; accepted: spgmr, amg\n",
            env);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

/* Factory: SPGMR backend (PRE-G0 path — must produce
 * bit-identical output bytes vs the pre-G0 baseline). Encapsulates
 * the existing SUNLinSol_SPGMR call at the L324 site below; the
 * argument list is preserved byte-for-byte (PREC_NONE +
 * get_spgmr_maxl_from_env() + sunctx). `y` is the same N_Vector
 * the caller previously passed as `udata` — variable renaming
 * does not affect the linker. */
static SUNLinearSolver create_spgmr_ls(N_Vector y, SUNContext sunctx);

/* Factory: AMG backend via SUNLinSol_Hypre wrapper. Probes the
 * HYPRE runtime via the wrapper constructor (which returns NULL +
 * stderr on missing dylib / runtime init failure). On NULL return
 * the caller fatal-exits per spec REQ-G0 "fatal-exit on AMG
 * factory failure". */
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

    /* S0-8a / openMP #10 — optional key=value persistence. stdout output
     * above is unchanged so existing log-scraping continues to work.
     * Field order is deterministic + machine-parseable; the six fields
     * required by the b0-archive spec (nfe, nfeLS, nni, nli, nsetups,
     * netf) lead the file, the remainder follow for completeness. */
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
        /* S5c-A (#173) — S5c diagnostic channel additions (master plan
         * §S5c L1365 + spec s5c-solver-diagnostics "接入 SUNDIALS CVODE
         * stats 7 个 API"). The 5 existing keys above (nst / nfe / netf
         * / nni / nli) plus these 2 satisfy the spec 7-key contract.
         *
         * `hlast` / `qlast` are gated behind SHUD_ENABLE_DIAGNOSTICS so
         * the default build emits the same 15-key snapshot PR-12 froze
         * (B1a-tag bitwise invariant). Both are SUNDIALS 6.0.0 public
         * API reads (post-solve, no RHS path mutation) — diagnostics-ON
         * dat outputs remain bitwise == B1a-tag; only this file's
         * trailing key set differs. */
        int qlast;
        realtype hlast;
        flag = CVodeGetLastStep(cvode_mem, &hlast);
        check_flag(&flag, "CVodeGetLastStep", 1);
        flag = CVodeGetLastOrder(cvode_mem, &qlast);
        check_flag(&flag, "CVodeGetLastOrder", 1);
        fprintf(fout, "hlast=%.17g\n", (double)hlast);
        fprintf(fout, "qlast=%d\n",    qlast);

        /* S5c-B (#174): RHS 7-bucket + forcing I/O wall-clock dump.
         * Sum of pct_rhs_* SHALL ∈ [99.5%, 100.5%] per spec scenario
         * "7 个 bucket 输出完整时间分布". Buckets 0-6 are defined in
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


/* p8tune-spgmr-maxl PR-C (capability spgmr-maxl-env-hook, GitHub #366).
 *
 * Runtime env-var hook for the SPGMR Krylov subspace dimension `maxl`
 * passed to SUNLinSol_SPGMR at L259. Default-unset / "" / "0" / "5" are
 * bit-identical to the prior SHUD 37be0fe production behavior (SUNDIALS
 * docs: maxl <= 0 collapses to the documented default 5; explicit 5 is
 * the default). Opt-in values {10, 15, 20, 30} flow through to SUNDIALS
 * unchanged + emit a stdout provenance line consumed by the PR-D 60-cell
 * sweep aggregator. Any other value (e.g. "7", "50", "foo", "-1") aborts
 * via myexit(ERRCVODE) BEFORE any SPGMR allocation so an invalid sbatch
 * cell fails fast rather than silently overwriting a baseline artifact.
 *
 * See openspec/changes/p8tune-spgmr-maxl/specs/spgmr-maxl-env-hook/spec.md
 * and design.md §D15 Invariant Matrix for the full contract. */
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
     * -> default 5). Preserve silent-default bit-identical contract per
     * design D15 Invariant Matrix regression rows: unset / "" / "0" all
     * suppress the provenance log line. */
    if (val == 0) {
        return 0;
    }

    /* val in {5, 10, 15, 20, 30}: emit provenance line so the PR-D
     * aggregator can attribute each cell to its maxl value. */
    fprintf(stdout, "[CVODE] SPGMR maxl=%ld pretype=PREC_NONE\n", val);
    fflush(stdout);
    return (int)val;
}

static SUNLinearSolver create_spgmr_ls(N_Vector y, SUNContext sunctx)
{
    /* Mirror the pre-G0 hardcoded SUNLinSol_SPGMR call at the
     * previous L324 site EXACTLY (PREC_NONE + maxl from env hook +
     * sunctx). The variable rename `udata` -> `y` is local to this
     * function and does not change generated code: SUNLinSol_SPGMR
     * sees the same N_Vector pointer the caller previously passed. */
    SUNLinearSolver LS = SUNLinSol_SPGMR(y, PREC_NONE,
                                         get_spgmr_maxl_from_env(),
                                         sunctx);
    check_flag((void *)LS, "SUNLinSol_SPGMR", 0);
    return LS;
}

static SUNLinearSolver create_amg_ls(N_Vector y, Model_Data *MD,
                                     SUNContext sunctx)
{
    /* G0: hardcoded (interp_type=6, coarsen_type=8). The wrapper
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

    /* P8-tune.G0 PR-0 — linsol selection parsed BEFORE CVodeCreate
     * so the fatal-exit-on-unknown invariant holds (no CVODE state
     * leaked on the error path). Marker line precedes any CVODE
     * output for log post-mortem ordering. */
    const linsol_selection sel = parse_linsol_env();
    fprintf(stdout, "[shud] linsol=%s source=%s\n",
            sel.sel_name, sel.source);
    fflush(stdout);

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

    flag = CVodeSStolerances(cvode_mem, MD->CS.reltol, MD->CS.abstol);
    check_flag(&flag, "CVodeSStolerances", 1);

    /* Factory dispatch (P8-tune.G0 PR-0). Default path
     * (LINSOL_SPGMR) invokes create_spgmr_ls which calls
     * SUNLinSol_SPGMR with byte-identical args to the pre-G0 site
     * — the G0-1 bit-identical anchor depends on this. */
    LS = (sel.sel == LINSOL_AMG)
       ? create_amg_ls(udata, MD, sunctx)
       : create_spgmr_ls(udata, sunctx);

    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    check_flag(&flag, "CVSpilsSetLinearSolver", 1);

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
