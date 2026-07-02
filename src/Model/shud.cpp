#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include "f_element.hpp"
//#include "f_River.hpp"
#include "f.hpp"
#include "IO.hpp"
#include "ModelConfigure.hpp"
#include "print.hpp"
#include "Macros.hpp"
#include "functions.hpp"
//#include "is_sm_et.hpp"
#include "cvode_config.hpp"
#include "Model_Data.hpp"
#include "TimeSeriesData.hpp"
#include "FloodAlert.hpp"
#include "CommandIn.hpp"
/* P11-osc PR-D1 (#434) — env-gated (strict `=1`) CVODE stepping diagnostics.
 * Default-off: when neither SHUD_DIAG_DT=1 nor SHUD_DIAG_OSC=1 is set, every
 * hook below is a no-op and the default hot path is byte-identical (keliya
 * B0 SHA gate). Reads ONLY the accepted CVODE state vector — never the uY*
 * RHS scratch globals (source-audit spec). */
#include "MD_osc_diag.hpp"
/* P12-nvec PR-N0 (#442) — env-gated (strict `=1`) NVector op-share
 * profiler. Default-off: when SHUD_NVEC_PROF!=1 the ops-table shims are
 * NEVER installed and the coupled-vector creation path is byte-identical
 * (keliya B0 SHA gate). When on, shims are pure delegation (counting +
 * monotonic-ns), so trajectories stay identical; nvec_prof.csv is dumped
 * at run end. Wraps ONLY the coupled udata/du family; the decoupled
 * 5-solver vectors (SHUD_uncouple below) are never wrapped. Composition
 * order: the profiler wraps LAST/OUTERMOST so a future PR-N1 hybrid
 * override (installed on the same vector BEFORE this call) is delegated
 * through. No code here lives on the RHS f() path (f.cpp / MD_rhs_core.cpp
 * untouched). */
#include "MD_nvec_prof.hpp"
/* P8-tune.G0 PR-B (#411) + PR-0 Phase 7 finding #4 cleanup —
 * SUNLinSol_Hypre_DrainTelemetry + SUNLinSolFree on shutdown. The
 * wrapper is link-always (cvode_config.cpp dispatches via SHUD_LINSOL);
 * the drain is a per-Solve telemetry ring drained to TSV at process
 * shutdown when $SHUD_TELEMETRY_TSV is set. SUNLinSolFree releases the
 * wrapper's HypreContent (Hypre handles + ring buffer) — addresses
 * PR-0 #414 Phase 7 finding #4 process-exit cleanup leak. On
 * SHUD_LINSOL=spgmr the drain is a no-op (wrapper not in use; the
 * SPGMR LS is a stock SUNDIALS object and SUNLinSolFree releases that
 * cleanly too — the wrapper's drain probes LS->ops->getid and returns
 * 0 if the handle is not the custom AMG wrapper). */
#include "sunlinsol_hypre.h"

/* S0-8a / openMP #10 — wall-clock profile timer infrastructure. Header
 * lives in the outer `tools/profile/` directory, sibling to SHUD/. We
 * #ifdef-guard the include so the SHUD submodule on `openmp-baseline`
 * stays self-contained at PROFILE=0 (the only consumer of the timer
 * API — `shud_profile::dump` below — is itself #ifdef-guarded, so the
 * no-op stubs in timer.h are dead weight at PROFILE=0). PROFILE=1
 * builds get the include + -I + impl source via Makefile injection. */
#ifdef SHUD_ENABLE_PROFILE
#include "timer.h"
#endif

double *uYsf;
double *uYus;
double *uYgw;
double *uYriv;
double *uYlake;
double *globalY;
double timeNow;
int dummy_mode = 0;
int global_fflush_mode = 0;
int global_implicit_mode = 1;
int global_verbose_mode = 1;
int lakeon = 0; /* Whether lake module ON(1), OFF(0) */
/* S5d.3 (#181) — NUMA first-touch gate. 1 = OMP_PROC_BIND was set at
 * startup so threads have a deterministic core affinity and Linux
 * first-touch page policy can route SoA pages to the consumer thread's
 * local NUMA node; 0 = OMP_PROC_BIND unset (design R3 mitigation #2:
 * skip parallel first-touch entirely to keep behaviour deterministic
 * and bitwise-identical to the serial baseline). Set ONCE at SHUD()
 * entry from getenv("OMP_PROC_BIND") and read by malloc_EleRiv() /
 * LoadIC() before each parallel first-touch loop. */
int g_numa_first_touch_enabled = 0;
using namespace std;

/* S5d.3 (#181) — deterministic NUMA log-token emitter, called once at
 * the start of every SHUD() entry point. Spec L75-77 grep ordering
 * assertion requires "[NUMA] OMP_PROC_BIND=" to appear BEFORE any
 * "[NUMA] first-touch begin" line — that ordering is satisfied by the
 * call sequence in SHUD()/SHUD_uncouple(): emit_numa_token() ->
 * MD->initialize() (which calls malloc_EleRiv with 3 first-touch
 * sites) -> MD->LoadIC() (1 first-touch site).
 *
 * S5d.4 (#182) — additional stderr WARNING emitted alongside the
 * stdout [NUMA] OMP_PROC_BIND=unset token, satisfying spec
 * s5d-data-layout-soa-numa Scenario "shud.cpp warning 路径生效"
 * (Requirement "线程绑定 run script 与 manifest 字段必填"). The
 * stdout token is preserved unchanged (PR #181 grep ordering gate
 * keeps consuming it from the run log); the new stderr warning is
 * the operator-facing channel matching design D5: program SHALL NOT
 * override OMP_PROC_BIND, only surface its absence as a warning that
 * points at tools/run_omp.sh as the fix. The wording embeds the
 * exact spec phrase so the spec grep assertion is verbatim. */
static void emit_numa_token(void){
    const char *bind = getenv("OMP_PROC_BIND");
    if (bind != NULL && bind[0] != '\0') {
        printf("[NUMA] OMP_PROC_BIND=%s\n", bind);
        g_numa_first_touch_enabled = 1;
    } else {
        printf("[NUMA] OMP_PROC_BIND=unset\n");
        printf("[NUMA] WARNING: OMP_PROC_BIND unset - skipping first-touch optimization for determinism guarantee.\n");
        /* S5d.4 (#182): stderr WARNING — operator-facing channel.
         * Spec phrase verbatim: "OMP_PROC_BIND not set, NUMA
         * first-touch may be ineffective." Append a pointer to the
         * canonical fix (tools/run_omp.sh) so the user gets a
         * one-step recovery without consulting docs. */
        fprintf(stderr,
                "[OMP] WARNING: OMP_PROC_BIND not set, NUMA "
                "first-touch may be ineffective. Use "
                "tools/run_omp.sh to set defaults.\n");
        fflush(stderr);
        g_numa_first_touch_enabled = 0;
    }
    fflush(stdout);
}

double SHUD(FileIn *fin, FileOut *fout){
    double ret = 0.;
    Model_Data  *MD;        /* Model Data                */
    N_Vector    udata;
    N_Vector    du;

    /* S5d.3 (#181) — emit deterministic NUMA log token + set
     * g_numa_first_touch_enabled BEFORE any malloc_EleRiv / LoadIC
     * call so spec L75-77 grep ordering assertion holds:
     *   min_line('[NUMA] OMP_PROC_BIND=') < min_line('[NUMA] first-touch begin')
     */
    emit_numa_token();

    SUNContext sunctx;
    ret = SUNContext_Create(NULL, &sunctx);
    check_flag(&ret, "SUNContext_Create", 1);

    void    *mem = NULL;
    SUNLinearSolver LS = NULL;
    int     flag;            /* flag to test return value */
    double  t, tnext;    /* stress period & step size */
    int NY = 0;
    int ierr = 0;
    /* allocate memory for model data structure */
    MD = new Model_Data(fin, fout);
    MD->loadinput();
    MD->initialize();
    MD->CheckInputData();
    fout->updateFilePath();
    NY = MD->NumY;
    globalY = new double[NY];
    /* S1d.2 (openMP #48) — N_Vector backend dispatch is now keyed on
     * SHUD_USE_OPENMP_NVECTOR (renamed from the legacy
     * three-concerns-conflated switch). Defaults to OFF → Serial
     * backend, which is the validation surface for Config A
     * (bitwise vs B0). When SHUD_USE_OPENMP_NVECTOR=1 the build
     * additionally links libsundials_nvecopenmp + pulls in
     * nvector_openmp.h via Macros.hpp. */
#ifdef SHUD_USE_OPENMP_NVECTOR
    omp_set_num_threads(MD->CS.num_threads);
    screeninfo("\nopenMP NVector: ON. No of Threads = %d\n", MD->CS.num_threads);
    udata = N_VNew_OpenMP(NY, MD->CS.num_threads, sunctx);
    du = N_VNew_OpenMP(NY, MD->CS.num_threads, sunctx);
#else
    screeninfo("\nopenMP NVector: OFF (Serial backend)\n");
    udata = N_VNew_Serial(NY, sunctx);
    du = N_VNew_Serial(NY, sunctx);
#endif
    /* P1e PR-G (#315 / design D3) — StrictOMP RHS startup single-point
     * thread-count setup. Independent of SHUD_USE_OPENMP_NVECTOR (Config
     * B's `omp_set_num_threads(MD->CS.num_threads)` above stays untouched
     * to preserve mode B/D historical NVector parity). Config C (Serial
     * NVec + StrictOMP RHS) has no other omp_set_num_threads call site
     * — without this block, `#pragma omp parallel` in MD_rhs_core.cpp
     * StrictOMP case would default to `omp_get_max_threads()` (TE-1
     * defect class). `SHUD_RHS_THREADS` env wins; if unset, fall back
     * to whatever `omp_get_max_threads()` reports (which is itself
     * driven by `OMP_NUM_THREADS` when set, else the OpenMP runtime
     * default). In Config D this re-set overrides the Config-B-style
     * MD->CS.num_threads set above, which is the intended behaviour:
     * the user wants `SHUD_RHS_THREADS` to be the canonical knob across
     * Config C/D, while NVector allocation in Config D has already
     * fixed the NVector-backend thread count at `N_VNew_OpenMP` time. */
#if defined(SHUD_ENABLE_OPENMP_RHS)
    {
        const char* shud_rhs_threads_env = getenv("SHUD_RHS_THREADS");
        int rhs_threads = shud_rhs_threads_env ? atoi(shud_rhs_threads_env) : 0;
        if (rhs_threads <= 0) {
            rhs_threads = omp_get_max_threads();
        }
        omp_set_num_threads(rhs_threads);
        fprintf(stdout, "P1e startup: SHUD_RHS_THREADS=%s -> omp_set_num_threads(%d); omp_get_max_threads=%d\n",
                shud_rhs_threads_env ? shud_rhs_threads_env : "(unset)",
                rhs_threads,
                omp_get_max_threads());
    }
#endif
    screeninfo("\nGlobal Implicit Mode: ON\n");
    MD->LoadIC();
    MD->SetIC2Y(udata);
    MD->initialize_output();
    MD->PrintInit(fout->Init_bak, 0);
    MD->InitFloodAlert(fout->floodout);
    /* P12-nvec PR-N0 (#442) — install the NVector op-share profiler shims
     * on the coupled vectors HERE: after creation (and after any future
     * SHUD-owned ops-table override, which a PR-N1 hybrid build inserts
     * above), and BEFORE SetCVODE (which calls CVodeInit → N_VClone
     * allocates every internal temporary, propagating the shim table).
     * No-op unless SHUD_NVEC_PROF=1. Wrap both udata and du (same family;
     * install is idempotent). The clone-propagation smoke assert runs once
     * on udata so the stdout PASS line is captured in the evidence log. */
    {
#ifdef SHUD_USE_OPENMP_NVECTOR
        const char *nvec_prof_backend = "openmp";
#else
        const char *nvec_prof_backend = "serial";
#endif
        nvec_prof_install(udata, nvec_prof_backend);
        nvec_prof_install(du, nvec_prof_backend);
        if (nvec_prof_is_on()) {
            nvec_prof_clone_carries_shims(udata);
        }
    }
    SetCVODE(mem, f, MD, udata, LS, sunctx);
    /* set start time */
    t = MD->CS.StartTime;
    tnext = t;
    //CheckInput(MD, &CS);
    /* start solver in loops */
//    getSecond();
    MD->modelSummary(0);
    MD->debugData(fout->outpath);
    MD->gc.write(fout->Calib_bak);
//    f(t, udata, du, MD); /* Initialized the status */
    /* P11-osc PR-D1 (#434) — construct + prime the env-gated diagnostics
     * BEFORE the profiled solver-loop scope so header-write / buffer alloc /
     * initial-state snapshot do not skew t_wall_total (mirrors the existing
     * "init lives outside this scope" intent). No-op unless SHUD_DIAG_DT=1
     * or SHUD_DIAG_OSC=1 (strict `=1`). State read via N_VGetArrayPointer
     * (accepted CVODE state), never the uY* RHS scratch globals. */
    OscDiag diag;
    if (diag.any_on()) {
        diag.begin(mem, udata, fout->projectname, MD->CS.SolverStep,
                   MD->NumEle, MD->NumRiv, fout->outpath);
    }
    {
#ifdef SHUD_ENABLE_PROFILE
        /* S0-10 / openMP #14 — t_wall_total wraps the main solver loop
         * (NumSteps iterations, each with forcing/ET/CVode/summary/
         * ExportResults). Used in dump() to derive t_other = wall_total
         * - (CVODE_raw + forcing + ET + output). Initialization /
         * cvode_stats persistence / profile dump itself live outside
         * this scope on purpose so they do not skew the loop wall. */
        shud_profile::Timer _t_wall("t_wall_total");
#endif
        for (int i = 0; i < MD->CS.NumSteps && !ierr; i++) {
            printDY(MD->file_debug);
#ifdef DEBUG
            printDY(MD->file_debug);
#endif
            flag = MD->ScreenPrint(t, i);
            MD->PrintInit(fout->Init_update, t);
            /* inner loops to next output points with ET step size control */
            tnext += MD->CS.SolverStep;
            while (t < tnext) {
                {
#ifdef SHUD_ENABLE_PROFILE
                    /* P2a: wrap forcing-file disk I/O + interp. */
                    shud_profile::Timer _t_fr("t_forcing_io");
#endif
                    MD->updateforcing(t);
                }
                /* calculate Interception Storage */
                {
#ifdef SHUD_ENABLE_PROFILE
                    /* P2a: ET / canopy / snow physics bucket. */
                    shud_profile::Timer _t_et("t_ET");
#endif
                    MD->ET(t, tnext);
                }
                if(dummy_mode){
                    t = tnext;  /* dummy mode only. */
                }else{
#ifdef SHUD_ENABLE_PROFILE
                    /* t_CVODE_raw includes the RHS sub-calls (CVode
                     * invokes f() internally). dump() subtracts the
                     * already-measured t_RHS_total to get the net
                     * t_CVODE_internal bucket. */
                    shud_profile::Timer _t_cvode("t_CVODE_raw");
#endif
                    flag = CVode(mem, tnext, udata, &t, CV_NORMAL);
                    check_flag(&flag, "CVode", 1);
                }
            }
            /* P11-osc PR-D1 (#434): per-interval diagnostic sample at the
             * accepted CVode-return boundary (t == tnext here; the inner
             * while runs exactly once per SolverStep in CV_NORMAL). Emits
             * one diag_dt_trace.csv row (counter deltas) and folds this
             * interval into the flip counters. Placed AFTER the while so it
             * is outside the t_CVODE_raw profile scope. No-op unless a
             * strict `=1` gate is set. */
            if (diag.any_on()) {
                diag.record(mem, udata, t);
            }
            //            CVODEstatus(mem, udata, t);
            {
#ifdef SHUD_ENABLE_PROFILE
                /* P2a: summary() flux post-proc to t_output bucket. */
                shud_profile::Timer _t_out("t_output");
#endif
                MD->summary(udata);
            }
            /* P1e PR-B (#310): SHUD_DUMP_CV_Y=1 env gates per-tout CV_Y state
             * vector dump. Used by tools/p1e_cv_y_hash to verify cross-build
             * (mode A/B/C/D) solver-state byte-equality at tout boundary.
             * Gated by env var (runtime-toggleable, no recompile required).
             * Disabled by default → no behavior change for normal builds. */
            if (getenv("SHUD_DUMP_CV_Y") != NULL) {
                double *cv_y_dump = N_VGetArrayPointer(udata);
                char cv_y_path[MAXLEN];
                snprintf(cv_y_path, sizeof(cv_y_path),
                         "%s/cv_y_%020.6f.bin", fout->outpath, t);
                FILE *cv_y_fp = fopen(cv_y_path, "wb");
                if (cv_y_fp != NULL) {
                    size_t nwritten = fwrite(cv_y_dump, sizeof(double),
                                             (size_t)NY, cv_y_fp);
                    fclose(cv_y_fp);
                    if (nwritten != (size_t)NY) {
                        fprintf(stderr,
                                "[CV_Y_DUMP] WARNING: short write %zu/%d "
                                "at t=%f for %s\n",
                                nwritten, NY, t, cv_y_path);
                    }
                } else {
                    fprintf(stderr,
                            "[CV_Y_DUMP] WARNING: cannot open %s for write "
                            "(errno preserved by fopen)\n", cv_y_path);
                }
            }
            /* P1e PR-B0 (#323): recompute river/lake/element flux caches
             * from Y(tout) before ExportResults fires PrintData. Fixes
             * non-determinism caused by PCtrl reading the side-effect
             * cache left by CVODE's last internal-step f() at
             * t_internal != tout (per docs/p1e/p1e_rivqdown_cache_audit.md
             * + spec p1e-strict-omp-rhs L260-285 + design D5 option 1). */
            MD->recompute_for_output(udata, t);
            {
#ifdef SHUD_ENABLE_PROFILE
                /* P2a: ExportResults disk-write to t_output (shared bucket). */
                shud_profile::Timer _t_out2("t_output");
#endif
                MD->CS.ExportResults(t);
            }
            MD->flood->FloodWarning(t);
        }
    }
    /* P11-osc PR-D1 (#434): dump the flip-counter CSVs at run end (no-op
     * unless SHUD_DIAG_OSC=1). Outside the t_wall_total scope by design. */
    if (diag.any_on()) {
        diag.finish();
    }
    /* P12-nvec PR-N0 (#442): dump nvec_prof.csv at run end (no-op unless
     * SHUD_NVEC_PROF=1). Reads only the per-op global counters, not the
     * vectors, so it is safe to call before/after N_VDestroy; placed here
     * (before free) to mirror the diag.finish() run-end pattern and stay
     * outside any profile-timer scope. */
    if (nvec_prof_is_on()) {
#ifdef SHUD_USE_OPENMP_NVECTOR
        const char *nvec_prof_backend = "openmp";
#else
        const char *nvec_prof_backend = "serial";
#endif
        /* Report the EFFECTIVE OpenMP thread count driving the run (== the
         * N knob: OMP_NUM_THREADS / SHUD_RHS_THREADS via omp_get_max_threads),
         * NOT MD->CS.num_threads (the cfg NUM_OPENMP value, which in Config C
         * only governs the unused Serial-NVector thread hint and can differ
         * from the run's N). Falls back to CS.num_threads when OpenMP is not
         * compiled in (plain serial `make shud`). */
#if defined(_OPENMP)
        int nvec_prof_nthreads = omp_get_max_threads();
#else
        int nvec_prof_nthreads = MD->CS.num_threads;
#endif
        nvec_prof_dump(fout->projectname, NY, nvec_prof_nthreads,
                       nvec_prof_backend, fout->outpath);
    }
    MD->ScreenPrint(t, MD->CS.NumSteps);
    MD->PrintInit(fout->Init_update, t);
    MD->modelSummary(1);
    /* Free memory.
     * S1d.2 (openMP #48) — the prior type-specific Serial destroy was
     * unsafe under SHUD_USE_OPENMP_NVECTOR=1 (it would receive an
     * N_VNew_OpenMP-allocated vector with a different content layout
     * and trigger UB; master plan §4.19). The generic `N_VDestroy`
     * dispatches via the N_Vector ops table and correctly routes to
     * whichever backend created `v`, so the same call works for both
     * Serial and OpenMP backends. */
    N_VDestroy(udata);
    N_VDestroy(du);

    /* S0-8a / openMP #10 — persist CVODE final stats next to the SHUD
     * output dir for the B0 archive script to pick up. stdout printout
     * inside PrintFinalStats is unchanged (back-compat). fopen failure
     * is non-fatal: PrintFinalStats(mem, NULL) still prints to stdout. */
    {
        char stats_path[MAXLEN];
        snprintf(stats_path, sizeof(stats_path), "%s/cvode_stats.txt",
                 fout->outpath);
        FILE *stats_fp = fopen(stats_path, "w");
        PrintFinalStats(mem, stats_fp);
        if (stats_fp != NULL) {
            fclose(stats_fp);
        } else {
            fprintf(stderr,
                    "[shud] WARN: cvode_stats.txt fopen failed at "
                    "'%s'; stdout-only fallback used.\n",
                    stats_path);
        }
    }

    /* S5c-C (#175): nFCall lives in its own file, NOT in cvode_stats.txt.
     * Per spec scenario "nFCall 独立 channel 上报" + "15-key snapshot 不包含
     * nFCall". This decouples the free-running SHUD counter from the 15-key
     * invariance gate. The file is small (1-2 lines) and read by the CI
     * workflow's nFCall column step + tools/cvode_stats_diff post-checks. */
    {
        char nfcall_path[MAXLEN];
        snprintf(nfcall_path, sizeof(nfcall_path), "%s/nfcall.txt",
                 fout->outpath);
        FILE *nfcall_fp = fopen(nfcall_path, "w");
        if (nfcall_fp != NULL) {
            fprintf(nfcall_fp, "nFCall=%lu\n", MD->nFCall);
            fclose(nfcall_fp);
        } else {
            fprintf(stderr,
                    "[shud] WARN: nfcall.txt fopen failed at '%s'; "
                    "stdout-only fallback used.\n", nfcall_path);
        }
    }

    /* P8-tune.G0 PR-B (#411): drain wrapper telemetry ring buffer to
     * $SHUD_TELEMETRY_TSV before freeing the LS. Env-var unset / drain
     * file fopen fail / SPGMR (non-AMG wrapper) — all benign no-ops:
     *   - getenv NULL ⇒ skip drain entirely.
     *   - fopen NULL  ⇒ emit one warn line, continue (do not abort).
     *   - SPGMR LS    ⇒ DrainTelemetry probes getid() and returns 0.
     * The wrapper drain writes TSV header + rows, then resets the ring
     * buffer head/tail (idempotent across re-invocations).
     *
     * PR-B #416 Phase 6 P1-6: drain BEFORE CVodeFree. CVODE 6.0
     * cvLsFree (CVodeFree teardown path) may touch LS->content during
     * integrator-side cleanup; draining first guarantees the wrapper's
     * HypreContent ring buffer is intact when DrainTelemetry walks it.
     * The reverse order (drain after CVodeFree) would also work for
     * SPGMR (no LS->content state to lose) but is unsafe for the AMG
     * wrapper whose HypreContent owns the telemetry buffer. */
    {
        const char *telemetry_tsv = getenv("SHUD_TELEMETRY_TSV");
        if (telemetry_tsv != NULL && telemetry_tsv[0] != '\0' && LS != NULL) {
            FILE *tsv_fp = fopen(telemetry_tsv, "w");
            if (tsv_fp != NULL) {
                int written = SUNLinSol_Hypre_DrainTelemetry(LS, tsv_fp);
                fclose(tsv_fp);
                fprintf(stdout,
                        "[shud-G0] SUNLinSol_Hypre_DrainTelemetry: "
                        "wrote %d telemetry rows to %s\n",
                        written, telemetry_tsv);
                fflush(stdout);
            } else {
                fprintf(stderr,
                        "[shud-G0] WARN: SHUD_TELEMETRY_TSV fopen failed "
                        "at '%s'; telemetry not drained (run continues).\n",
                        telemetry_tsv);
            }
        }
    }

    /* Free integrator memory */
    CVodeFree(&mem);

    /* P8-tune.G0 PR-B (#411) + PR-0 #414 Phase 7 finding #4 — release
     * the SUNLinearSolver (wrapper releases its HypreContent which holds
     * the Hypre AMG handle + IJ matrix/vector handles + ring buffer;
     * stock SPGMR LS releases its workspace too). Guarded against the
     * pre-G0 NULL-default LS case (the wrapper or SPGMR ctor populates
     * LS via SetCVODE; if SetCVODE bailed early LS remains NULL). */
    if (LS != NULL) {
        SUNLinSolFree(LS);
    }

#ifdef SHUD_ENABLE_PROFILE
    /* S0-8a / openMP #10 — dump profile bucket skeleton. #10 ships
     * infrastructure only; bucket values are all 0.0 until S0-10
     * adds the actual instrumentation hook points. */
    {
        char prof_path[MAXLEN];
        snprintf(prof_path, sizeof(prof_path), "%s/profile_B0.yaml",
                 fout->outpath);
        shud_profile::dump(prof_path);
    }
#endif

    SUNContext_Free(&sunctx);
    delete MD;
    return ret;
}


double SHUD_uncouple(FileIn *fin, FileOut *fout){
    double ret = 0.;
    Model_Data  *MD;        /* Model Data                */
    N_Vector    u1, u2, u3, u4, u5;
    N_Vector    du1, du2, du3, du4, du5;

    /* S5d.3 (#181) — same NUMA token + gate emit as in SHUD() (kept
     * symmetric so uncouple-path runs also satisfy the L75-77 grep
     * ordering assertion). */
    emit_numa_token();

    SUNContext sunctx1, sunctx2, sunctx3, sunctx4, sunctx5;
    ret = SUNContext_Create(NULL, &sunctx1);check_flag(&ret, "SUNContext_Create", 1);
    ret = SUNContext_Create(NULL, &sunctx2);check_flag(&ret, "SUNContext_Create", 1);
    ret = SUNContext_Create(NULL, &sunctx3);check_flag(&ret, "SUNContext_Create", 1);
    ret = SUNContext_Create(NULL, &sunctx4);check_flag(&ret, "SUNContext_Create", 1);
    ret = SUNContext_Create(NULL, &sunctx5);check_flag(&ret, "SUNContext_Create", 1);
    
    void    *mem1 = NULL, *mem2 = NULL, *mem3 = NULL, *mem4 = NULL, *mem5 = NULL;
    SUNLinearSolver LS1 = NULL, LS2 = NULL, LS3 = NULL, LS4 = NULL, LS5 = NULL;
    int     flag;            /* flag to test return value */
    double  t = 0, dt = 0, tout = 0;    /* stress period & step size */
    int NY = 0;
    int N1, N2, N3, N4, N5;
    int ierr = 0;
    /* allocate memory for model data structure */
    MD = new Model_Data(fin, fout);
    MD->loadinput();
    MD->initialize();
    MD->CheckInputData();
    fout->updateFilePath();
    NY = MD->NumY;
    N1 = MD->NumEle;
    N2 = MD->NumEle;
    N3 = MD->NumEle;
    N4 = MD->NumRiv;
    N5 = MD->NumLake;

    screeninfo("\nopenMP: OFF\n");
    screeninfo("\nGlobal Implicit Mode: OFF\n");
    u1 = N_VNew_Serial(N1,sunctx1);
    u2 = N_VNew_Serial(N2,sunctx2);
    u3 = N_VNew_Serial(N3,sunctx3);
    u4 = N_VNew_Serial(N4,sunctx4);
    u5 = N_VNew_Serial(N5,sunctx5);
    
    du1 = N_VNew_Serial(N1,sunctx1);
    du2 = N_VNew_Serial(N2,sunctx2);
    du3 = N_VNew_Serial(N3,sunctx3);
    du4 = N_VNew_Serial(N4,sunctx4);
    du5 = N_VNew_Serial(N5,sunctx5);

    MD->LoadIC();
    MD->SetIC2Y(u1, u2, u3, u4, u5);
    MD->initialize_output();
    MD->PrintInit(fout->Init_bak, 0);
    MD->InitFloodAlert(fout->floodout);
    
    SetCVODE(mem1, f_surf,  MD, u1, LS1, sunctx1);
    SetCVODE(mem2, f_unsat, MD, u2, LS2, sunctx2);
    SetCVODE(mem3, f_gw,    MD, u3, LS3, sunctx3);
    SetCVODE(mem4, f_river, MD, u4, LS4, sunctx4);
    SetCVODE(mem5, f_lake,  MD, u5, LS5, sunctx5);
    
//    flag = CVodeSetMaxStep(mem1, max(MD->CS.MaxStep/4., 1.) );
//    check_flag(&flag, "CVodeSetMaxStep", 1);
    
    /* set start time */
    t = MD->CS.StartTime;
    double tnext = t;
    //CheckInput(MD, &CS);
    /* start solver in loops */
//    getSecond();
    MD->modelSummary(0);
    MD->debugData(fout->outpath);
    MD->gc.write(fout->Calib_bak);
    
//    FILE *fp1, *fp2, *fp3, *fp4;
//    fp1=fopen("y1.txt", "w");
//    fp2=fopen("y2.txt", "w");
//    fp3=fopen("y3.txt", "w");
//    fp4=fopen("y4.txt", "w");
    double t0 = t, tnext_et = tnext;
    for (int i = 0; i < MD->CS.NumSteps && !ierr; i++) {
        /* inner loops to next output points with ET step size control */
        tnext += MD->CS.SolverStep;
        while (t < tnext ) {
//            if (t + MD->CS.ETStep >=tnext) {
                tout = tnext;
//            } else {
//                tout = t + MD->CS.ETStep;
//            }
            dt = tout - t;
            MD->updateforcing(t);
//            if(t >= tnext_et){
                /* calculate Interception Storage */
                MD->ET(t, tnext);
//                tnext_et += MD->CS.ETStep;
//            }
            
            t=t0;
            MD->t0=t0; MD->t1=tout;
            Global2Sub(MD->NumEle, MD->NumRiv, MD->NumLake);
            flag = CVode(mem1, tout, u1, &t, CV_NORMAL);
            check_flag(&flag, "CVode1 SURF", 1);
            
            t=t0;
            Global2Sub(MD->NumEle, MD->NumRiv, MD->NumLake);
            flag = CVode(mem2, tout, u2, &t, CV_NORMAL);
            check_flag(&flag, "CVode2 UNSAT", 1);
            
            t=t0;
            Global2Sub(MD->NumEle, MD->NumRiv, MD->NumLake);
            flag = CVode(mem3, tout, u3, &t, CV_NORMAL);
            check_flag(&flag, "CVode3 GW", 1);
            
            t=t0;
            Global2Sub(MD->NumEle, MD->NumRiv, MD->NumLake);
            flag = CVode(mem4, tout, u4, &t, CV_NORMAL);
            check_flag(&flag, "CVode4 RIV", 1);
            
            if(lakeon && N5 > 0){
                t=t0;
                Global2Sub(MD->NumEle, MD->NumRiv, MD->NumLake);
                flag = CVode(mem5, tout, u2, &t, CV_NORMAL);
                check_flag(&flag, "CVode5 LAKE", 1);
            }
        }
        t0 = t;
        MD->summary(u1, u2, u3, u4, u5);
        MD->CS.ExportResults(t);
        flag = MD->ScreenPrintu(t, i);
        MD->PrintInit(fout->Init_update, t);
//        printVector(fp1, globalY, 0, N1, t);
//        printVector(fp2, globalY, N1, N2, t);
//        printVector(fp3, globalY, N1*2, N3, t);
//        printVector(fp4, globalY, N1*3, N4, t);
        MD->flood->FloodWarning(t);
    }
//    fclose(fp1);
//    fclose(fp2);
//    fclose(fp3);
//    fclose(fp4);
    MD->modelSummary(1);
    /* Free memory — generic N_VDestroy dispatch (see Coupled-path
     * comment above for the §4.19 backend-mismatch UB rationale).
     * Uncouple path currently only ever allocates Serial vectors
     * (N_VNew_Serial calls below) so the destroy was never strictly
     * broken; we migrate for consistency + to make future OMP-backend
     * uncouple support a zero-touch change. */
    N_VDestroy(u1);
    N_VDestroy(u2);
    N_VDestroy(u3);
    N_VDestroy(u4);
    N_VDestroy(u5);

    N_VDestroy(du1);
    N_VDestroy(du2);
    N_VDestroy(du3);
    N_VDestroy(du4);
    N_VDestroy(du5);

    /* S0-8a / openMP #10 — persist CVODE final stats from the surface
     * solver (mem1) as the representative. We pick mem1 because it is
     * the driving solver in the uncouple loop and its counters cover
     * the longest model-time span; mem2..mem5 stats remain on stdout
     * only (a future schema rev can split per-mem if needed). */
    {
        char stats_path[MAXLEN];
        snprintf(stats_path, sizeof(stats_path), "%s/cvode_stats.txt",
                 fout->outpath);
        FILE *stats_fp = fopen(stats_path, "w");
        PrintFinalStats(mem1, stats_fp);
        if (stats_fp != NULL) {
            fclose(stats_fp);
        } else {
            fprintf(stderr,
                    "[shud] WARN: cvode_stats.txt fopen failed at "
                    "'%s'; stdout-only fallback used.\n",
                    stats_path);
        }
    }

    /* S5c-C (#175): nFCall channel emission in the uncouple path. Same
     * rationale as the main SHUD() path above — single nfcall.txt next
     * to cvode_stats.txt, NOT a 15-key column. MD->nFCall is the single
     * global free-running counter (incremented inside f() at Model/f.cpp:61);
     * f_surf/f_unsat/f_gw/f_river/f_lake increment the alt counters
     * nFCall1..5 which are NOT shipped per spec (out of scope). */
    {
        char nfcall_path[MAXLEN];
        snprintf(nfcall_path, sizeof(nfcall_path), "%s/nfcall.txt",
                 fout->outpath);
        FILE *nfcall_fp = fopen(nfcall_path, "w");
        if (nfcall_fp != NULL) {
            fprintf(nfcall_fp, "nFCall=%lu\n", MD->nFCall);
            fclose(nfcall_fp);
        } else {
            fprintf(stderr,
                    "[shud] WARN: nfcall.txt fopen failed at '%s'; "
                    "stdout-only fallback used.\n", nfcall_path);
        }
    }

    /* P8-tune.G0 PR-B (#411): drain wrapper telemetry from LS1 (the
     * surface solver — same representative-mem1 rationale as the
     * cvode_stats.txt emission above). The other LS{2..5} are released
     * below but their telemetry is not drained (uncouple path is not a
     * G0 target; G0 covers the implicit / coupled path). */
    {
        const char *telemetry_tsv = getenv("SHUD_TELEMETRY_TSV");
        if (telemetry_tsv != NULL && telemetry_tsv[0] != '\0' && LS1 != NULL) {
            FILE *tsv_fp = fopen(telemetry_tsv, "w");
            if (tsv_fp != NULL) {
                int written = SUNLinSol_Hypre_DrainTelemetry(LS1, tsv_fp);
                fclose(tsv_fp);
                fprintf(stdout,
                        "[shud-G0] SUNLinSol_Hypre_DrainTelemetry "
                        "(uncouple LS1): wrote %d telemetry rows to %s\n",
                        written, telemetry_tsv);
                fflush(stdout);
            } else {
                fprintf(stderr,
                        "[shud-G0] WARN: SHUD_TELEMETRY_TSV fopen failed "
                        "at '%s'; telemetry not drained (run continues).\n",
                        telemetry_tsv);
            }
        }
    }

    /* Free integrator memory */
    CVodeFree(&mem1);
    CVodeFree(&mem2);
    CVodeFree(&mem3);
    CVodeFree(&mem4);
    CVodeFree(&mem5);

    /* P8-tune.G0 PR-B (#411) + PR-0 #414 Phase 7 finding #4 — release
     * all 5 SUNLinearSolvers (each holds either an AMG wrapper context
     * or stock SPGMR workspace). */
    if (LS1 != NULL) SUNLinSolFree(LS1);
    if (LS2 != NULL) SUNLinSolFree(LS2);
    if (LS3 != NULL) SUNLinSolFree(LS3);
    if (LS4 != NULL) SUNLinSolFree(LS4);
    if (LS5 != NULL) SUNLinSolFree(LS5);

#ifdef SHUD_ENABLE_PROFILE
    /* S0-8a / openMP #10 — profile bucket dump (uncouple path). */
    {
        char prof_path[MAXLEN];
        snprintf(prof_path, sizeof(prof_path), "%s/profile_B0.yaml",
                 fout->outpath);
        shud_profile::dump(prof_path);
    }
#endif

    SUNContext_Free(&sunctx1);
    SUNContext_Free(&sunctx2);
    SUNContext_Free(&sunctx3);
    SUNContext_Free(&sunctx4);
    SUNContext_Free(&sunctx5);
    
    delete MD;
    return ret;
}

int SHUD(int argc, char *argv[]){
    CommandIn CLI;
    FileIn *fin = new FileIn;
    FileOut *fout = new FileOut;
    CLI.parse(argc, argv);
    CLI.setFileIO(fin, fout);
    if(global_implicit_mode){
        SHUD(fin, fout);
    }else{
        SHUD_uncouple(fin, fout);
    }
    delete fin;
    delete fout;
    return 0;
}

