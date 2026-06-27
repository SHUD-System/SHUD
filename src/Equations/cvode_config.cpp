
#include "cvode_config.hpp"
/* S5c-B (#174): RHS 7-bucket + forcing I/O wall-clock timer dumps.
 * Header is empty under default (SHUD_ENABLE_DIAGNOSTICS undefined). */
#include "../Model/MD_diagnostics.hpp"
/* P8-precond-0 (#345): identity preconditioner stub used by
 * CVodeSetPreconditioner below. Wires CVLS PREC_LEFT call path
 * (nps/npe stat accumulation) while keeping P^{-1} = I so B1b
 * bitwise neutrality is preserved. */
#include "MD_precond_identity.h"

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


void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS, SUNContext &sunctx){
    
    int flag;
    
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
    
    //    LS = SUNSPGMR(udata, 0, 0); //v3.x
    /* P8-precond-0 (#345): pretype PREC_NONE → PREC_LEFT so SPGMR
     * invokes the preconditioner installed below (identity stub).
     * Bitwise neutral vs B1b because P^{-1} = I, but the CVLS
     * preconditioner call path (nps/npe stats, t_precond_setup
     * timer) is now exercised. */
    LS = SUNLinSol_SPGMR(udata, PREC_LEFT, 0, sunctx);
    check_flag((void *)LS, "SUNLinSol_SPGMR", 0);

    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    check_flag(&flag, "CVSpilsSetLinearSolver", 1);

    /* P8-precond-0 (#345): register identity preconditioner pair +
     * setup-frequency cap. CVLS preconditioner registration requires
     * the linear-solver memory to be attached first (above), per
     * SUNDIALS 6.0.0 cvode_ls.h docs. LSetupFrequency=50 matches
     * SUNDIALS default and provides explicit-knob evidence for the
     * spike. */
    flag = CVodeSetPreconditioner(cvode_mem, PSetupIdentity, PSolveIdentity);
    check_flag(&flag, "CVodeSetPreconditioner", 1);
    flag = CVodeSetLSetupFrequency(cvode_mem, 50);
    check_flag(&flag, "CVodeSetLSetupFrequency", 1);

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
