#include "f.hpp"
#ifdef SHUD_ENABLE_PROFILE
#include "timer.h"
#endif
/* S1d.1 (openMP #47): MD_rhs_core.hpp is now unconditionally pulled in
 * for the default (LEGACY_RHS=0) serial dispatch path. LEGACY_RHS=1
 * builds also include it harmlessly — the header is `Model_Data`-only
 * and the ExecPolicy enum is data-only; no runtime cost. */
#include "MD_rhs_core.hpp"
int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
#ifdef SHUD_ENABLE_PROFILE
    /* S0-10 / openMP #14 — t_RHS_total wraps the entire outer RHS
     * callback: vector unwrap + sub-call dispatch + nFCall accounting +
     * (debug-only) DY snapshot. The inner t_RHS_kernel timer below
     * scopes ONLY the three flux-computation sub-calls (update / loop /
     * applyDY), so kernel ≤ total; the difference is "outer overhead"
     * (NV_DATA dereference, counter bump, etc.). Both #ifdef-guarded so
     * PROFILE=0 builds are bitwise-equivalent to the un-instrumented
     * baseline (timer.h header-only no-ops, but the RAII object would
     * still emit an extra ctor/dtor frame at -O0 — guards keep the
     * release build clean too). */
    shud_profile::Timer _t_rhs_total("t_RHS_total");
#endif
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
    timeNow = t;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
    {
#ifdef SHUD_ENABLE_PROFILE
        shud_profile::Timer _t_rhs_kernel("t_RHS_kernel");
#endif
        MD->f_update_omp(Y, DY, t);
        MD->f_loop_omp(Y, DY, t);
        MD->f_applyDY_omp(DY, t);
    }
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
    /* Debug Code
    N_VectorContent_Serial x;
    x =(N_VectorContent_Serial)(CV_Y->content);
    printf("%f\n", x->data[0 + 2 * MD->NumEle]);
    printf("%f\n", x->data[0 + 3 * MD->NumEle]);
    printf("%f\n", x->data[0 + 3 * MD->NumEle + MD->NumRiv]);
     */
    {
#ifdef SHUD_ENABLE_PROFILE
        shud_profile::Timer _t_rhs_kernel("t_RHS_kernel");
#endif
        /* S1d.1 (openMP #47): LEGACY_RHS = 0 (default) routes to the
         * extracted B1a path via `rhs_core(ExecPolicy::Serial)`;
         * LEGACY_RHS = 1 routes to the original `f_update/f_loop/
         * f_applyDY` chain (B0 binary path). The prior S1a scaffold
         * macro has been retired in the SAME atomic commit that
         * introduces LEGACY_RHS, per spec exec-policy-enum Scenario
         * "S1d.1 Step 4.4 + 4.5 atomic landing enforced at review
         * time" — keeping legacy reachable via a recompile and
         * preventing any intermediate state where the scaffold is
         * removed but LEGACY_RHS is not yet wired up.
         *
         * The `#ifdef _OPENMP_ON` branch above is UNCHANGED on
         * purpose: legacy `_omp` dispatch belongs to S2 / #48 and S1
         * is forbidden from touching `f_*_omp` source per master plan
         * §C1 (`_omp` path frozen through S1d). */
#ifdef LEGACY_RHS
        MD->f_update(Y, DY, t);
        MD->f_loop(t);
        MD->f_applyDY(DY, t);
#else
        MD->rhs_core(Y, DY, t, ExecPolicy::Serial);
#endif
    }
#endif
    MD->nFCall++;
#ifdef DEBUG
    printDY(MD->file_debug, DY, MD->NumY, t);
#endif
    return 0;
}

int f_surf(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    timeNow = t;
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
#endif
//printf("f_surf t0=%f, t1=%f, t=%f\n", MD->t0, MD->t1, t);
    MD->f_updatei(Y, DY, t, 1);
    MD->f_loopET(t);
    MD->f_loop1(t);
    MD->f_applyDYi(DY, t, 1);
    MD->nFCall1++;
    return 0;
}
int f_unsat(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    timeNow = t;
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
#endif
    MD->f_updatei(Y, DY, t, 2);
    MD->f_loop2(t);
    MD->f_applyDYi(DY, t, 2);
    MD->nFCall2++;
    return 0;
}
int f_gw(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    timeNow = t;
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
#endif
    MD->f_updatei(Y, DY, t, 3);
    MD->f_loop3(t);
    MD->f_applyDY_gw(DY, t);
    MD->nFCall3++;
    return 0;
}
int f_river(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    timeNow = t;
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
#endif
    MD->f_updatei(Y, DY, t, 4);
    MD->f_loop4(t);
    MD->f_applyDYi(DY, t, 4);
    MD->nFCall4++;
    return 0;
}
int f_lake(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    timeNow = t;
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
#else
    Y = NV_DATA_S(CV_Y);
    DY = NV_DATA_S(CV_Ydot);
#endif
    MD->f_updatei(Y, DY, t, 5);
    MD->f_loop5(t);
    MD->f_applyDYi(DY, t, 5);
    MD->nFCall5++;
    return 0;
}
