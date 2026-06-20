#include "f.hpp"
#ifdef SHUD_ENABLE_PROFILE
#include "timer.h"
#endif
/* S2 capstone (PR-8): MD_rhs_core.hpp is the only serial dispatch path;
 * the prior legacy-vs-rhs_core fork has been retired. */
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
    /* S1d.2 (openMP #48) — generic N_Vector data accessor. SUNDIALS 6
     * `N_VGetArrayPointer` dispatches on the N_Vector's ops table at
     * runtime, so the same source compiles + works against
     * nvector_serial OR nvector_openmp backends (selected by
     * SHUD_USE_OPENMP_NVECTOR via the N_VNew_* dispatch in shud.cpp).
     * The 6 prior backend-specific data-accessor blocks (f / f_surf /
     * f_unsat / f_gw / f_river / f_lake) have all been collapsed to
     * this generic form. See design.md D5 (N_VGetArrayPointer
     * rationale — generic ops-table dispatch, not the type-specific
     * NV_Ith). */
    Y = N_VGetArrayPointer(CV_Y);
    DY = N_VGetArrayPointer(CV_Ydot);
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
        /* S2 capstone (PR-8): f() always routes to rhs_core (Serial); the
         * legacy `_omp` RHS receivers (MD_f_omp.cpp) and the legacy/rhs_core
         * fork have been retired. PURE CARRY-OVER `rhs_update/rhs_flux/
         * rhs_apply` are byte-for-byte copies of `f_update/f_loop/
         * f_applyDY`. */
        MD->rhs_core(Y, DY, t, ExecPolicy::Serial);
    }
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
    Y = N_VGetArrayPointer(CV_Y);
    DY = N_VGetArrayPointer(CV_Ydot);
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
    Y = N_VGetArrayPointer(CV_Y);
    DY = N_VGetArrayPointer(CV_Ydot);
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
    Y = N_VGetArrayPointer(CV_Y);
    DY = N_VGetArrayPointer(CV_Ydot);
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
    Y = N_VGetArrayPointer(CV_Y);
    DY = N_VGetArrayPointer(CV_Ydot);
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
    Y = N_VGetArrayPointer(CV_Y);
    DY = N_VGetArrayPointer(CV_Ydot);
    MD->f_updatei(Y, DY, t, 5);
    MD->f_loop5(t);
    MD->f_applyDYi(DY, t, 5);
    MD->nFCall5++;
    return 0;
}
