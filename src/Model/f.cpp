#include "f.hpp"
#ifdef SHUD_ENABLE_PROFILE
#include "timer.h"
#endif
/* MD_rhs_core.hpp declares the coupled RHS (rhs_core) that f()
 * dispatches to. */
#include "MD_rhs_core.hpp"
int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
#ifdef SHUD_ENABLE_PROFILE
    /* t_RHS_total wraps the entire RHS callback: vector unwrap +
     * rhs_core dispatch + nFCall accounting + (debug-only) DY
     * snapshot. The inner t_RHS_kernel timer below scopes ONLY the
     * rhs_core call (update / flux / apply), so kernel ≤ total; the
     * difference is callback overhead (array-pointer lookup, counter
     * bump, etc.). Both are #ifdef-guarded so that builds without
     * SHUD_ENABLE_PROFILE contain no timer code at all (not even an
     * empty RAII ctor/dtor frame at -O0). */
    shud_profile::Timer _t_rhs_total("t_RHS_total");
#endif
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
    timeNow = t;
    /* Generic N_Vector data accessor. SUNDIALS `N_VGetArrayPointer`
     * dispatches on the N_Vector's ops table at runtime, so the same
     * source compiles + works against nvector_serial OR
     * nvector_openmp backends (selected by SHUD_USE_OPENMP_NVECTOR
     * via the N_VNew_* dispatch in shud.cpp). All six RHS callbacks
     * in this file (f / f_surf / f_unsat / f_gw / f_river / f_lake)
     * use this form; do not use the type-specific NV_Ith_* / NV_DATA_*
     * macros, which are only valid for one backend. */
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
        /* f() routes to rhs_core (rhs_update -> rhs_flux -> rhs_apply).
         * The execution policy is selected by a build flag:
         *   - SHUD_ENABLE_OPENMP_RHS=0 (default for `make shud`):
         *     ExecPolicy::Serial.
         *   - SHUD_ENABLE_OPENMP_RHS=1 (default for `make shud_omp`):
         *     ExecPolicy::StrictOMP, the single-region OpenMP
         *     implementation in MD_rhs_core.cpp (3 phases with implicit
         *     barriers; the Makefile adds -fopenmp for this setting). */
#ifdef SHUD_ENABLE_OPENMP_RHS
        MD->rhs_core(Y, DY, t, ExecPolicy::StrictOMP);
#else
        MD->rhs_core(Y, DY, t, ExecPolicy::Serial);
#endif
    }
    /* nFCall is SHUD's entry counter for the coupled RHS (the nFCall1..5
     * counters belong to the uncoupled callbacks below). Free-running; it is
     * written to nfcall.txt, separately from the CVODE statistics, and may
     * legitimately differ from CVODE's own nfe count. */
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
