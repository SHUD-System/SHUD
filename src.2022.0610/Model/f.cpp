#include "f.hpp"
int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS){
    double       *Y, *DY;
    Model_Data      * MD;
    MD = (Model_Data *) DS;
    timeNow = t;
#ifdef _OPENMP_ON
    Y = NV_DATA_OMP(CV_Y);
    DY = NV_DATA_OMP(CV_Ydot);
    MD->f_update_omp(Y, DY, t);
    MD->f_loop_omp(Y, DY, t);
    MD->f_applyDY_omp(DY, t);
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
    MD->f_update(Y, DY, t);
    MD->f_loop(t);
    MD->f_applyDY(DY, t);
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
