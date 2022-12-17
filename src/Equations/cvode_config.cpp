
#include "cvode_config.hpp"

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
void PrintFinalStats(void *cvode_mem)
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
    LS = SUNLinSol_SPGMR(udata, 0, 0, sunctx);
    check_flag((void *)LS, "SUNLinSol_SPGMR", 0);
    
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
