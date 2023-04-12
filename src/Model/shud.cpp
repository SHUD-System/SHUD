#include <stdio.h>
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
using namespace std;
double SHUD(FileIn *fin, FileOut *fout){
    double ret = 0.;
    Model_Data  *MD;        /* Model Data                */
    N_Vector    udata;
    N_Vector    du;
    
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
#ifdef _OPENMP_ON
    omp_set_num_threads(MD->CS.num_threads);
    screeninfo("\nopenMP: ON. No of Threads = %d\n", MD->CS.num_threads);
    udata = N_VNew_OpenMP(NY, MD->CS.num_threads, sunctx);
    du = N_VNew_OpenMP(NY, MD->CS.num_threads, sunctx);
#else
    screeninfo("\nopenMP: OFF\n");
    udata = N_VNew_Serial(NY, sunctx);
    du = N_VNew_Serial(NY, sunctx);
#endif
    screeninfo("\nGlobal Implicit Mode: ON\n");
    MD->LoadIC();
    MD->SetIC2Y(udata);
    MD->initialize_output();
    MD->PrintInit(fout->Init_bak, 0);
    MD->InitFloodAlert(fout->floodout);
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
            MD->updateforcing(t);
            /* calculate Interception Storage */
            MD->ET(t, tnext);
            if(dummy_mode){
                t = tnext;  /* dummy mode only. */
            }else{
                flag = CVode(mem, tnext, udata, &t, CV_NORMAL);
                check_flag(&flag, "CVode", 1);
            }
        }
        //            CVODEstatus(mem, udata, t);
        MD->summary(udata);
        MD->CS.ExportResults(t);
        MD->flood->FloodWarning(t);
    }
    MD->ScreenPrint(t, MD->CS.NumSteps);
    MD->PrintInit(fout->Init_update, t);
    MD->modelSummary(1);
    /* Free memory */
    N_VDestroy_Serial(udata);
    N_VDestroy_Serial(du);
    /* Free integrator memory */
    CVodeFree(&mem);

    SUNContext_Free(&sunctx);
    delete MD;
    return ret;
}


double SHUD_uncouple(FileIn *fin, FileOut *fout){
    double ret = 0.;
    Model_Data  *MD;        /* Model Data                */
    N_Vector    u1, u2, u3, u4, u5;
    N_Vector    du1, du2, du3, du4, du5;
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
    /* Free memory */
    N_VDestroy_Serial(u1);
    N_VDestroy_Serial(u2);
    N_VDestroy_Serial(u3);
    N_VDestroy_Serial(u4);
    N_VDestroy_Serial(u5);
    
    N_VDestroy_Serial(du1);
    N_VDestroy_Serial(du2);
    N_VDestroy_Serial(du3);
    N_VDestroy_Serial(du4);
    N_VDestroy_Serial(du5);
    /* Free integrator memory */
    CVodeFree(&mem1);
    CVodeFree(&mem2);
    CVodeFree(&mem3);
    CVodeFree(&mem4);
    CVodeFree(&mem5);

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

