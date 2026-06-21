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
                MD->updateforcing(t);
                /* calculate Interception Storage */
                MD->ET(t, tnext);
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
            //            CVODEstatus(mem, udata, t);
            MD->summary(udata);
            MD->CS.ExportResults(t);
            MD->flood->FloodWarning(t);
        }
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

    /* Free integrator memory */
    CVodeFree(&mem);

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

    /* Free integrator memory */
    CVodeFree(&mem1);
    CVodeFree(&mem2);
    CVodeFree(&mem3);
    CVodeFree(&mem4);
    CVodeFree(&mem5);

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

