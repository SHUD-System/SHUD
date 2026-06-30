//  Model_Control.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef Model_Control_hpp
#define Model_Control_hpp

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Macros.hpp"
#include "ModelConfigure.hpp"

class Print_Ctrl{
private:
    char    filename[MAXLEN] = {};
    char    header[1024] = {};
    int     Interval = NA_VALUE;
    int     NumVar = NA_VALUE;
    int     NumUpdate = 0;
    double  *icol = NULL;
    double  **PrintVar = NULL;
    double  *buffer = NULL;
    int     Binary = 1;
    int     Ascii = 0;
    double  tau = 1440.;    // time unit in calculation. [min]
    FILE    *fid_bin = NULL;
    FILE    *fid_asc = NULL;
    char    filea[MAXLEN] = {};
    char    fileb[MAXLEN] = {};
    long    StartTime = 0;
public:
    Print_Ctrl();
    Print_Ctrl(const Print_Ctrl&) = delete;
    Print_Ctrl& operator=(const Print_Ctrl&) = delete;
    ~Print_Ctrl();
    void    open_file(int a, int b);
    void    PrintData (double dt, double t);
    void    setHeader(const char *s);
    void    Init(long st, int n, const char *s, int dt, double *x, int iFlux);
    void    Init(long st, int n, const char *s, int dt, double *x, int iFlux, int *flag);
    /* S5d.2-5a (#179) — flat-array InitIJ overloads for the
     * QeleSurf/QeleSub flattening. `x_flat` is a contiguous
     * `double[n*3]` block in row-major order at(i,j) ↔
     * x_flat[3*i + j]; `j` (0..2) picks the column the PrintCtrl
     * slot reads. Used by MD_initialize.cpp for ele_Q_sub{0,1,2}
     * and ele_Q_surf{0,1,2} files. The legacy `double**` InitIJ
     * overloads were removed in PR #197 (review A-S1): grep showed
     * zero remaining callers after the jagged→flat refactor. */
    void    InitIJ(long st, int n, const char *s, int dt, double *x_flat, int j, int iFlux);
    void    InitIJ(long st, int n, const char *s, int dt, double *x_flat, int j, int iFlux, int *flag);
private:
    void    fun_printASCII(double t, double dt);
    void    fun_printBINARY(double t, double dt);
    void    close_file();
};
class PrintOutDt {
public:
    /* Element storage */
    int dt_ye_gw = 0;
    int dt_ye_surf = 0;
    int dt_ye_snow = 0;
    int dt_ye_ic = 0;
    int dt_ye_unsat = 0;
    
    /* Element Fluxes */
    int dt_qe_prcp = 1440; // default output PRCP.
    int dt_qe_infil = 0;
    int dt_qe_et = 0;
    int dt_qe_rech = 0;
    int dt_qe_etp = 0;
    int dt_qe_eta = 0;
    
    /* Element volume Fluxes */
    int dt_Qe_sub = 0;
    int dt_Qe_subx = 0;
    int dt_Qe_surf = 0;
    int dt_Qe_surfx = 0;
    int dt_Qe_rsub = 0;
    int dt_Qe_rsurf = 0;
    
    /* River Stage */
    int dt_yr_stage = 0;
    /* River volume Fluxes */
    int dt_Qr_up = 0;
    int dt_Qr_down = 0;
    int dt_Qr_sub = 0;
    int dt_Qr_surf = 0;
    
    /* Lake  stage */
    int dt_lake     = 1440;
    
    void calibmode(int dt);
    void defaultmode();
};
class Control_Data : public PrintOutDt{
private:
    double DayStart = 0;    /* START Day [Day] */
    double DayEnd = 10;     /* END Day [Day] */
public:
    int Verbose = 0;
    int CloseBoundary = 1; /* Whether the close boundary exist. 1 = no boundary flux. 0 = Default boundary flux. [bool]*/
    int Ascii = 0;      /* Whether export result as ASCII File [bool]*/
    int Binary = 1;     /* Whether export result as Binary File [bool]*/
    int Spinup = 0;     /* Number of days for spinup */
    int screenIntv = 1440;
    //    int Solver;    /* Solver type */
    unsigned long NumSteps;    /* Number of external time steps
                      * (when results can be printed) for
                      * the whole simulation [-]*/
    int num_threads;    /* Number of Threads in OPENmp only [-]*/
    int init_type = 3;    /* initialization mode [-]*/
    int cryosphere = 0;
    double abstol = 1.0e-4;    /* absolute tolerance [-]*/
    double reltol = 1.0e-3;    /* relative tolerance [-]*/
    double InitStep = 1.e-2;    /* initial step size [min]*/
    double MaxStep = 30;       /* Maximum step size [min] */
    double SolverStep = 2;       /* Maximum step size [min] */
    int UpdateICStep = 1440;       /* Maximum step size [min] */
    double ETStep = 60;         /* Step for et from interception [min]*/
    
    int     ET_Mode = 0; /* 0 - PET_Penman_Monteith
                            1 - PET_Hargreaves (Rn, Tmax, Tmin)
                            2 - PET_Priestley_Taylor
                          */
    
    double StartTime = 0.;      /* Start time of simulation [min]*/
    double EndTime = 14400;     /* End time of simulation [min]*/
    double dt = 1;
    /* #401 sub-task 1 — NSDMI nullptr default so the dtor's `delete[]`
     * is a defined no-op when this field is never allocated. As of
     * 2026-06-30 `Tout` is a dead field (no assignment anywhere in
     * src/), but the NSDMI + delete[] pair is kept defensive for any
     * future writer. */
    double *Tout = nullptr;
    int NumPrint = 0;;
    int exfiltration = 0;
    Print_Ctrl PCtrl[100];
    
    /* Methods */
    Control_Data();
    ~Control_Data();
    void ExportResults(double t);
    void updateSimPeriod(double day0, double day1);
    void read(const char *fn);
    void write(const char *fn);
    void getValue(const char *varname);
private:
    void updateSimPeriod();
} ;

#endif /* Model_Control_hpp */
