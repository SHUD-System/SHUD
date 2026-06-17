//  cvode_config.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef cvode_config_h
#define cvode_config_h

#include <stdio.h>
#include "ModelConfigure.hpp"
#include "Model_Data.hpp"
#include "functions.hpp"
 /******* SUNDIAL 3.0 and above ****************/
#include "cvode/cvode.h"	/* prototypes for CVODE fcts., consts.  */
#include "nvector/nvector_serial.h"	/* access to serial N_Vector            */
#include "sunlinsol/sunlinsol_spgmr.h"	/* access to SPGMR SUNLinearSolver      */
#include "cvode/cvode_spils.h"	/* access to CVSpils interface          */
#include "sundials/sundials_dense.h"	/* use generic dense solver in precond. */
#include "sundials/sundials_types.h"	/* defs. of realtype, sunindextype      */
#include "sundials/sundials_math.h"	/* contains the macros ABS, SUNSQR, EXP */


/*==========cvode flags===============*/
int check_flag(void *flagvalue, const char *funcname, int opt);

/* PrintFinalStats — query SUNDIALS for the final CVODE stat counters
 * (nfe, nfeLS, nni, nli, nsetups, netf, ...) and print them.
 *
 * The original signature (`void *cvode_mem`) prints to stdout only and
 * is kept as an inline back-compat wrapper.
 *
 * The S0-8a signature (`void *cvode_mem, FILE *fout`) ALSO writes a
 * machine-parsable key=value file to `fout` when `fout != NULL`. stdout
 * output is unchanged. When `fout == NULL` the two signatures behave
 * identically — this lets B0 callers persist stats next to the SHUD
 * output dir (see SHUD/src/Model/shud.cpp ~line 114) without disturbing
 * builds that don't ship the new caller. */
void PrintFinalStats(void *cvode_mem, FILE *fout);
inline void PrintFinalStats(void *cvode_mem) {
    PrintFinalStats(cvode_mem, NULL);
}

// void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS); // CVODE 5.X
void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS, SUNContext &sunctx); // CVODE 6.X
void CVODEstatus(void *cvode_mem, N_Vector u, realtype t);

#endif				/* cvode_config_h */
