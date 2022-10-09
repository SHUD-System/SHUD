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
void PrintFinalStats(void *cvode_mem);

// void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS); // CVODE 5.X
void SetCVODE(void * &cvode_mem, CVRhsFn f, Model_Data *MD,  N_Vector udata, SUNLinearSolver &LS, SUNContext &sunctx); // CVODE 6.X
void CVODEstatus(void *cvode_mem, N_Vector u, realtype t);

#endif				/* cvode_config_h */
