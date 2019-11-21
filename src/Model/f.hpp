#ifndef f_hpp
#define f_hpp
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ModelConfigure.hpp"
#include "Macros.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"

int f(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);

int f_surf(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);
int f_unsat(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);
int f_gw(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);
int f_river(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);
int f_lake(double t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);

//int f_explicit(double t, N_Vector CV_Y, void *DS, double dt);
//void f_applyDY(Model_Data * MD, double *DY);
//void f_loop(Model_Data * MD, double  *Y, double  *DY, double t);
//void f_DY2Y(double *Y, double *DY, double dt);
//void fSolver(Model_Data *md, N_Vector udata, double t, double tout);
#endif /* f_h */
