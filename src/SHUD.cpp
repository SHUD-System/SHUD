/*******************************************************************************
 * File        : main.cpp                                                      *
 * Version     : Nov, 2018 (SHUD v1.0)                                         *
 * Function    : SHUD (Solver for Hydrologic Unstructured Domains)            *
 * Developer of SHUD 1.0:        Lele Shu (lele.shu@gmail.com)                 *
 *                                                                             *
 * SHUD inherits the fundamental idea of solving hydrological variables in     *
 * CVODE, but it is imcompatible with PIHM input/output as well as algorithm.  *
 *-----------------------------------------------------------------------------*
 * Developer of PIHM 3.0:        Gopal Bhatt (gopal.bhatt@psu.edu)             *
 * Developer of PIHM 2.2:        Xuan Yu                                       *
 * Developer of PIHM 2.0:        Mukesh Kumar                                  *
 * Developer of PIHM-hydro:      Shuangcai Li                                  *
 * Developer of PIHM 1.0:        Yizhong Qu   (quyizhong@gmail.com)            *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in SHUD V1.0 ..........................*
 * 0) Change the language and structure of code from C to C++;
 * 1) Update the CVODE from v2.2 to v5.0
 * 2) Support OpenMP Parrallel computing
 * 3) Change the input/output format. Check the Manual of SHUD on github;
 * 4) Change the structure of River
 * 5) The functions to handle the time-series data, including forcing, LAI,
 *    Roughness Length, Boundary Condition, Melt factor
 * 6) Add Lakes into the hydrological process
 *******************************************************************************/
#include "shud.hpp"
#include "print.hpp"
int main(int argc, char *argv[]){
    SHUDlogo();
    SHUD(argc, argv);
    return 0;
}

