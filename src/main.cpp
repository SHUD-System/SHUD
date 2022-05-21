/*******************************************************************************
 * File        : main.cpp                                                      *
 * Version     : June, 2022 (SHUD v2.0)                                         *
 * Function    : SHUD (Simulator for Hydrologic Unstructured Domains)          *
 * Website     : https://www.shud.xyz/
 * Maintainer  : Lele Shu (lele.shu@gmail.com)
 * Developer of SHUD 2.0 (2022):        Lele Shu (lele.shu@gmail.com)          *
 * Developer of SHUD 1.0 (2019):        Lele Shu (lele.shu@gmail.com)          *
 *                                                                             *
 * SHUD inherits the fundamental idea of solving hydrological variables in     *
 * CVODE, but it is imcompatible with PIHM input/output, neither the algorithm.*
 *-----------------------------------------------------------------------------*
 * Developer of PIHM 3.0:        Gopal Bhatt (gopal.bhatt@psu.edu)             *
 * Developer of PIHM 2.2:        Xuan Yu                                       *
 * Developer of PIHM 2.0:        Mukesh Kumar                                  *
 * Developer of PIHM-hydro:      Shuangcai Li                                  *
 * Developer of PIHM 1.0:        Yizhong Qu   (quyizhong@gmail.com)            *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *******************************************************************************/
#include "shud.hpp"
#include "print.hpp"
int main(int argc, char *argv[]){
    SHUDlogo();
    SHUD(argc, argv);
    return 0;
}

