//  MD_f.cpp
//
//  Created by Lele Shu on 1/27/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif
/* The coupled-mode RHS is not in this file. It lives in
 * MD_rhs_core.cpp as Model_Data::rhs_update / rhs_flux / rhs_apply,
 * which CVODE reaches through MD->rhs_core(...). Their SHUD_DUMP_RHS
 * site names are "f_update", "f_loop", "f_loop_before_passvalue" and
 * "f_applyDY"; snapshot files are selected by these names, so they
 * must not be renamed.
 *
 * All gather work (segment->river/element, downstream river, lake
 * river-in, lake-bank surf/sub) is done by
 * Model_Data::rhs_deterministic_gather(), also in MD_rhs_core.cpp.
 * It walks the 7 adjacency lists from MD_adjacency.hpp, which are
 * built in ascending array-index order, so each accumulator sees the
 * same += sequence as a plain serial loop and results are bitwise
 * identical to the serial build. The qLakeEvap/qLakePrcp
 * per-element->per-lake gather is NOT part of it: it stays in
 * rhs_flux because the lake clamp pass reads those sums BEFORE
 * rhs_deterministic_gather() is called. */

void Model_Data::applyBCSS(double *DY, int i){
    /* Ele[i].{iBC, QBC, area, iSS, QSS} are read from the SoA
     * mirror. Read-only; no AoS write. */
    if(hot.iBC[i] > 0){ // Fix head of GW.
        DY[iGW] = 0;
    }else if(hot.iBC[i] < 0){ // Fix flux in GW
        DY[iGW] += hot.QBC[i] / hot.area[i];
    }else{}

    if(hot.iSS[i] > 0){ // SS in Landusrface
        DY[iSF] += hot.QSS[i] / hot.area[i];
    }else if(hot.iSS[i] < 0){ // SS in GW
        DY[iGW] += hot.QSS[i] / hot.area[i];
    }else{}
}
