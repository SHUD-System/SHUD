//  MD_f.cpp
//
//  Created by Lele Shu on 1/27/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
/* S3c.3 (PR-11 #155): PassValue_legacy() retired; f_loop's gather call now
 * dispatches to Model_Data::rhs_deterministic_gather() defined in
 * MD_rhs_core.cpp (member declaration in Model_Data.hpp). f_loop is
 * legacy dead code preserved as a mirror of rhs_flux. */
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif
/* P1d.2.0 PR-C0 (#291): Model_Data f_loop legacy carry-over deleted.
 * Live counterpart is Model_Data::rhs_flux in MD_rhs_core.cpp, which
 * CVODE invokes via f.cpp:54 -> MD->rhs_core(...). The "f_loop" and
 * "f_loop_before_passvalue" SHUD_DUMP_RHS tag strings at
 * MD_rhs_core.cpp:431 are preserved as the golden-file dump contract. */

/* P1d.2.0 PR-C0 (#291): Model_Data f_applyDY legacy carry-over deleted.
 * Live counterpart is Model_Data::rhs_apply in MD_rhs_core.cpp, which
 * CVODE invokes via f.cpp:54 -> MD->rhs_core(...). The "f_applyDY"
 * SHUD_DUMP_RHS tag string at MD_rhs_core.cpp:676 is preserved as
 * the golden-file dump contract. */

/* S3c.3 (PR-11 #155): Model_Data::PassValue_legacy() retired. Its body has
 * been refactored into Model_Data::rhs_deterministic_gather() in
 * MD_rhs_core.cpp (per design.md D12 -- gather lives WITH the RHS
 * core, not in a separate MD_gather.cpp file). The new function
 * consumes the 7 S4 adjacency lists (PR-10) to perform all gather
 * work (segment->river/element, downstream river, lake river-in,
 * lake-bank surf/sub). Bitwise neutrality vs B0 is preserved because
 * the S4 lists are built in B0 ascending array-index order, so the
 * per-accumulator += sequence is identical to the legacy serial
 * iteration. The B.4 qLakeEvap/qLakePrcp per-element->per-lake
 * gather REMAINS in rhs_flux (MD_rhs_core.cpp ~L214-L226) because
 * the lake clamp pass reads those gathered values BEFORE the
 * call to rhs_deterministic_gather() -- ordering constraint
 * documented in PR-9 / S3b.4.
 *
 * Comment retained verbatim about the original f_loop NumEle
 * neighbor sanity loop (commented-out by upstream long before PR-11):
 *     for (i = 0; i < NumEle; i++) { ... } */

void Model_Data::applyBCSS(double *DY, int i){
    /* S5d.1 (#178) — Ele[i].{iBC, QBC, area, iSS, QSS} reads rerouted
     * to SoA mirror. Read-only; no AoS write. */
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
