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
void Model_Data:: f_loop(double t){
    /* S5d.1 (#178) — all Ele[i].<hot-field> reads rerouted to
     * hot.<field>[i]. The four _Element writer methods
     * (updateLakeElement / updateElement / Flux_Infiltration via
     * fun_Ele_Infiltraion / Flux_Recharge via fun_Ele_Recharge) keep
     * AoS dispatch; each is followed by sync_hot_dynamic(i) so the
     * dynamic SoA subset (u_qi, u_qex, u_effKH, u_satn) reflects the
     * post-write AoS value before subsequent reads. */
    int i;
    for (i = 0; i < NumEle; i++) {
        if(lakeon && hot.iLake[i] > 0){
            /* Lake elements */
            Ele[i].updateLakeElement();
            sync_hot_dynamic(i);
            fun_Ele_lakeVertical(i, t);
            /* S3b.4 (PR-9): shared writes extracted to per-element slots.
             * See MD_rhs_core.cpp::rhs_flux for full rationale (mirror
             * change here so the dead-code f_loop stays semantically
             * identical to the active rhs_flux). */
            qEleEvapo_lake[i] = qEleEvapo[i] / lake[hot.iLake[i] - 1].NumEleLake;
            qElePrep_lake[i]  = qElePrep[i]  / lake[hot.iLake[i] - 1].NumEleLake;
        }else{
            f_etFlux(i, t);
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            sync_hot_dynamic(i);
            fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
            fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
        }
    }
    for (i = 0; i < NumEle; i++) {
        if(lakeon && hot.iLake[i] > 0){
            /* Lake elements */
            fun_Ele_lakeHorizon(i, t);
        }else{
            /*========surf/gw flow Function==============*/
            fun_Ele_surface(i, t);  // AFTER infiltration, do the lateral flux. ESP for overland flow.
            fun_Ele_sub(i, t);
        }
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    /* S3b.4 (PR-9): transitional gather; mirrors rhs_flux. Must run
     * BEFORE the lake clamp. */
    if(lakeon){
        for (i = 0; i < NumLake; i++) {
            qLakeEvap[i] = 0.;
            qLakePrcp[i] = 0.;
        }
        for (i = 0; i < NumEle; i++) {
            if(hot.iLake[i] > 0){
                int ilake = hot.iLake[i] - 1;
                qLakeEvap[ilake] += qEleEvapo_lake[i];
                qLakePrcp[ilake] += qElePrep_lake[i];
            }
        }
    }
    for (i = 0; i < NumLake; i++) {
        qLakeEvap[i] = min(qLakeEvap[i], qLakePrcp[i] + yLakeStg[i]);
        qLakeEvap[i] = max(0, qLakeEvap[i]);
    }
    /* #43 (S1-pre-B): before-PassValue_legacy probe. Dumps Qe2r_Surf
     * (length NumEle), a PassValue_legacy() write-set member that carries
     * the previous iteration's element-to-river surface flux state.
     * PassValue_legacy (see body at L182-205) zero-resets Qe2r_Surf[0..NumEle-1]
     * and then accumulates QsegSurf over NumSegmt segments; capturing
     * Qe2r_Surf HERE gives a deterministic snapshot of the value that
     * is about to be cleared + re-derived by PassValue_legacy. PR #54 round-1
     * fix F4 replaced the prior QeleSurfTot probe payload (which
     * f_update zero-resets so the snapshot was always all zeros and
     * thus useless as a before-vs-after PassValue_legacy diff) with this
     * write-set member. Site tag "f_loop_before_passvalue" is distinct
     * from the no-op "f_loop" hook below + the "f_update" hook in
     * MD_update.cpp:151; the writer SHUD_DUMP_FNAME_SUFFIX env
     * disambiguates output files (`snapshot_t<v>_before_passvalue.bin`
     * vs `snapshot_t<v>.bin`). SHUD_DUMP_RHS=0 builds emit zero code
     * (compile-switch neutrality contract). */
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop_before_passvalue", t, Qe2r_Surf, NumEle);
#endif
    /* Shared for both OpenMP and Serial, to update.
     * S3c.3 (PR-11 #155): PassValue_legacy() retired, replaced by
     * rhs_deterministic_gather() (MD_rhs_core.cpp); f_loop is
     * legacy dead code kept as a mirror of rhs_flux. */
    rhs_deterministic_gather();
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop", t, NULL, 0);
#endif
}

void Model_Data::f_applyDY(double *DY, double t){
    /* S5d.1 (#178) — Ele[i].{area, iBC, QBC, iSS, QSS, Sy, iLake} reads
     * rerouted to SoA mirror. No writes to AoS here. */
    double area;
    int isf, ius, igw;
    for (int i = 0; i < NumEle; i++) {
        isf = iSF; ius = iUS; igw = iGW;
        area = hot.area[i];
        QeleSurfTot[i] = Qe2r_Surf[i];
        QeleSubTot[i] = Qe2r_Sub[i];
        /* S5d.2-5a (#179) — flat read via accessor; bitwise equivalent
         * to former `QeleSurf[i][j]` because the underlying storage is
         * one contiguous double[NumEle*3] with at(i,j) ↔ `_flat[3*i + j]`. */
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurfAt(i, j);
            QeleSubTot[i] += QeleSubAt(i, j);
            CheckNANij(QeleSurfAt(i, j), i, "QeleSurfAt(i, j)");
            CheckNANij(QeleSubAt(i, j), i, "QeleSubAt(i, j)");
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area - qEs[i];
        DY[ius] = qEleInfil[i] - qEleRecharge[i] - qEu[i] - qTu[i];
        DY[igw] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area - qEg[i] - qTg[i];
//        if(uYgw[i] > hot.WetlandLevel[i] ){ /* IF GW above surface, exfil from GW, OR, from UNSAT*/
//            DY[igw] += - qEleExfil[i];
//        }else{
//            if(uYus[i] > 0.){
//                DY[ius] += - qEleExfil[i];
//            }else{
//                DY[igw] += - qEleExfil[i];
//            }
//        }
        /* Boundary condition and Source/Sink */
        if(hot.iBC[i] == 0){
        }else if(hot.iBC[i] > 0){ // Fix head of GW.
            DY[igw] = 0;
        }else if(hot.iBC[i] < 0){ // Fix flux in GW
            DY[igw] += hot.QBC[i] / area;
        }

        if(hot.iSS[i] == 0){
        }else if(hot.iSS[i] > 0){ // SS in Landusrface
            DY[isf] += hot.QSS[i] / area;
        }else if(hot.iSS[i] < 0){ // SS in GW
            DY[igw] += hot.QSS[i] / area;
        }
        /* Convert with specific yield */
        DY[ius] /= hot.Sy[i];
        DY[igw] /= hot.Sy[i];

//        if(i+1==ID_ELE && DY[i]+uYsf[i] > 0.1 ){  // debug only
//            printf("%.3f, %d: %.2e, %.2e, %.2e | (%.2e, %.2e, %.2e, %.2e, %.2e )\n",
//                   t, i+1,
//                   uYsf[i], uYus[i], uYgw[i], qEleInfil[i], - qEleRecharge[i], - qEu[i],  - qTu[i], QeleSurfTot[i] / area);
//            printf("\n");
//        }
//        if(i +1== ID_ELE){  // debug only
//            printf("%.3f, %d: %f, %f, %f | (%f, %f, %f), %.2e, %.2e\n",
//                   t, i+1,
//                   DY[i], DY[ius], DY[igw], QeleSurf[i][0] / area, QeleSurf[i][1] / area, QeleSurf[i][2] / area,
//                   -QeleSurfTot[i] / area, qEleRecharge[i]);
//            printf("\n");
//        }
        if(hot.iLake[i] > 0){
            DY[i] = 0.;
            DY[ius] = 0.;
            DY[igw] = 0.;
        }
#ifdef DEBUG
        CheckNANi(DY[i], i, "DY[i] (Model_Data::f_applyDY)");
        CheckNANi(DY[ius], i, "DY[ius] (Model_Data::f_applyDY)");
        CheckNANi(DY[igw], i, "DY[igw] (Model_Data::f_applyDY)");
#endif
    }
    for (int i = 0; i < NumRiv; i++) {
        if(Riv[i].BC > 0){
//            Newmann condition.
            DY[iRIV] = 0.;
        }else{
            DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].Length; // dA on CS
            if(DY[iRIV] < -1. * Riv[i].u_CSarea){ /* The negative dA cannot larger then Availalbe Area. */
                DY[iRIV] = -1. * Riv[i].u_CSarea;
            }
            DY[iRIV] = fun_dAtodY(DY[iRIV], Riv[i].u_topWidth, Riv[i].bankslope);
            
//            if(i+1 == ID_RIV && fabs(DY[iRIV])> 0.001){
//                printf("%d, up:%f, down:%f, sub:%f, surf:%f, dy=%f\n", i+1,
//                       - QrivUp[i] / Riv[i].u_TopArea, - QrivDown[i] / Riv[i].u_TopArea,
//                       - QrivSub[i] / Riv[i].u_TopArea, - QrivSurf[i] / Riv[i].u_TopArea,
//                       DY[iRIV]);
//                i=i;
//            }
        }
#ifdef DEBUG
        CheckNANi(DY[i + 3 * NumEle], i, "DY[i] of river (Model_Data::f_applyDY)");
#endif
    }
    for(int i = 0; i < NumLake; i++){
//        DY[i + 3 * NumEle + NumRiv]
        DY[iLAKE] = qLakePrcp[i] - qLakeEvap[i]  +
                    (QLakeRivIn[i] - QLakeRivOut[i] + QLakeSub[i] + QLakeSurf[i] ) / y2LakeArea[i] ;        
//        if(fabs(DY[iLAKE]) > 1.0e-4){
//            printf("%f: %g + %g\n",t, yLakeStg[i], DY[iLAKE]);
//            i=i;
//        }
#ifdef DEBUG
        CheckNANi(DY[iLAKE], i, "DY[i] of LAKE (Model_Data::f_applyDY)");
#endif
    }
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_applyDY", t, DY, 3 * NumEle + NumRiv + NumLake);
#endif
}

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
