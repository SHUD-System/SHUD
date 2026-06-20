//  MD_f.cpp
//
//  Created by Lele Shu on 1/27/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif
void Model_Data:: f_loop(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            Ele[i].updateLakeElement();
            fun_Ele_lakeVertical(i, t);
            qLakeEvap[Ele[i].iLake - 1] += qEleEvapo[i] / lake[Ele[i].iLake - 1].NumEleLake;
            qLakePrcp[Ele[i].iLake - 1] += qElePrep[i] / lake[Ele[i].iLake - 1].NumEleLake;
        }else{
            f_etFlux(i, t);
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
            fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
        }
    }
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
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
    for (i = 0; i < NumLake; i++) {
        qLakeEvap[i] = min(qLakeEvap[i], qLakePrcp[i] + yLakeStg[i]);
        qLakeEvap[i] = max(0, qLakeEvap[i]);
    }
    /* #43 (S1-pre-B): before-PassValue probe. Dumps Qe2r_Surf
     * (length NumEle), a PassValue() write-set member that carries
     * the previous iteration's element-to-river surface flux state.
     * PassValue (see body at L182-205) zero-resets Qe2r_Surf[0..NumEle-1]
     * and then accumulates QsegSurf over NumSegmt segments; capturing
     * Qe2r_Surf HERE gives a deterministic snapshot of the value that
     * is about to be cleared + re-derived by PassValue. PR #54 round-1
     * fix F4 replaced the prior QeleSurfTot probe payload (which
     * f_update zero-resets so the snapshot was always all zeros and
     * thus useless as a before-vs-after PassValue diff) with this
     * write-set member. Site tag "f_loop_before_passvalue" is distinct
     * from the no-op "f_loop" hook below + the "f_update" hook in
     * MD_update.cpp:151; the writer SHUD_DUMP_FNAME_SUFFIX env
     * disambiguates output files (`snapshot_t<v>_before_passvalue.bin`
     * vs `snapshot_t<v>.bin`). SHUD_DUMP_RHS=0 builds emit zero code
     * (compile-switch neutrality contract). */
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop_before_passvalue", t, Qe2r_Surf, NumEle);
#endif
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop", t, NULL, 0);
#endif
}

void Model_Data::f_applyDY(double *DY, double t){
    double area;
    int isf, ius, igw;
    for (int i = 0; i < NumEle; i++) {
        isf = iSF; ius = iUS; igw = iGW;
        area = Ele[i].area;
        QeleSurfTot[i] = Qe2r_Surf[i];
        QeleSubTot[i] = Qe2r_Sub[i];
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurf[i][j];
            QeleSubTot[i] += QeleSub[i][j];
            CheckNANij(QeleSurf[i][j], i, "QeleSurf[i][j]");
            CheckNANij(QeleSub[i][j], i, "QeleSub[i][j]");
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area - qEs[i];
        DY[ius] = qEleInfil[i] - qEleRecharge[i] - qEu[i] - qTu[i];
        DY[igw] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area - qEg[i] - qTg[i];
//        if(uYgw[i] > Ele[i].WetlandLevel ){ /* IF GW above surface, exfil from GW, OR, from UNSAT*/
//            DY[igw] += - qEleExfil[i];
//        }else{
//            if(uYus[i] > 0.){
//                DY[ius] += - qEleExfil[i];
//            }else{
//                DY[igw] += - qEleExfil[i];
//            }
//        }
        /* Boundary condition and Source/Sink */
        if(Ele[i].iBC == 0){
        }else if(Ele[i].iBC > 0){ // Fix head of GW.
            DY[igw] = 0;
        }else if(Ele[i].iBC < 0){ // Fix flux in GW
            DY[igw] += Ele[i].QBC / area;
        }

        if(Ele[i].iSS == 0){
        }else if(Ele[i].iSS > 0){ // SS in Landusrface
            DY[isf] += Ele[i].QSS / area;
        }else if(Ele[i].iSS < 0){ // SS in GW
            DY[igw] += Ele[i].QSS / area;
        }
        /* Convert with specific yield */
        DY[ius] /= Ele[i].Sy;
        DY[igw] /= Ele[i].Sy;
        
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
        if(Ele[i].iLake > 0){
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

void Model_Data::PassValue(){
    int i, ie, ir;
    for (i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    for (i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }
    for (i = 0; i < NumSegmt; i++) {
        ie = RivSeg[i].iEle-1;
        ir = RivSeg[i].iRiv-1;
        QrivSurf[ir] += QsegSurf[i]; // Positive from River to Element
        QrivSub[ir] += QsegSub[i];
        Qe2r_Surf[ie] += -QsegSurf[i]; // Positive from Element to River
        Qe2r_Sub[ie] += -QsegSub[i];
    }
    for (i = 0; i < NumRiv; i++) {
        if(iDownStrm >= 0 && Riv[i].toLake <= 0){
            QrivUp[iDownStrm] += - QrivDown[i];
        }
    }
    /* S3b.1 (PR-9): transitional lake gather river -> per-lake.
     * Replaces the shared `QLakeRivIn[Riv[i].toLake] += QrivDown[i]`
     * write that was previously in Flux_RiverDown (MD_RiverFlux.cpp).
     * Reads per-river QrivDown (just computed in Flux_RiverDown loop)
     * and aggregates to per-lake QLakeRivIn. Will be replaced by
     * rhs_deterministic_gather() in S3c (PR-11). */
    if(lakeon){
        for (i = 0; i < NumLake; i++) {
            QLakeRivIn[i] = 0.;
        }
        for (i = 0; i < NumRiv; i++) {
            if(Riv[i].toLake >= 0){
                QLakeRivIn[Riv[i].toLake] += QrivDown[i]; // toLake is 0-indexed
            }
        }
    }
    //    for (i = 0; i < NumEle; i++) { /*Check flux A->B  = Flux B->A*/
    //        for (j = 0; j < 3; j++) {
    //            inabr = Ele[i].nabr[j] - 1;
    //            if (inabr >= 0) {
    //                jnabr = Ele[i].nabrToMe[j];
    //                if(Ele[inabr].iupdSF[jnabr]){
    //                    //void
    //                }else{
    //                    QeleSurf[inabr][jnabr] = - QeleSurf[i][j];
    //                    Ele[inabr].iupdSF[jnabr] = 1;
    //                    QeleSub[inabr][jnabr] = - QeleSub[i][j];
    //                    Ele[inabr].iupdGW[jnabr] = 1;
    //                }
    //            }
    //        }
    //    }
}

void Model_Data::applyBCSS(double *DY, int i){
    if(Ele[i].iBC > 0){ // Fix head of GW.
        DY[iGW] = 0;
    }else if(Ele[i].iBC < 0){ // Fix flux in GW
        DY[iGW] += Ele[i].QBC / Ele[i].area;
    }else{}
    
    if(Ele[i].iSS > 0){ // SS in Landusrface
        DY[iSF] += Ele[i].QSS / Ele[i].area;
    }else if(Ele[i].iSS < 0){ // SS in GW
        DY[iGW] += Ele[i].QSS / Ele[i].area;
    }else{}
}
