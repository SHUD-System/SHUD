//  MD_f.cpp
//
//  Created by Lele Shu on 1/27/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
void Model_Data:: f_loop(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        f_etFlux(i, t);
        /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
        /*========infiltration/Recharge Function==============*/
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
        fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
        fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
    }
    for (i = 0; i < NumEle; i++) {
        /*========surf/gw flow Function==============*/
        fun_Ele_surface(i, t);  // AFTER infiltration, do the lateral flux. ESP for overland flow.
        fun_Ele_sub(i, t);
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
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
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area - qEs[i];
        DY[ius] = qEleInfil[i] - qEleRecharge[i] - qEu[i] - qTu[i];
        DY[igw] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area - qEg[i] - qTg[i];

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
//        if(DY[ius] < 0. && uYus[i] < -10.){  // debug only
//            printf("%.3f, %d: %.2f, %.2e, %.2f | (%.2e, %.2e, %.2e), %.2e, %.2e\n", t, i+1,
//                   uYus[i], DY[ius], uYgw[i], qEleE_IC[i], qEleTrans[i], qEleEvap[i],
//                   qEleInfil[i], qEleRecharge[i]);
//            printf("\n");
//        }
//                DY[isf] =0.0;  // debug only.
//                DY[ius] =0.0;
//                DY[igw] =0.0;
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
            DY[iRIV] = fun_dAtodY(DY[iRIV], Riv[i].u_topWidth, Riv[i].bankslope);
//            DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].u_TopArea;
        }
//        DY[iRIV] = 0.0;
#ifdef DEBUG
        CheckNANi(DY[i + 3 * NumEle], i, "DY[i] of river (Model_Data::f_applyDY)");
#endif
    }
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
        if(iDownStrm >= 0){
            QrivUp[iDownStrm] += - QrivDown[i];
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
