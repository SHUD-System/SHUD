//  MD_uncouple.cpp
//
//  Created by Lele Shu on 9/6/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"

void Model_Data::f_loopET(double t){
    for (int i = 0; i < NumEle; i++) {
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] );
        f_etFlux(i, t);
    }
}
void Model_Data::f_loop1(double t){
    int i;
    int ie, ir;
    for (i = 0; i < NumEle; i++) {
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
        fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
    }
    /* The loop must be run twice becaue require to update infiltration.*/
    for (i = 0; i < NumEle; i++) {
        fun_Ele_surface(i, t);  // AFTER infiltration, do the lateral flux. ESP for overland flow.
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
    }
    for (i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
    }
    for (i = 0; i < NumSegmt; i++) {
        ie = RivSeg[i].iEle-1;
        ir = RivSeg[i].iRiv-1;
        QrivSurf[ir] += QsegSurf[i];
        Qe2r_Surf[ie] += -QsegSurf[i]; // Positive from Element to River
    }
}

void Model_Data:: f_loop2(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
        /*========ET Function==============*/
        //        f_etFlux(i, t); // et is moved out of f_loop()
        /*========infiltration/Recharge Function==============*/
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
        fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
        fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
    }
}

void Model_Data:: f_loop3(double t){
    int i;
    int ie, ir;
//    for (i = 0; i < NumEle; i++) {
//        /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
//        /*========ET Function==============*/
//        //        f_etFlux(i, t); // et is moved out of f_loop()
//        /*========infiltration/Recharge Function==============*/
//        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
//        fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
//    }
    for (i = 0; i < NumEle; i++) {
        fun_Ele_sub(i, t);
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumEle; i++) {
        Qe2r_Sub[i] = 0.;
    }
    for (i = 0; i < NumRiv; i++) {
        QrivSub[i] = 0.;
    }
    for (i = 0; i < NumSegmt; i++) {
        ie = RivSeg[i].iEle-1;
        ir = RivSeg[i].iRiv-1;
        QrivSub[ir] += QsegSub[i];
        Qe2r_Sub[ie] += -QsegSub[i];
    }
}

void Model_Data:: f_loop4(double t){
    int i;
//    for (i = 0; i < NumEle; i++) {
//        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
//    }
//    for (i = 0; i < NumSegmt; i++) {
//        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
//        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
//    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    /* Shared for both OpenMP and Serial, to update */
    for (i = 0; i < NumRiv; i++) {
        if(iDownStrm >= 0){
            QrivUp[iDownStrm] += - QrivDown[i];
        }
    }
//    PassValue();
}
void Model_Data:: f_loop5(double t){
}

void Model_Data::f_applyDY_gw(double *DY, double t){
    double area = 0.;
    for (int i = 0; i < NumEle; i++) {
        area = Ele[i].area;
        QeleSubTot[i] = Qe2r_Sub[i];
        for (int j = 0; j < 3; j++) {
//            checkExchangeValue(QeleSub, i, j, Ele[i].nabr[j]-1, Ele[i].nabrToMe[j]-1);
            QeleSubTot[i] += QeleSub[i][j];
        }
        DY[i] = (qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area );
        if(uYsf[i] < EPSILON){ /* NO ponding water*/
            if (uYgw[i] < Ele[i].WetlandLevel){
            }else{
                /*Evaporate from Ground water*/
                DY[i] += -qEleEvapo[i];
            }
            
        }else{
            
        }
        if (uYgw[i] > Ele[i].RootReachLevel) {
            /*Vegetation sucks water from Ground water*/
            DY[i] += - qEleTrans[i];
        } else {
        }
        /* Boundary condition and Source/Sink */
        if(Ele[i].iBC > 0){ // Fix head of GW.
            DY[i] = 0;
        }else if(Ele[i].iBC < 0){ // Fix flux in GW
            DY[i] += Ele[i].QBC / area;
        }else{ /* Void */}

        if(Ele[i].iSS > 0){
        }else if(Ele[i].iSS < 0){ // SS in GW
            DY[i] += Ele[i].QSS / area;
        }else{}
        /* Convert with specific yield */
        DY[i] /= Ele[i].Sy;
#ifdef DEBUG
        CheckNANi(uYgw[i], i, "uYgw[i] of Groundwater (Model_Data::f_applyDYi)");
        CheckNANi(DY[i], i, "DY[i] of Groundwater (Model_Data::f_applyDYi)");
#endif
        
    }
}
void Model_Data::f_applyDYi(double *DY, double t, int flag){
//    printf("f_applyDYi t=%f\n", t);
    double area = 0.;
    if(flag ==1){
        for (int i = 0; i < NumEle; i++) {
            area = Ele[i].area;
            QeleSurfTot[i] = Qe2r_Surf[i];
            for (int j = 0; j < 3; j++) {
                QeleSurfTot[i] += QeleSurf[i][j];
            }
            DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area;
            if(Ele[i].iSS > 0){ // SS in Landusrface
                DY[i] += Ele[i].QSS / area;
            }
#ifdef DEBUG
            CheckNANi(uYsf[i], i, "uYsf[i] of Surface (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of Surface (Model_Data::f_applyDYi)");
#endif
        }
    }else if(flag ==2){
        for (int i = 0; i < NumEle; i++) {
            DY[i] = qEleInfil[i] - qEleRecharge[i];
            DY[i] += -qEleEvapo[i];
            if (uYgw[i] > Ele[i].RootReachLevel) {
            } else {
                DY[i] += - qEleTrans[i];
            }
            /* Convert with specific yield */
            DY[i] /= Ele[i].Sy;
#ifdef DEBUG
            CheckNANi(uYus[i], i, "uYus[i] of Unsat (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of Unsat (Model_Data::f_applyDYi)");
#endif
        }
    }else if(flag ==3){
        
    }else if(flag ==4){
        for (int i = 0; i < NumRiv; i++) {
            if(Riv[i].BC > 0){
                //            Newmann condition.
                DY[i] = 0.;
            }else{
                DY[i] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].u_TopArea;
            }
#ifdef DEBUG
CheckNANi(uYriv[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
CheckNANi(DY[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
#endif
        }
    }else if(flag ==5){
        for (int i = 0; i < NumLake; i++) {
            
#ifdef DEBUG
            CheckNANi(uYlake[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
#endif
        }
    }else{
        fprintf(stderr, "ERROR: Wrong flag in function Model_Data::f_applyDYi\n");
        myexit(ERRCONSIS);
    }
} // END OF function Model_Data::f_applyDYi.


