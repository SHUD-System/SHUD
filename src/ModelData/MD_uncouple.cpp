//  MD_uncouple.cpp
//
//  Created by Lele Shu on 9/6/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"

void Model_Data::f_loop1(double t){
    int i;
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
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
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
    /* Shared for both OpenMP and Serial, to update */
    //    PassValue();
}

void Model_Data:: f_loop3(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
        /*========ET Function==============*/
        //        f_etFlux(i, t); // et is moved out of f_loop()
        /*========infiltration/Recharge Function==============*/
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
        fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
    }
    for (i = 0; i < NumEle; i++) {
        fun_Ele_sub(i, t);
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
}

void Model_Data:: f_loop4(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
    }
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
void Model_Data:: f_loop5(double t){
}
void Model_Data::f_applyDYi(double *DY, double t, int flag){
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
            if(abs(DY[i])>1.5){
                CheckNANi(uYsf[i], i, "uYsf[i] of Surface (Model_Data::f_applyDYi)");
                CheckNANi(DY[i], i, "DY[i] of Surface (Model_Data::f_applyDYi)");
            }
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
            CheckNANi(uYus[i], i, "uYus[i] of Unsat (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of Unsat (Model_Data::f_applyDYi)");
        }
    }else if(flag ==3){
        for (int i = 0; i < NumEle; i++) {
            area = Ele[i].area;
            QeleSubTot[i] = Qe2r_Sub[i];
            for (int j = 0; j < 3; j++) {
                QeleSubTot[i] += QeleSub[i][j];
            }
            DY[i] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area;
            
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
//            DY[i] = 0;
            CheckNANi(uYgw[i], i, "uYgw[i] of Groundwater (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of Groundwater (Model_Data::f_applyDYi)");
        }
    }else if(flag ==4){
        for (int i = 0; i < NumRiv; i++) {
            if(Riv[i].BC > 0){
                //            Newmann condition.
                DY[i] = 0.;
            }else{
                DY[i] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].u_TopArea;
            }
            CheckNANi(uYriv[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
#ifdef DEBUG
#endif
        }
    }else if(flag ==5){
        for (int i = 0; i < NumLake; i++) {
            CheckNANi(uYlake[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
            CheckNANi(DY[i], i, "DY[i] of river (Model_Data::f_applyDYi)");
        }
    }else{
        fprintf(stderr, "ERROR: Wrong flag in function Model_Data::f_applyDYi\n");
        myexit(ERRCONSIS);
    }
} // END OF function Model_Data::f_applyDYi.


