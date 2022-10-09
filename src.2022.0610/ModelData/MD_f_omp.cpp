//  MD_f_omp.cpp
//
//  Created by Lele Shu on 9/5/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"

void Model_Data::f_applyDY_omp(double *DY, double t){
    double area;
    int isf, ius, igw, i;
#pragma omp parallel  default(shared) private(i) num_threads(CS.num_threads)
    {
#pragma omp for
        for (i = 0; i < NumEle; i++) {
            isf = iSF; ius = iUS; igw = iGW;
            area = Ele[i].area;
            QeleSurfTot[i] = Qe2r_Surf[i];
            QeleSubTot[i] = Qe2r_Sub[i];
            for (int j = 0; j < 3; j++) {
                QeleSurfTot[i] += QeleSurf[i][j];
                QeleSubTot[i] += QeleSub[i][j];
                //            CheckNANi(QeleSubTot[i], 1, "QeleSubTot[i]");
                //            CheckNANi(QeleSurfTot[i], 1, "QeleSurfTot[i]");
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
            DY[iUS] /= Ele[i].Sy;
            DY[iGW] /= Ele[i].Sy;
#ifdef DEBUG
            CheckNANi(DY[i], i, "DY[i] (Model_Data::f_applyDY)");
            CheckNANi(DY[iUS], i, "DY[iUS] (Model_Data::f_applyDY)");
            CheckNANi(DY[iGW], i, "DY[iGW] (Model_Data::f_applyDY)");
#endif
        }// end of for 1:NumEle
#pragma omp for
        for (i = 0; i < NumRiv; i++) {
            if(Riv[i].BC > 0){
                //            Newmann condition.
                DY[iRIV] = 0.;
            }else{
                DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].u_TopArea;
            }
            //        DY[iRIV] = 0.0;
#ifdef DEBUG
            CheckNANi(DY[i + 3 * NumEle], i, "DY[i] of river (Model_Data::f_applyDY)");
#endif
        } // end of for 1:NumRiv
    } // end of parallel
}

void Model_Data:: f_loop_omp( double  *Y, double  *DY, double t){
    int i;
    tnow = t;
#pragma omp parallel  default(shared) private(i) num_threads(CS.num_threads)
    {
#pragma omp for
        for (i = 0; i < NumEle; i++) {
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
            fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
        }
#pragma omp for
        for (i = 0; i < NumEle; i++) {
            /*========surf/gw flow Function==============*/
            fun_Ele_surface(i, t);  // AFTER infiltration, do the lateral flux. ESP for overland flow.
            fun_Ele_sub(i, t);
        } //end of for loop.
#pragma omp for
        for (i = 0; i < NumSegmt; i++) {
            fun_Seg_surface( RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
            fun_Seg_sub( RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        }
#pragma omp for
        for (i = 0; i < NumRiv; i++) {
            Flux_RiverDown(t, i);
        }
    }
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
}



void Model_Data::f_update_omp(double  *Y, double *DY, double t){
    /* Initialization of temporary state variables */
    int i;
#pragma omp parallel  default(shared) private(i) num_threads(CS.num_threads)
    {
#pragma omp for
        for (i = 0; i < NumY; i++) {
            DY[i] = 0.;
        }
        
#pragma omp for
        for (i = 0; i < NumEle; i++) {
            uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
            uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
            
            if(Ele[i].iBC == 0){ // NO BC
                uYgw[i] = max(0.0, Y[iGW]);
                Ele[i].QBC = 0.;
            }else if(Ele[i].iBC > 0){ // BC fix head
                Ele[i].yBC = tsd_eyBC.getX(t, Ele[i].iBC);
                uYgw[i] = Ele[i].yBC;
                Ele[i].QBC = 0.;
            }else{ // BC fix flux to GW
                Ele[i].QBC = tsd_eqBC.getX(t, -Ele[i].iBC);
            }
            /***** SS and BC *****/
//            for(int j = 0; j<3;j++){
//                Ele[i].iupdSF[j] = 0;
//                Ele[i].iupdGW[j] = 0;
//            }
//            Qe2r_Surf[i] = 0.;
//            Qe2r_Sub[i] = 0.;
            qEleExfil[i] = 0.;
            qEleInfil[i] = 0.;
            /********* Below are remove because the bass-balance issue. **********/
            //        for (int j = 0; j < 3; j++) {
            //            if(Ele[i].nabr[j] > 0){
            //                Ele[i].surfH[j] = (Ele[Ele[i].nabr[j] - 1].zmax + uYsf[Ele[i].nabr[j] - 1]);
            //            }else{
            //                Ele[i].surfH[j] = (Ele[i].zmax + uYsf[i]);
            //            }
            //        }
            //        Ele[i].dhBYdx = dhdx(Ele[i].surfX, Ele[i].surfY, Ele[i].surfH);
            //        Ele[i].dhBYdy = dhdy(Ele[i].surfX, Ele[i].surfY, Ele[i].surfH);
            //        Ele[i].Avg_Sf = sqpow2(Ele[i].dhBYdx, Ele[i].dhBYdy);
        }//end of for j=1:NumEle
        
#pragma omp for
        for (i = 0; i < NumRiv; i++ ){
            uYriv[i] = (Y[iRIV] >= 0.) ? Y[iRIV] : 0.;
            /* qrivsurf and qrivsub are calculated in Element fluxes.
             qrivDown and qrivUp are calculated in River fluxes. */
//            QrivDown[i] = 0.;
            Riv[i].updateRiver(uYriv[i]);
            /***** SS and BC *****/
            Riv[i].qBC = 0.0;
            if(Riv[i].BC == 0){
                /* Void */
            }else if(Riv[i].BC < 0){ // Fixed Flux INTO river Reaches.
                Riv[i].qBC = tsd_rqBC.getX(t, -Riv[i].BC);
            }else if (Riv[i].BC > 0){ // Fixed Stage of river reach.
                Riv[i].yBC = tsd_ryBC.getX(t, Riv[i].BC);
                uYriv[i] = Riv[i].yBC;
            }
        }
    } // end omp parallel.
}

