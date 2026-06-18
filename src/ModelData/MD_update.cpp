#include "Model_Data.hpp"
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif

void Model_Data::f_updatei(double  *Y, double *DY, double t, int flag){
    switch (flag) {
        case 1:
            for (int i = 0; i < NumEle; i++) {
                uYsf[i] = (Y[i] >= 0.) ? Y[i] : 0.;
            }
            break;
        case 2:
            for (int i = 0; i < NumEle; i++) {
                uYus[i] = (Y[i] >= 0.) ? Y[i] : 0.;
            }
            break;
        case 3:
            for (int i = 0; i < NumEle; i++) {
                uYgw[i] = (Y[i] >= 0.) ? Y[i] : 0.;
                if(Ele[i].iBC == 0){ // NO BC
                    uYgw[i] = max(0.0, Y[i]);
                    Ele[i].QBC = 0.;
                }else if(Ele[i].iBC > 0){ // BC fix head
                    Ele[i].yBC = tsd_eyBC.getX(t, Ele[i].iBC);
                    uYgw[i] = Ele[i].yBC;
                    Ele[i].QBC = 0.;
                }else{ // BC fix flux to GW
                    Ele[i].QBC = tsd_eqBC.getX(t, -Ele[i].iBC);
                }
            }
            break;
        case 4:
            for (int i = 0; i < NumRiv; i++) {
                uYriv[i] = (Y[i] >= 0.) ? Y[i] : 0.;
//                uYriv[i] = Y[i];
//                QrivSurf[i] = 0.;
//                QrivSub[i] = 0.;
                QrivUp[i] = 0.;
                QrivDown[i] = 0.;
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
            break;
        case 5:
            for (int i = 0; i < NumLake; i++) {
                uYlake[i] = (Y[i] >= 0.) ? Y[i] : 0.;
            }
            break;
        default:
            break;
    }
}
void Model_Data::f_update(double  *Y, double *DY, double t){
    for (int i = 0; i < NumEle; i++) {
//        uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
//        uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
        for(int j = 0; j < 3; j++){
            QeleSub[i][j] = 0.;
            QeleSurf[i][j] = 0.;
            QeleSubTot[i] = 0.;
            QeleSurfTot[i] = 0.;
        }
        uYsf[i] = Y[iSF];
        uYus[i] = Y[iUS];
        if(Ele[i].iBC == 0){ // NO BC
//            uYgw[i] = max(0.0, Y[iGW]);
            uYgw[i] = Y[iGW];
            Ele[i].QBC = 0.;
        }else if(Ele[i].iBC > 0){ // BC fix head
            Ele[i].yBC = tsd_eyBC.getX(t, Ele[i].iBC);
            uYgw[i] = Ele[i].yBC;
            Ele[i].QBC = 0.;
        }else{ // BC fix flux to GW
            Ele[i].QBC = tsd_eqBC.getX(t, -Ele[i].iBC);
        }
        qEleExfil[i] = 0.;
        qEleInfil[i] = 0.;
        /***** SS and BC *****/
//        for(int j = 0; j<3;j++){
//            Ele[i].iupdSF[j] = 0;
//            Ele[i].iupdGW[j] = 0;
//        }
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
    
    for (int i = 0; i < NumRiv; i++ ){
        uYriv[i] = Y[iRIV];
        /* qrivsurf and qrivsub are calculated in Element fluxes.
         qrivDown and qrivUp are calculated in River fluxes. */
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
#ifdef DEBUG
        CheckNANi(uYriv[i], i, "uYriv in f_update.");
#endif
    }
    
    for (int i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    for (int i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }
    for (int i = 0; i < NumLake; i++) {
        yLakeStg[i] = Y[iLAKE];
        lake[i].yStage = yLakeStg[i];
        lake[i].update();
        y2LakeArea[i] = lake[i].u_toparea;
        QLakeSub[i] = 0.;
        QLakeSurf[i] = 0.;
        qLakeEvap[i] = 0.;
        qLakePrcp[i] = 0.;
        QLakeRivIn[i] = 0.;
        QLakeRivOut[i] = 0.;
    }
    for (int i = 0; i < NumY; i++){
        DY[i] = 0.;
    }
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_update", t, DY, NumY);
#endif
}
void Model_Data::summary (N_Vector udata){
    double  *Y;
    /* S1d.2 (openMP #48) — generic N_Vector data accessor (see f.cpp
     * + Macros.hpp comments). Replaces the prior backend-specific
     * NV_DATA_OMP / NV_DATA_S branch. */
    Y = N_VGetArrayPointer(udata);
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y[iSF];
        yEleUnsat[i] = Y[iUS];
        
        if(Ele[i].iBC > 0){
            yEleGW[i] = Ele[i].yBC;
        }else{
            yEleGW[i] = Y[iGW];
        }
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y[iRIV];
        //        uYriv[i] = Y[iRIV];
        if(Riv[i].BC > 0){
            yRivStg[i] = Riv[i].yBC;
        }else{
            yRivStg[i] = Y[iRIV];
        }
    }
}
void Model_Data::summary (N_Vector u1, N_Vector u2, N_Vector u3, N_Vector u4, N_Vector u5){

    /* S1d.2 (openMP #48) — generic N_VGetArrayPointer collapse of the
     * prior backend-split blocks (NV_Ith_OMP / NV_Ith_S). Hoist the
     * five pointer fetches out of the inner loops so the per-element
     * indexing stays a single load-and-store; bitwise-equivalent to
     * `NV_Ith_S(v, i)` because nvector_serial's NV_Ith_S expands to
     * `((NV_DATA_S(v))[i])` (verified against sundials 6.0.0). */
    double *Y1 = N_VGetArrayPointer(u1);
    double *Y2 = N_VGetArrayPointer(u2);
    double *Y3 = N_VGetArrayPointer(u3);
    double *Y4 = N_VGetArrayPointer(u4);
    double *Y5 = N_VGetArrayPointer(u5);
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y1[i];
        yEleUnsat[i] = Y2[i];
        yEleGW[i] = Y3[i];
        if(Ele[i].iBC > 0){
            yEleGW[i] = Ele[i].yBC;
        }
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y4[i];
        if(Riv[i].BC > 0){
            yRivStg[i] = Riv[i].yBC;
        }
    }
    for (int i = 0; i < NumLake; i++){
        yLakeStg[i] = Y5[i];
    }
    Sub2Global(yEleSurf, yEleUnsat, yEleGW, yRivStg, yLakeStg, NumEle, NumRiv, NumLake);
//    printVector(stdout, yEleSurf, 0, NumEle, 0);
//    printVector(stdout, yEleUnsat, 0, NumEle, 0);
//    printVector(stdout, yEleGW, 0, NumEle, 0);
//    printVector(stdout, yRivStg, 0, NumRiv, 0);
//
//    printVector(stdout, globalY, 0, NumEle, 0);
//    printVector(stdout, globalY, NumEle, NumEle, 0);
//    printVector(stdout, globalY, NumEle*2, NumEle, 0);
//    printVector(stdout, globalY, NumEle*3, NumRiv, 0);
}

int Model_Data::PrintInit (const char *fn, double t){
    unsigned long t_long = (long) t;
    if( t_long % CS.UpdateICStep ){
        return 0;
    }
    FILE           *fp;
    fp = fopen (fn, "w");
    CheckFile(fp, fn);
    /************* Element status **************/
    fprintf (fp, "%d\t %d \t%lf\n", NumEle, 6, t);
    fprintf (fp, "%s\t%s\t%s\t%s\t%s\t%s\n","Index",
             "Canopy", "Snow", "Surface", "Unsat", "GW");
    for (int i = 0; i < NumEle; i++){
        fprintf (fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", i+1, yEleIS[i], yEleSnow[i], yEleSurf[i], yEleUnsat[i], yEleGW[i]);
    }
    /************* River Reach status **************/
    fprintf (fp, "%d\t%d\n", NumRiv, 2);
    fprintf (fp, "%s\t%s\n", "Index", "Stage");
    for (int i = 0; i < NumRiv; i++){
        fprintf (fp, "%d\t%lf\n", i+1, yRivStg[i]);
    }
    /************* Lake status **************/
    if(NumLake > 0){
        fprintf (fp, "%d\t%d\n", NumLake, 2);
        fprintf (fp, "%s\t%s\n", "Index", "LakeStage");
        for (int i = 0; i < NumLake; i++){
            fprintf (fp, "%d\t%lf\n", i+1,yLakeStg[i]);
        }
    }
    fclose (fp);
    return 1;
}
