#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ModelConfigure.hpp"
#include "IO.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"

void Model_Data::LoadIC(FileIn *fid){
    for (int i = 0; i < NumEle; i++) {
        yEleWetFront[i] = 0.;
    }
    switch (CS.init_type) {
        case 0:
            /* groundwater relief case */
            for (int i = 0; i < NumEle; i++) {
                yEleIS[i] = 0.;
                yEleSnow[i] = 0.;
                yEleSurf[i] = 0.;
                yEleUnsat[i] = 0.;
                yEleGW[i] = Ele[i].AquiferDepth;
            }
            for (int i = 0; i < NumRiv; i++) {
                yRivStg[i] = 0.;
            }
            for (int i = 0; i < NumLake; i++) {
                yLakeStg[i] = 0.;
            }
            break;
        case 1:
            /* Default case */
            for (int i = 0; i < NumEle; i++) {
                yEleIS[i] = 0.;
                yEleSnow[i] = 0.;
                yEleSurf[i] = 0.;
                yEleUnsat[i] = 0.;
                yEleGW[i] = 0.;
            }
            for (int i = 0; i < NumRiv; i++) {
                yRivStg[i] = 0.;
            }
            for (int i = 0; i < NumLake; i++) {
                yLakeStg[i] = 0.;
            }
            break;
        case 2:
            /* default guess case */
            for (int i = 0; i < NumEle; i++) {
                yEleIS[i] = 0.;
                yEleSnow[i] = 0.;
                yEleSurf[i] = 0.;
                yEleUnsat[i] = 0.3 * Ele[i].AquiferDepth;
                yEleGW[i] = 0.4 * Ele[i].AquiferDepth;
            }
            for (int i = 0; i < NumRiv; i++) {
                yRivStg[i] = 0.2 * Riv[i].depth;
            }
            for (int i = 0; i < NumLake; i++) {
                yLakeStg[i] = 0.;
            }
            break;
        default:
            /* reading from  */
            FILE * fp =  fopen(fid->file_init, "r");
            CheckFile(fp, fid->file_init);
            
            TabularData tb;
            int nr = tb.read(fp);
            if( nr != NumEle){
                fprintf(stderr, "\nWARNING: number of rows (%d) in .ic DOSE NOT match the number of cell in .mesh file (%d).\n Press anykey to continue ...\n", nr, NumEle);
//                getchar();
            }
            for (int i = 0; i < NumEle; i++) {
                yEleIS[i]   = tb.x[i][1];
                yEleSnow[i] = tb.x[i][2];
                yEleSurf[i] = tb.x[i][3];
                yEleUnsat[i] = tb.x[i][4];
                yEleGW[i]   = tb.x[i][5];
            }
            
            nr = tb.read(fp);
            if( nr != NumRiv){
                fprintf(stderr, "\nWARNING: number of rows (%d) in .ic DOSE NOT match the number of cell in .riv file (%d).\n Press anykey to continue ...\n", nr, NumRiv);
//                getchar();
            }
            for (int i = 0; i < NumRiv; i++) {
                yRivStg[i] = tb.x[i][1];
            }
            if (NumLake > 0){
                nr = tb.read(fp);
                if(nr != NumLake){
                    fprintf(stderr, "\nWARNING: number of rows (%d) in .ic DOSE NOT match the number of cell in .riv file (%d).\n Press anykey to continue ...\n", nr, NumEle);
//                    getchar();
                }
                for (int i = 0; i < NumLake; i++) {
                    yLakeStg[i] = tb.x[i][1];
                }
            }
            fclose(fp);
            break;
    } /* End of Switch*/
    
    for (int i = 0; i < NumEle; i++) {
        yEleSnowGrnd[i] = (1 - Ele[i].VegFrac) * yEleSnow[i];
        yEleSnowCanopy[i] = Ele[i].VegFrac * yEleSnow[i];
    }
    Sub2Global(yEleSurf, yEleUnsat, yEleGW, yRivStg, yLakeStg, NumEle, NumRiv, NumLake);
}
void Model_Data::SetIC2Y(N_Vector udata){
    /* PUT the values into CV_Y */
    for (int i = 0; i < NumEle; i++) {
        SET_VALUE(udata, iSF) = globalY[iSF];
        SET_VALUE(udata, iUS) = globalY[iUS];
        SET_VALUE(udata, iGW) = globalY[iGW];
    }
    for (int i = 0; i < NumRiv; i++) {
        SET_VALUE(udata, iRIV) = globalY[iRIV];
    }
    for (int i = 0; i < NumLake; i++) {
        SET_VALUE(udata, iLAKE) = globalY[iLAKE];
    }
}
void Model_Data::SetIC2Y(N_Vector udata1, N_Vector udata2,
                         N_Vector udata3, N_Vector udata4,
                         N_Vector udata5){
    /* PUT the values into CV_Y */
    for (int i = 0; i < NumEle; i++) {
        SET_VALUE(udata1, i) = globalY[iSF];
        SET_VALUE(udata2, i) = globalY[iUS];
        SET_VALUE(udata3, i) = globalY[iGW];
    }
    for (int i = 0; i < NumRiv; i++) {
        SET_VALUE(udata4, i) = globalY[iRIV];
    }
    for (int i = 0; i < NumLake; i++) {
        SET_VALUE(udata5, i) = globalY[iLAKE];
    }
}

void Model_Data:: initialize(){
    screeninfo("\nInitializing data structure ... \n");
    allocateMemory();
    copyCalib(); /* Put the calibration parameter into classes*/
    
    for (int i = 0; i < NumNode; i++) {
        /* Update the aquifer depth. */
        Node[i].Init(gc.cAqD);
    }
    for (int i = 0; i < NumEle; i++) {
        Ele[i].applyGeometry(Node);
        Ele[i].copyGeol(Geol);
        Ele[i].copySoil(Soil);
        Ele[i].copyLandc(LandC);
        Ele[i].windH = HeightWindMeasure;
        Ele[i].InitElement();
        Ele[i].infKsatV *= 1 - Ele[i].Landcover::SoilDgrd;
        Ele[i].macKsatV *= 1 - Ele[i].Landcover::SoilDgrd;
        Ele[i].VegFrac *= 1 - Ele[i].Landcover::ImpAF;
    }
    for (int i = 0; i < NumSegmt; i++){
        Ele[RivSeg[i].iEle - 1].RivID = RivSeg[i].iRiv;
        Ele[RivSeg[i].iEle - 1].RivSegID = i + 1;
    }
    rmSinks();
    for (int i = 0; i < NumEle; i++) {
        Ele[i].applyNabor(Node, Ele);
    }
#ifdef DEBUG
    /* Topological relationship between Elements*/
    for (int i = 0; i < NumEle; i++) {
        printf("%d: ", i+1);
        for(int j = 0; j < 3; j++){
            printf("%d\t", Ele[i].nabr[j]);
        }
        for(int j = 0; j < 3; j++){
            printf("%d\t", Ele[i].nabrToMe[j]);
        }
        printf("\n");
    }
#endif
//    rmSinks();
    /* Correct river bed elevation */
//    CorrectRiver(0.05);
    for (int i = 0; i < NumRiv; i++) {
        // loop 1 for river
        Riv[i].initialRiver(Riv_Type);
        Riv[i].BedSlope = max(MINRIVSLOPE, Riv[i].BedSlope);
    }
    for (int i = 0; i < NumRiv; i++){
        // loop 2 for river. update the values for downstream.
        Riv[i].updateFrDownstream(Riv);
    }
    for (int i = 0; i < NumSegmt; i++){
        RivSeg[i].Cwr = Riv_Type[Riv[RivSeg[i].iRiv - 1].type - 1 ].Cwr;
        RivSeg[i].KsatH = Riv_Type[Riv[RivSeg[i].iRiv - 1].type - 1 ].KsatH;
        RivSeg[i].eqDistance = Ele[RivSeg[i].iEle - 1].area / RivSeg[i].length * .5;
        CheckNonZero(RivSeg[i].Cwr, i, "River Segment Cwr");
//        CheckNonZero(RivSeg[i].KsatH, i, "River Segment KsatH");
    }
#ifndef _CALIBMODE
    flood = new FloodAlert();
#endif
    updateArea(); /* Calculate the Total Area of the watershed. */
}
void Model_Data:: initialize_output (FileOut *fout){
    int ip = 0;
    /* Storage */
    if (CS.dt_ye_ic > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_y_ic, CS.dt_ye_ic,yEleIS, 0);
    if (CS.dt_ye_snow > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_y_snow, CS.dt_ye_snow, yEleSnow, 0);
    if (CS.dt_ye_surf > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_y_surf, CS.dt_ye_surf, yEleSurf, 0);
    if (CS.dt_ye_unsat > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_y_unsat, CS.dt_ye_unsat, yEleUnsat, 0);
    
    if (CS.dt_ye_gw > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_y_gw, CS.dt_ye_gw, yEleGW, 0);
    /* Fluxes */
    if (CS.dt_qe_prcp > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_prcp, CS.dt_qe_prcp, qElePrep, 1);
    if (CS.dt_qe_prcp > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_netprcp, CS.dt_qe_prcp, qEleNetPrep, 1);
    if (CS.dt_qe_etp > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_ETP, CS.dt_qe_etp, qEleETP, 1);
    if (CS.dt_qe_eta > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_ETA, CS.dt_qe_eta, qEleETA, 1);
    if (CS.dt_qe_rech > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_rech, CS.dt_qe_rech, qEleRecharge, 1);
    if (CS.dt_Qe_sub > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_Q_subTot, CS.dt_Qe_sub, QeleSubTot, 1);
    }
    if (CS.dt_Qe_surf > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_Q_surfTot, CS.dt_Qe_surf, QeleSurfTot, 1);
    }
    if (CS.dt_Qe_rsub > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_Q_rsub, CS.dt_Qe_rsub, Qe2r_Sub, 1);
    }
    if (CS.dt_Qe_rsurf > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_Q_rsurf, CS.dt_Qe_rsurf, Qe2r_Surf, 1);
    }
    if (CS.dt_qe_infil > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_infil, CS.dt_qe_infil, qEleInfil, 1);
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, fout->ele_q_exfil, CS.dt_qe_infil, qEleExfil, 1);
    }
    
    /* ET  2d array */
    for (int i = 0; i < 3; i++){
        if (CS.dt_qe_et[i] > 0)
            CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, fout->ele_q_ET[i], CS.dt_qe_et[i], qEleET, i, 1);
    }
    /* Print control for River */
    if (CS.dt_Qr_up > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, fout->riv_Q_up, CS.dt_Qr_up, QrivUp, 1);
    if (CS.dt_Qr_down > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, fout->riv_Q_down, CS.dt_Qr_down, QrivDown, 1);
    if (CS.dt_Qr_sub > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, fout->riv_Q_sub, CS.dt_Qr_sub, QrivSub, 1);
    if (CS.dt_Qr_surf > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, fout->riv_Q_surf, CS.dt_Qr_surf, QrivSurf, 1);
    if (CS.dt_yr_stage > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, fout->riv_y_stage, CS.dt_yr_stage, yRivStg, 0);
    
    /* Print control for Lake */
    if (CS.dt_yl_stage > 0 || NumLake > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, fout->lake_y_stage, CS.dt_yl_stage, yLakeStg, 0);
    if (CS.dt_ql_et > 0 || NumLake > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, fout->lake_q_evap, CS.dt_ql_et, qLakeEvap, 1);
    if (CS.dt_ql_prcp > 0 || NumLake > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, fout->lake_q_prcp, CS.dt_ql_prcp, qLakePrcp, 1);
    if (CS.dt_Ql_chn > 0 || NumLake > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, fout->lake_Q_riv, CS.dt_Ql_chn, QLakeRiv, 1);
    if (CS.dt_Ql_surf > 0 || NumLake > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, fout->lake_Q_surf, CS.dt_Ql_surf, QLakeSurf, 1);
    if (CS.dt_Ql_sub > 0 || NumLake > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, fout->lake_Q_sub, CS.dt_Ql_sub, QLakeSub, 1);
    
    CS.NumPrint = ip;
    
    for (int i = 0; i < CS.NumPrint; i++)
    {
        CS.PCtrl[i].open_file(CS.Ascii, CS.Binary);
    }
}

