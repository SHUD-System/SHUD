#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ModelConfigure.hpp"
#include "IO.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"

void Model_Data::LoadIC(){
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
                yLakeStg[i] = 0.3 * (lake[i].bathymetry.yi[1] - lake[i].bathymetry.yi[0]);
            }
            break;
        default:
            /* reading from  */
            FILE * fp =  fopen(pf_in->file_init, "r");
            CheckFile(fp, pf_in->file_init);
            
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
                fprintf(stderr, "\nWARNING: number of RIVER rows (%d) in .ic DOSE NOT match the number of cell in .riv file (%d).\n Press anykey to continue ...\n", nr, NumRiv);
//                getchar();
            }
            for (int i = 0; i < NumRiv; i++) {
                yRivStg[i] = tb.x[i][1];
            }
            if (NumLake > 0){
                nr = tb.read(fp);
                if(nr != NumLake){
                    fprintf(stderr, "\nWARNING: number of LAKE rows (%d) in .ic DOSE NOT match the number of cell in .lake.bathy file (%d).\n Press anykey to continue ...\n", nr, NumEle);
//                getchar();
                    for (int i = 0; i < NumLake; i++) {
                        yLakeStg[i] = 2.;
                    }
                }else{
                    for (int i = 0; i < NumLake; i++) {
                        yLakeStg[i] = tb.x[i][1];
                    }
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
//    N_VectorContent_Serial x;
    for (int i = 0; i < NumLake; i++) {
//    N_VectorContent_Serial x;
//        x =(N_VectorContent_Serial)(udata->content);
//        printf("%f\n", x->data[i + 3 * NumEle + NumRiv]);
        SET_VALUE(udata, iLAKE) = globalY[iLAKE];
//        printf("%f\n", x->data[i + 3 * NumEle + NumRiv]);
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
void Model_Data::initializeLake(){
    if(lakeon){
        NumLake = LakeUniqueID();
        LakeInitialize();
        lake_readBathy(pf_in->file_lake_bathy);
        //        lake_readIC(fin->file_lake_ic);
        for(int i = 0; i < NumEle; i++){
            if(Ele[i].iLake > 0){
                if(Ele[i].z_surf < lake[Ele[i].iLake - 1].zmin){
                    printf("WARNING: lake(%d) zmin is higher than element(%d) zmax\n", Ele[i].iLake, i+1);
                }
            }
        }
    }
    
}
void Model_Data::initialize(){
    screeninfo("\nInitializing data structure ... \n");
    malloc_EleRiv();
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
            printf("%d\t", Ele[i].lakenabr[j]);
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
    flood = new FloodAlert();
    updateArea(); /* Calculate the Total Area of the watershed. */
    
    AccT_sub_max = gc.cfrozen.FT_sub_max;
    AccT_sub_min = gc.cfrozen.FT_sub_min;
    AccT_surf_max = gc.cfrozen.FT_surf_max;
    AccT_surf_min = gc.cfrozen.FT_surf_min;
    for (int i = 0; i < NumEle; i++) {
        AccT_sub[i].setLength(gc.cfrozen.FT_sub_Day);
        AccT_surf[i].setLength(gc.cfrozen.FT_surf_Day);
    }
    initializeLake();
    malloc_Y();
    read_cfgout(pf_in->file_cfgout);
}
void Model_Data:: initialize_output (){
    int ip = 0;
    /* Storage */
    if (CS.dt_ye_ic > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_y_ic, CS.dt_ye_ic,yEleIS, 0, io_ele);
    if (CS.dt_ye_snow > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_y_snow, CS.dt_ye_snow, yEleSnow, 0, io_ele);
    if (CS.dt_ye_surf > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_y_surf, CS.dt_ye_surf, yEleSurf, 0, io_ele);
    if (CS.dt_ye_unsat > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_y_unsat, CS.dt_ye_unsat, yEleUnsat, 0, io_ele);
    
    if (CS.dt_ye_gw > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_y_gw, CS.dt_ye_gw, yEleGW, 0, io_ele);
    /* Fluxes */
    if (CS.dt_qe_prcp > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_prcp, CS.dt_qe_prcp, qElePrep, 1, io_ele);
    if (CS.dt_qe_prcp > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_netprcp, CS.dt_qe_prcp, qEleNetPrep, 1, io_ele);
    if (CS.dt_qe_etp > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_ETP, CS.dt_qe_etp, qEleETP, 1, io_ele);
    if (CS.dt_qe_eta > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_ETA, CS.dt_qe_eta, qEleETA, 1, io_ele);
    if (CS.dt_qe_rech > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_rech, CS.dt_qe_rech, qEleRecharge, 1, io_ele);
    if (CS.dt_Qe_sub > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_Q_subTot, CS.dt_Qe_sub, QeleSubTot, 1, io_ele);
    }
    if (CS.dt_Qe_subx > 0){
        CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, pf_out->ele_Q_sub0, CS.dt_Qe_sub, QeleSub, 0, 1, io_ele);
        CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, pf_out->ele_Q_sub1, CS.dt_Qe_sub, QeleSub, 1, 1, io_ele);
        CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, pf_out->ele_Q_sub2, CS.dt_Qe_sub, QeleSub, 2, 1, io_ele);
    }
    if (CS.dt_Qe_surf > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_Q_surfTot, CS.dt_Qe_surf, QeleSurfTot, 1, io_ele);
    }
    if (CS.dt_Qe_surfx > 0){
        CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, pf_out->ele_Q_surf0, CS.dt_Qe_surf, QeleSurf, 0, 1, io_ele);
        CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, pf_out->ele_Q_surf1, CS.dt_Qe_surf, QeleSurf, 1, 1, io_ele);
        CS.PCtrl[ip++].InitIJ(ForcStartTime, NumEle, pf_out->ele_Q_surf2, CS.dt_Qe_surf, QeleSurf, 2, 1, io_ele);
    }
    if (CS.dt_Qe_rsub > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_Q_rsub, CS.dt_Qe_rsub, Qe2r_Sub, 1, io_ele);
    }
    if (CS.dt_Qe_rsurf > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_Q_rsurf, CS.dt_Qe_rsurf, Qe2r_Surf, 1, io_ele);
    }
    if (CS.dt_qe_infil > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_infil, CS.dt_qe_infil, qEleInfil, 1, io_ele);
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_exfil, CS.dt_qe_infil, qEleExfil, 1, io_ele);
    }
    
    /* ET  2d array */
    if (CS.dt_qe_et > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_ET[0], CS.dt_qe_et, qEleE_IC, 1, io_ele);
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_ET[1], CS.dt_qe_et, qEleTrans, 1, io_ele);
        CS.PCtrl[ip++].Init(ForcStartTime, NumEle, pf_out->ele_q_ET[2], CS.dt_qe_et, qEleEvapo, 1, io_ele);
    }
    
    
    /* Print control for River */
    if (CS.dt_Qr_up > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, pf_out->riv_Q_up, CS.dt_Qr_up, QrivUp, 1, io_riv);
    if (CS.dt_Qr_down > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, pf_out->riv_Q_down, CS.dt_Qr_down, QrivDown, 1, io_riv);
    if (CS.dt_Qr_sub > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, pf_out->riv_Q_sub, CS.dt_Qr_sub, QrivSub, 1, io_riv);
    if (CS.dt_Qr_surf > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, pf_out->riv_Q_surf, CS.dt_Qr_surf, QrivSurf, 1, io_riv);
    if (CS.dt_yr_stage > 0)
        CS.PCtrl[ip++].Init(ForcStartTime, NumRiv, pf_out->riv_y_stage, CS.dt_yr_stage, yRivStg, 0, io_riv);
    
    /* Print control for Lake */
    if (CS.dt_lake > 0 && NumLake > 0){
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_y_stage, CS.dt_lake, yLakeStg, 0, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_a_area, CS.dt_lake, y2LakeArea, 0, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_q_evap, CS.dt_lake, qLakeEvap, 1, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_q_prcp, CS.dt_lake, qLakePrcp, 1, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_Q_rivin, CS.dt_lake, QLakeRivIn, 1, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_Q_rivout, CS.dt_lake, QLakeRivOut, 1, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_Q_surf, CS.dt_lake, QLakeSurf, 1, io_lake);
        CS.PCtrl[ip++].Init(ForcStartTime, NumLake, pf_out->lake_Q_sub, CS.dt_lake, QLakeSub, 1, io_lake);
    }
    CS.NumPrint = ip;
    
    for (int i = 0; i < CS.NumPrint; i++)
    {
        CS.PCtrl[i].open_file(CS.Ascii, CS.Binary);
    }
}

