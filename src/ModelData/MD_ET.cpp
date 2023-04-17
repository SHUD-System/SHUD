//
//  MD_ET.cpp
//  SHUD
//
//  Created by Lele Shu on 10/26/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
void Model_Data::updateforcing(double t){
    int i;
#ifdef _OPENMP_ON
#pragma omp for
#endif
    for (i = 0; i < NumForc; i++){
        tsd_weather[i].movePointer(t);
    }
    tsd_MF.movePointer(t);
    tsd_LAI.movePointer(t);
//    tsd_RL.movePointer(t);
    for(i = 0; i < NumEle; i++){
        Ele[i].updateElement(uYsf[i], uYus[i], uYgw[i]);
        tReadForcing(t,i);
    }
}
void Model_Data::tReadForcing(double t, int i){
    int idx = Ele[i].iForc - 1;
    double etp, ra, rs, t0, hc, U2, Uz, Zmeasure, lai;
    double GroundHeatFlux, RG;
    t_prcp[i] = tsd_weather[idx].getX(t, i_prcp) * gc.cPrep;
    t0= tsd_weather[idx].getX(t, i_temp);
    t_temp[i] = TemperatureOnElevation(t0, Ele[i].z_surf, tsd_weather[idx].xyz[2]) +  gc.cTemp;
    t_lai[i] = tsd_LAI.getX(t, Ele[i].iLC) * gc.cLAItsd ;
    lai = t_lai[i];
    t_mf[i] = tsd_MF.getX(t, Ele[i].iMF) * gc.cMF / 1440.;  /*  [m/day/C] to [m/min/C].
                                                            1.6 ~ 6.0 mm/day/C is typical value in USDA book
                                                            Input is 1.4 ~ 3.0 mm/d/c */
    t_rn[i] = tsd_weather[idx].getX(t, i_rn) * (1 - Ele[i].Albedo);
    Uz = t_wind[i] = (fabs(tsd_weather[idx].getX(t, i_wind) ) + 0.001); // +.001 voids ZERO.
    t_rh[i] = tsd_weather[idx].getX(t, i_rh);
//    t_hc[i] = tsd_RL.getX(t, Ele[i].iLC);
//    t_hc[i] = max(t_hc[i], CONSt_hc);
    /* Precipitation  */
    t_prcp[i]   = t_prcp[i] * 0.001 / 1440. ; // [mm d-1] to [m min-1]
    /* Potential ET */
    t_rn[i]     = t_rn[i] * 1.0e-6;  // [W m-2] to [MJ m-2 s-1]
    /*
     t_wind [m s-1] ;
     t_temp  [C] ;
     t_rh  [0-1] ;
     */
    t_rh[i]     = min(max(t_rh[i], CONST_RH), 1.0); // [value is b/w 0~1 ]
    
    qElePrep[i] = t_prcp[i];
    double lambda = LatentHeat(t_temp[i]);                      // eq 4.2.1  [MJ/kg]
    double Gamma = PsychrometricConstant(Ele[i].FixPressure, lambda); // eq 4.2.28  [kPa C-1]
    double es = VaporPressure_Sat(t_temp[i]);                   // eq 4.2.2 [kpa]
    double ea = es * t_rh[i];   // [kPa]
    double ed = es - ea ;  // [kPa]
    double Delta = SlopeSatVaporPressure(t_temp[i], es);        // eq 4.2.3 [kPa C-1]
    double rho = AirDensity(Ele[i].FixPressure, t_temp[i]);;    // eq 4.2.4 [kg m-3]
    /* R - G in the PM equation.*/
    if(Ele[i].iLake > 0 ){
        GroundHeatFlux = 0.;
        RG = t_rn[i];
    }else{
        if(lai > 0){
            GroundHeatFlux = 0.4 * exp(-0.5 * lai) * t_rn[i];
        }else{
            GroundHeatFlux = 0.1 * t_rn[i];
        }
    }
    RG = t_rn[i] - GroundHeatFlux;
    U2 = WindProfile(2.0, t_wind[i], Ele[i].windH, 0., ROUGHNESS_WATER); // [m s-1]
    qPotEvap[i] = gc.cETP * PET_PM_openwater(Delta, Gamma, lambda, RG, ed, U2) * 60.; // eq 4.2.30
    if(Ele[i].iLake > 0){        /* Open-water */
        qPotTran[i] = gc.cETP * 0.;
        etp = qPotEvap[i];
    }else if(lai <= 0.){        /* Bare soiln */
        qPotTran[i] = gc.cETP * 0.;
        etp = qPotEvap[i];
    }else{
//        hc = lai2hc(lai);
        hc = lai * 0.5;
        Zmeasure = hc*1.3333; /* When hc > Zm, Zm = hc + 5.0m */
        ra = AerodynamicResistance(Uz, hc, Zmeasure, Zmeasure); // eq 4.2.25  [s m-1]
//        if( Zmeasure > HeightWindMeasure){ /* Veg Height > Wind Measure Height */
//            ra = AerodynamicResistance(Uz, hc, Zmeasure, Zmeasure); // eq 4.2.25  [s m-1]
//        }else{
//            ra = AerodynamicResistance(Uz, hc, HeightWindMeasure, 2.0); // eq 4.2.25  [s m-1]
//        }
//        ra = min(300., ra);
//        if(ra < 0){
//            ra = AerodynamicResistance(Uz, hc, Zmeasure, Zmeasure); // eq 4.2.25  [s m-1]
//            ra = AerodynamicResistance(Uz, hc, HeightWindMeasure, 2.0); // eq 4.2.25  [s m-1]
//        }
        CheckNonZero(ra, i, "Aerodynamic Resistance");
//        CheckNANi(ra, i, "Aerodynamic Resistance");
        rs = BulkSurfaceResistance(lai);  // eq 4.2.22 & 4.2.25  [s m-1]
        qPotTran[i] = gc.cETP * PET_Penman_Monteith(RG, rho, ed, Delta, ra, rs, Gamma, lambda) * 60.;// eq 4.2.27
        etp = qPotTran[i] * Ele[i].VegFrac + qPotEvap[i] * (1. - Ele[i].VegFrac);
        CheckNANi(qPotTran[i], i, "qPotTran[i]");
    }
    qEleETP[i] = etp;
}
void Model_Data::ET(double t, double tnext){
    double  T=NA_VALUE,  LAI=NA_VALUE, MF =NA_VALUE, prcp = NA_VALUE;
    double  snFrac, snAcc, snMelt, snStg;
    double  icAcc, icEvap, icStg, icMax, vgFrac;
    double  DT_min = tnext - t;
    double  ta_surf, ta_sub;
    int i;
#ifdef _OPENMP_ON
#pragma omp for
#endif
    for(i = 0; i < NumEle; i++) {
        T = t_temp[i];
        prcp = t_prcp[i];
        /* Snow Accumulation */
        MF = t_mf[i];
        snStg = yEleSnow[i];
        /* Snow Accumulation/Melt Calculation*/
        snFrac  = FrozenFraction(T, Train, Tsnow);
        
        if(CS.cryosphere){
            AccT_surf[i].push(T, t);
            AccT_sub[i].push(T, t);
            ta_surf = AccT_surf[i].getACC();
            ta_sub  = AccT_sub[i].getACC();
            fu_Sub[i] = 1. - FrozenFraction(ta_sub, AccT_sub_max, AccT_sub_min);
            fu_Surf[i] = 1. - FrozenFraction(ta_surf, AccT_surf_max, AccT_surf_min);
        }else{
            fu_Sub[i] = 1.;
            fu_Surf[i] = 1.;
        }
        
        snAcc = snFrac * prcp;
        snMelt = (T > To ? (T - To) * MF : 0.);    /* eq. 7.3.14 in Maidment */
        snMelt = min(max(0., snStg / DT_min), max(0., snMelt));
//        CheckNonNegative(snMelt, i, "Snow Melting");
        snStg += (snAcc - snMelt) * DT_min;
        
        /* Interception */
        LAI = t_lai[i];
        icStg = yEleIS[i];
        vgFrac = Ele[i].VegFrac;
        if(LAI > ZERO){
            icMax = gc.cISmax * IC_MAX * LAI;
            icAcc = min(prcp - snAcc, max(0., (icMax - icStg) / DT_min) );
            icEvap = min(max(0., icStg / DT_min), qPotEvap[i]);
        }else{
            icAcc = 0.;
            icEvap = 0.;
        }
        icStg += (icAcc - icEvap) * DT_min;
        
        /* Update the storage value and net precipitaion */
        yEleIS[i] = icStg * vgFrac;
        yEleSnow[i] = snStg;
        qEleE_IC[i] = icEvap * vgFrac;
        qEleNetPrep[i] = (1. - snFrac) * prcp + snMelt - icAcc * vgFrac ;

//        CheckNonNegative(qEleNetPrep[i], i, "Net Precipitation");
//        CheckNonNegative(qEleE_IC[i], i, "qEleE_IC");
    }
}
void Model_Data::f_etFlux(int i, double t){
    double Es = 0., Eu = 0., Tu = 0., Eg = 0., Tg = 0.;
    double va = Ele[i].VegFrac, vb = 1. - Ele[i].VegFrac;
    double pj = 1. - Ele[i].ImpAF;
    iBeta[i] = SoilMoistureStress(Soil[(Ele[i].iSoil - 1)].ThetaS, Soil[(Ele[i].iSoil - 1)].ThetaR, Ele[i].u_satn);
    /* Evaporation from SURFACE ponding water */
    Es = min(max(0., uYsf[i]), qPotEvap[i]) * vb;
    if(Es < qPotEvap[i]){
        /* Some PET is extracted by surface Evaporation, so PET - Es is the effective PET now. */
        if(uYgw[i] > Ele[i].WetlandLevel){
            /* Evporation from GroundWater, ONLY when gw above wetland level*/
            Eg = min(max(0., uYgw[i]), qPotEvap[i] - Es) * pj * vb;
            Eu = 0.;
        }else{
            Eg = 0.;
            /* Evaporation from Unsaturated Zone. */
            Eu = min(max(0., uYus[i]), iBeta[i] * (qPotEvap[i] - Es)) * pj * vb;
        }
    }else{
        /* All evporation is from land surface ONLY */
        Eg = 0.;
        Eu = 0.;
    }
    /* Vegetation Transpiration */
    if(t_lai[i] > ZERO){
        if(qEleE_IC[i] >= qPotTran[i]){
            Tg = Tu = 0.;
            qEleE_IC[i] = qPotTran[i] * pj * va;
        }else{
            if(uYgw[i] > Ele[i].RootReachLevel){
                Tg = min(max(0., uYgw[i]), (qPotTran[i] - qEleE_IC[i]) ) * pj * va;
                Tu = 0.;
            }else{
                Tg = 0.;
                Tu = min(max(0., uYus[i]), iBeta[i] * (qPotTran[i] - qEleE_IC[i]) ) * pj * va;
            }
        }
    }else{
        Tg = Tu = qEleE_IC[i] = 0.;
    }
    qEs[i] = Es;
    qEu[i] = Eu;
    qEg[i] = Eg;
    qTu[i] = Tu;
    qTg[i] = Tg;
    qEleTrans[i] = Tg + Tu;
    qEleEvapo[i] = Eu + Eg + Es;  
    qEleETA[i] = qEleE_IC[i] + qEleEvapo[i] + qEleTrans[i];
    if(qEleETA[i] > qEleETP[i] * 2.){
        printf("Warning: More AET(%.3E) than PET(%.3E) on Element (%d).", qEleETA[i], qEleETP[i], i+1);
    }
    CheckNonNegative(Es, i, "Es"); // Debug Only
    CheckNonNegative(Eu, i, "Eu");
    CheckNonNegative(Eg, i, "Eg");
    CheckNonNegative(Tu, i, "Tu");
    CheckNonNegative(Tg, i, "Tg");
    CheckNANi(qEleETA[i], i, "Potential ET (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleEvapo[i], i, "Transpiration (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleTrans[i], i, "Soil Evaporation (Model_Data::EvapoTranspiration)");
#ifdef DEBUG
#endif
}

