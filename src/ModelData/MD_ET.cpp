//
//  MD_ET.cpp
//  SHUD
//
//  Created by Lele Shu on 10/26/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
/* updateforcing() and ET() carry no profiling Timer of their own.
 * The `t_forcing_io` / `t_ET` buckets are timed by the RAII Timers
 * that wrap `MD->updateforcing(t)` and `MD->ET(t, tnext)` in the
 * shud.cpp main loop; a second Timer in here would add the same span
 * to the same bucket again and report more time than the wall clock. */
void Model_Data::updateforcing(double t){
    int i;
    /* Plain serial loop: updateforcing() is called from the
     * single-threaded main loop, outside any `omp parallel` region,
     * so a bare `#pragma omp for` here would have no effect. */
    for (i = 0; i < NumForc; i++){
        tsd_weather[i].movePointer(t);
    }
    tsd_MF.movePointer(t);
    tsd_LAI.movePointer(t);
//    tsd_RL.movePointer(t);
    for(i = 0; i < NumEle; i++){
        /* updateElement writes Ele[i].u_effKH / u_satn
         * / Kmax / u_deficit / u_theta / u_satKr / u_phius / u_effkInfi
         * (member method on AoS); sync_hot_dynamic(i) refreshes the SoA
         * mirror of the subset that the RHS reads afterwards. */
        Ele[i].updateElement(uYsf[i], uYus[i], uYgw[i]);
        sync_hot_dynamic(i);
        tReadForcing(t,i);
    }
}
void Model_Data::tReadForcing(double t, int i){
    /* Ele[i].{iForc, z_surf, iLC, iMF, Albedo, FixPressure, iLake,
     * windH} are read from the SoA mirror. tsd_weather / tsd_LAI /
     * tsd_MF are not per-element state and are accessed directly. */
    int idx = hot.iForc[i] - 1;
    double etp, ra, rs, t0, hc, U2, Uz, Zmeasure, lai;
    double GroundHeatFlux, RG;
    t_prcp[i] = tsd_weather[idx].getX(t, i_prcp) * gc.cPrep;
    t0= tsd_weather[idx].getX(t, i_temp);
    t_temp[i] = TemperatureOnElevation(t0, hot.z_surf[i], tsd_weather[idx].xyz[2]) +  gc.cTemp;
    t_lai[i] = tsd_LAI.getX(t, hot.iLC[i]) * gc.cLAItsd ;
    lai = t_lai[i];
    t_mf[i] = tsd_MF.getX(t, hot.iMF[i]) * gc.cMF / 1440.;  /*  [m/day/C] to [m/min/C].
                                                            1.6 ~ 6.0 mm/day/C is typical value in USDA book
                                                            Input is 1.4 ~ 3.0 mm/d/c */
    t_rn[i] = tsd_weather[idx].getX(t, i_rn) * (1 - hot.Albedo[i]);
    Uz = t_wind[i] = (fabs(tsd_weather[idx].getX(t, i_wind) ) + 0.001); // +.001 voids ZERO.
    t_rh[i] = tsd_weather[idx].getX(t, i_rh);
//    t_hc[i] = tsd_RL.getX(t, hot.iLC[i]);
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
    double Gamma = PsychrometricConstant(hot.FixPressure[i], lambda); // eq 4.2.28  [kPa C-1]
    double es = VaporPressure_Sat(t_temp[i]);                   // eq 4.2.2 [kpa]
    double ea = es * t_rh[i];   // [kPa]
    double ed = es - ea ;  // [kPa]
    double Delta = SlopeSatVaporPressure(t_temp[i], es);        // eq 4.2.3 [kPa C-1]
    double rho = AirDensity(hot.FixPressure[i], t_temp[i]);;    // eq 4.2.4 [kg m-3]
    /* R - G in the PM equation.*/
    if(hot.iLake[i] > 0 ){
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
    U2 = WindProfile(2.0, t_wind[i], hot.windH[i], 0., ROUGHNESS_WATER); // [m s-1]
    qPotEvap[i] = gc.cETP * PET_PM_openwater(Delta, Gamma, lambda, RG, ed, U2) * 60.; // eq 4.2.30
    if(hot.iLake[i] > 0){        /* Open-water */
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
        etp = qPotTran[i] * hot.VegFrac[i] + qPotEvap[i] * (1. - hot.VegFrac[i]);
        CheckNANi(qPotTran[i], i, "qPotTran[i]");
    }
    qEleETP[i] = etp;
}
void Model_Data::ET(double t, double tnext){
    /* Not timed here; the caller in shud.cpp owns the t_ET Timer (see
     * the comment above updateforcing()). */
    double  DT_min = tnext - t;
    /* Plain serial loop, called outside any `omp parallel` region.
     * All element-local scalars (T, LAI, MF, prcp, snFrac, snAcc,
     * snMelt, snStg, icAcc, icEvap, icStg, icMax, vgFrac, ta_surf,
     * ta_sub, and the loop index i) are declared at use inside the
     * for body, so each iteration touches only its own element.
     * DT_min is the only shared value and is loop-invariant. */
    for(int i = 0; i < NumEle; i++) {
        double T = t_temp[i];
        double prcp = t_prcp[i];
        /* Snow Accumulation */
        double MF = t_mf[i];
        double snStg = yEleSnow[i];
        /* Snow Accumulation/Melt Calculation*/
        double snFrac = FrozenFraction(T, Train, Tsnow);

        if(CS.cryosphere){
            AccT_surf[i].push(T, t);
            AccT_sub[i].push(T, t);
            double ta_surf = AccT_surf[i].getACC();
            double ta_sub  = AccT_sub[i].getACC();
            fu_Sub[i] = 1. - FrozenFraction(ta_sub, AccT_sub_max, AccT_sub_min);
            fu_Surf[i] = 1. - FrozenFraction(ta_surf, AccT_surf_max, AccT_surf_min);
        }else{
            fu_Sub[i] = 1.;
            fu_Surf[i] = 1.;
        }

        double snAcc = snFrac * prcp;
        double snMelt = (T > To ? (T - To) * MF : 0.);    /* eq. 7.3.14 in Maidment */
        snMelt = min(max(0., snStg / DT_min), max(0., snMelt));
//        CheckNonNegative(snMelt, i, "Snow Melting");
        snStg += (snAcc - snMelt) * DT_min;

        /* Interception */
        double LAI = t_lai[i];
        double icStg = yEleIS[i];
        /* Ele[i].VegFrac SoA read. */
        double vgFrac = hot.VegFrac[i];
        double icAcc, icEvap;
        if(LAI > ZERO){
            double icMax = gc.cISmax * IC_MAX * LAI;
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
    /* Ele[i].{VegFrac, ImpAF, iSoil, u_satn, WetlandLevel,
     * RootReachLevel} are read from the SoA mirror. The Soil[idx]
     * lookup uses the Soil array itself (not per-element state); only
     * the `iSoil - 1` selector comes from the SoA. */
    double Es = 0., Eu = 0., Tu = 0., Eg = 0., Tg = 0.;
    double va = hot.VegFrac[i], vb = 1. - hot.VegFrac[i];
    double pj = 1. - hot.ImpAF[i];
    iBeta[i] = SoilMoistureStress(Soil[(hot.iSoil[i] - 1)].ThetaS, Soil[(hot.iSoil[i] - 1)].ThetaR, hot.u_satn[i]);
    /* Evaporation from SURFACE ponding water */
    Es = min(max(0., uYsf[i]), qPotEvap[i]) * vb;
    if(Es < qPotEvap[i]){
        /* Some PET is extracted by surface Evaporation, so PET - Es is the effective PET now. */
        if(uYgw[i] > hot.WetlandLevel[i]){
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
            if(uYgw[i] > hot.RootReachLevel[i]){
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
    /* The per-element AET/PET warning is compiled only under
     * `#ifdef DEBUG`: f_etFlux runs once per element inside rhs_flux,
     * and an unconditional `printf` here would put a stdout write on
     * the RHS hot path. Default builds (no -DDEBUG) emit no code
     * here. The warning only ever goes to stdout, never to the model
     * output files, so results are the same with or without DEBUG.
     * This matches how `CheckNANi` is guarded in MD_rhs_core.cpp. */
#ifdef DEBUG
    if(qEleETA[i] > qEleETP[i] * 2.){
        printf("Warning: More AET(%.3E) than PET(%.3E) on Element (%d).", qEleETA[i], qEleETP[i], i+1);
    }
#endif
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

