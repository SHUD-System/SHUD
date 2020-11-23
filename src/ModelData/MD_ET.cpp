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
    tsd_RL.movePointer(t);
    for(i = 0; i < NumEle; i++){
        Ele[i].updateElement(uYsf[i], uYus[i], uYgw[i]);
        tReadForcing(t,i);
    }
}
void Model_Data::tReadForcing(double t, int i){
    t_prcp[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_prcp) * gc.cPrep;
    t_temp[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_temp) +  gc.cTemp;
    t_lai[i] = tsd_LAI.getX(t, Ele[i].iLC) * gc.cLAItsd ;
    t_mf[i] = tsd_MF.getX(t, Ele[i].iMF) * gc.cMF;
    t_rn[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_rn) * (1 - Ele[i].Albedo);
    t_wind[i] = (fabs(tsd_weather[Ele[i].iForc - 1].getX(t, i_wind) ) + 1); // +1 voids ZERO.
    t_rh[i] = tsd_weather[Ele[i].iForc - 1].getX(t, i_rh);
    t_rl[i] = tsd_RL.getX(t, Ele[i].iLC);
    t_rl[i] = max(t_rl[i], CONST_RL);
    t_rh[i] = min(max(t_rh[i], CONST_RH), 1.0); // [-]
    t_prcp[i]   = t_prcp[i] / 1440.; // [m min-1]
    t_rn[i]     = t_rn[i] / 1440.;  // J/m2/day to [MJ m-2 min-1]
    t_wind[i]   = t_wind[i] / 1440. ; // m/d =>> m/min  [m min-1]
    
    qElePrep[i] = t_prcp[i];
    
    double lambda = LatentHeat(t_temp[i]);                      // eq 4.2.1  [MJ/kg]
    double es = VaporPressure_Sat(t_temp[i]);                   // eq 4.2.2 [kpa]
    double ea = es * t_rh[i];   // [kPa]
    double ed = es - ea ;  // [kPa]
    double Delta = SlopeSatVaporPressure(t_temp[i], es);        // eq 4.2.3 [kPa C-1]
    double rho = AirDensity(Ele[i].FixPressure, t_temp[i]);;    // eq 4.2.4 [kg m-3]
    double r_s = 0.;                                            // (bulk) surface resistances. [min m-1]
    if(t_lai[i] > 0.){
        r_s = BulkSurfaceResistance(t_lai[i]);                  // eq 4.2.22 [min m-1]
    }
    double r_ab   = AerodynamicResistance(t_wind[i], CONST_RL, Ele[i].windH, 10.);     // eq 4.2.25  [min m-1]
    double r_av   = AerodynamicResistance(t_wind[i], t_rl[i], Ele[i].windH, 10.);     // eq 4.2.25  [min m-1]
    
    double Gamma = PsychrometricConstant(Ele[i].FixPressure, lambda); // eq 4.2.28  [kPa C-1]
    qPotEvap[i] = gc.cETP * PET_Penman_Monteith(Ele[i].FixPressure, (t_rn[i] - 0) * 1.0e-6, rho,
                                            ed, Delta, r_ab, 0., Gamma, lambda);                // eq 4.2.27
    qPotTran[i] = gc.cETP * PET_Penman_Monteith(Ele[i].FixPressure, (t_rn[i] - 0) * 1.0e-6, rho,
                                            ed, Delta, r_av, r_s, Gamma, lambda);                // eq 4.2.27
    qEleETP[i] = qPotTran[i] * Ele[i].VegFrac + qPotEvap[i] * (1. - Ele[i].VegFrac);
}
void Model_Data::ET(double t, double tnext){
    double  T=NA_VALUE,  LAI=NA_VALUE, MF =NA_VALUE, prcp = NA_VALUE;
    double  snFrac, snAcc, snMelt, snStg;
    double  icAcc, icEvap, icStg, icMax, vgFrac;
    double  DT_min = tnext - t;
    int i;
#ifdef _OPENMP_ON
#pragma omp for parallel default(shared) private(i) num_threads(CS.num_threads)
#endif
    for(i = 0; i < NumEle; i++) {
        T = t_temp[i];
        prcp = t_prcp[i];
        /* Snow Accumulation */
        MF = t_mf[i];
        snStg = yEleSnow[i];
        /* Snow Accumulation/Melt Calculation*/
        if( T < Ts){
            snFrac = 1.;
        }else if(T > Tr){
            snFrac = 0.;
        }else{
            snFrac = (Tr - T) / (Tr - Ts);
        }
        snAcc = snFrac * prcp;
        snMelt = (T > To ? (T - To) * MF : 0.);    /* eq. 7.3.14 in Maidment */
        snMelt = min(max(0., snStg / DT_min), max(0., snMelt));
        snStg += (snAcc - snMelt) * DT_min;
        
        /* Interception */
        LAI = t_lai[i];
        icStg = yEleIS[i];
        vgFrac = Ele[i].VegFrac;
        if(LAI > ZERO){
            icMax = gc.cISmax * IC_MAX * LAI;
            icAcc = min(prcp - snAcc, max(0., (icMax - icStg) / DT_min) );
            icEvap = min(max(0., icStg / DT_min), qEleEvapo[i]);
        }else{
            icAcc = 0.;
            icEvap = 0.;
        }
        icStg += (icAcc - icEvap) * DT_min;
        
        if(icAcc>0){
            i=i;
        }
        /* Update the storage value and net precipitaion */
        yEleIS[i] = icStg * vgFrac;
        yEleSnow[i] = snStg;
        qEleE_IC[i] = icEvap * vgFrac;
        qEleNetPrep[i] = prcp + snMelt - snAcc - icAcc * vgFrac ;
        CheckNonNegative(qEleNetPrep[i], i, "Net Precipitation");
    }
}
void Model_Data::f_etFlux(int i, double t){
    double  LAI = t_lai[i];
    double Es = 0., Eu = 0., Tu = 0., Eg = 0., Tg = 0.;
    double va = Ele[i].VegFrac, vb = 1. - Ele[i].VegFrac;
    double pj = 1. - Ele[i].ImpAF;
    iBeta[i] = SoilMoistureStress(Soil[(Ele[i].iSoil - 1)].ThetaS, Soil[(Ele[i].iSoil - 1)].ThetaR, Ele[i].u_satn);
    /* Evaporation from SURFACE ponding water */
    Es = min(max(0., uYsf[i]), qPotEvap[i]);
    /* Evaporation from Unsaturated Zone. */
    Eu = min(max(0., uYus[i]), iBeta[i] * qPotEvap[i]) * pj * vb;
    /* Evporation from GroundWater, ONLY when gw above wetland level*/
    if(uYgw[i] > Ele[i].WetlandLevel){
        Eg = min(max(0., uYgw[i]), 1. * qPotEvap[i]) * pj * vb;
    }
    /* Vegetation Transpiration */
    if(LAI > ZERO){
        Tu = min(max(0., uYus[i]), iBeta[i] * qPotTran[i]) * pj * va;
        if(uYgw[i] > Ele[i].RootReachLevel){
            Tg = min(max(0., uYgw[i]), 1. * qPotTran[i]) * pj * va;
        }
    }else{
        /* VOID */
    }
    qEs[i] = Es;
    qEu[i] = Eu;
    qEg[i] = Eg;
    qTu[i] = Tu;
    qTg[i] = Tg;
    qEleE_sf[i] = Es;
    qEleTrans[i] = Tg + Tu;
    qEleEvapo[i] = Eu + Eg;
    qEleETA[i] = qEleE_IC[i] + qEleE_sf[i] + qEleEvapo[i] + qEleTrans[i];
#ifdef DEBUG
    CheckNonNegative(Es, i, "Es"); // Debug Only
    CheckNonNegative(Eu, i, "Eu");
    CheckNonNegative(Eg, i, "Eg");
    CheckNonNegative(Tu, i, "Tu");
    CheckNonNegative(Tg, i, "Tg");
    CheckNANi(ETp, i, "Potential ET (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleEvapo[i], i, "Transpiration (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleTrans[i], i, "Soil Evaporation (Model_Data::EvapoTranspiration)");
#endif
}

