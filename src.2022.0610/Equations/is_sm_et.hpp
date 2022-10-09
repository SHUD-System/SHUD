//  is_sm_et.h
//
//  Created by Lele Shu on 6/25/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef is_sm_et_h
#define is_sm_et_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "ModelConfigure.hpp"
#include "Macros.hpp"
#include "functions.hpp"

//void is_sm_et(double t, double stepsize, void *DS, N_Vector VY);


double PET_Hargreaves(double Rad, double Tmax, double Tmin,
                      double Tavg, double lambda, double x6,
                      double x7, double x8, double x9);
double PET_Priestley_Taylor(double alpha, double Delta, double Gamma,
                            double Hnet, double G, double lambda,
                            double x7, double x8, double x9);
double PET_Penman_Monteith(double Rn, double rho,
                       double ed, double Delta, double r_a, double r_s,
                       double Gamma, double Lambda);
double PET_PM_openwater(double Delta, double Gamma, double lambda,
                        double Rad, double ed, double U2);
//double PET_Penman_Monteith(double RH,
//                       double Temp, double Vel, double Press,
//                       double Rn, double rl, double windH,
//                       double Gamma, double lai, double R_ref);

double PlantCoeff(double Rs_ref, double Rmin, double LAI, double Temp, double r_a,
                  double Rn, double es, double ea,
                  double beta_s, double Gamma, double Delta);
double SoilMoistureStress(double ThetaS, double ThetaR, double SatRatio);
double Eact_Interception(double etp, double yis, double ymax);
double Eact_Vegetation(double etp, double yis, double ymax);
double ActualEvaporation(double etp, double ThetaS, double ThetaR, double SatRatio);

double VaporPressure_Sat(double T_in_C);
double VaporPressure_Act(double esat, double rh);
double AerodynamicResistance(double Uz, double hc, double Z_u, double Z_e);
double SlopeSatVaporPressure(double Temp, double Es);
double AirDensity(double Temp, double Es);

double CanopyResistance(double rmin, double lai);
double PsychrometricConstant(double Pressure);
double PressureElevation(double z);
double BulkSurfaceResistance(double R_ref, double lai);



inline double Eact_Interception(double etp, double yis, double ymax){
    return etp * pow23(yis / ymax); // Deardorff (1978)
    /* Mukeksh(2009) eq 4.7c*/
}
inline double Eact_Vegetation(double etp, double yis, double ymax){
    double aet=0;
    return aet;
}
inline double ActualEvaporation(double etp, double ThetaS, double ThetaR, double SatRatio){
    return SoilMoistureStress(ThetaS, ThetaR, SatRatio) * etp;
}
inline double relativeHumidity(double sph, double Temp, double pres){
    return 0.;
}
inline double LatentHeat(double Temp){
    /* Eq 4.2.1 in David R Maidment, Handbook of Hydrology */
    /* Temp in [C], lambda in [MJ/kg]*/
    /* Latent heat of vaporization */
    return 2.501 - 0.002361 * Temp;
}
inline double BulkSurfaceResistance(double R_ref, double lai){
    // R_ref     bulk stomatal resistance of the well-illuminated leaf [s m-1],
    return R_ref * 2. / lai; /* Allen(1998) */
}
inline double BulkSurfaceResistance(double lai){
    /* Eq 4.2.22 in David R Maidment, Handbook of Hydrology */
    return 200. / lai;  /* [s m-1] */
}
inline double lai2hc(double lai){
    /* Eq 4.2.23 in David R Maidment, Handbook of Hydrology */
    /* From LAI to Height of crop */
    return exp((lai - 5.5)/1.5 );  /* [m2 m-2] */
}
inline double PressureElevation(double z){
    return 101.325 * pow((293. - 0.0065 * z) / 293, 5.26); /*Pressure based on Elevation*/
    //    P atmospheric pressure [kPa],
    //    z elevation above sea level [m],
    /* Allen(1998) Eq(7) */
}
inline double PsychrometricConstant(double Pressure, double lambda){
    /* Eq 4.2.28 in David R Maidment, Handbook of Hydrology */
    /* Presser in [kPa], lambda in [MJ/kg]. */
    return 0.0016286 * Pressure / lambda;
    //    γ psychrometric constant [kPa °C-1],
}
inline double VaporPressure_Sat(double T_in_C){
    /* Eq 4.2.2 in David R Maidment, Handbook of Hydrology*/
    return 0.6108 * exp( 17.27 * T_in_C / (T_in_C + 237.3));
}
inline double VaporPressure_Act(double esat, double rh){
    return esat * rh;
}
//double VaporPressure_Act(double RH, double TC, double P){
//    double VP = 611.2 * exp(17.67 * TC / (TC + 243.5)) * RH;
//    double qv_sat = 0.622 * (VP / RH) / P;
//    return qv_sat;
//}
inline double AerodynamicResistance(double Uz, double hc, double Z_u, double Z_e){
    /*
     Eq 4.2.25 in David R Maidment, Handbook of Hydrology
     Or:
     Allen, R. G., S, P. L., Raes, D., & Martin, S. (1998).
     Crop evapotranspiration : Guidelines for computing crop water requirements
     by Richard G. Allen ... [et al.].
     FAO irrigation and drainage paper: 56.
     Eq(4) in reference
        ra aerodynamic resistance [s m-1],
        zm height of wind measurements [m],
        zh height of humidity measurements [m],
        d zero plane displacement height [m],
        zom roughness length governing momentum transfer [m],
        zoh roughness length governing transfer of heat and vapour [m],
        k von Karman's constant, 0.41 [-],
        uz wind speed at height z [m s-1].
     */
    double  r_a, d, Z_om, Z_ov;
    d = 0.67 * hc;;
    Z_om = 0.123 * hc;
    Z_ov = 0.0123 * hc;
    r_a = log( fabs(Z_u - d) / Z_om ) * log( fabs(Z_e - d) / (Z_ov)) / (VON_KARMAN * VON_KARMAN* Uz);
    return r_a;
}
inline double AerodynamicResistance(double U2, double Zu, double Z0){
    /* Eq 4.2.29 in David R Maidment, Handbook of Hydrology */
    /* for open watet. Only the wind-measure-height availalble */
    double  temp;
    temp = log(Zu / Z0);
    return 4.72 * temp * temp / (1 + 0.536 * U2);
}

inline double WindProfile(double Zx, double Um, double Zm, double d, double Z0){
    /* Equation 2.3 in [@abtew2012evaporation] */
    return Um * log( (Zx - d) / Z0 ) / log ((Zm - d) /  Z0);
}
inline double AirDensity (double P, double Temp){
    /* Eq 4.2.4 in David R Maidment, Handbook of Hydrology*/
    /*  P  in [kP]
        Temp in [C]
        rho  -- Density of Air. kg/m3
     */
    return 3.486 * P / (275. + Temp);
}
inline double SlopeSatVaporPressure(double Temp, double ES){
    double tt = (Temp + 237.3);
    double delta = 4098. * ES /  ( tt * tt);
    return delta;
}

inline double CanopyResistance(double rmin, double lai){
    double rc;
    rc = rmin * 2. / lai;
    return rc;
}


#endif /* is_sm_et_h */
