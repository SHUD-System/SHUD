#include "is_sm_et.hpp"

double PET_Hargreaves(double Rad, double Tmax, double Tmin,
                      double Tavg, double lambda, double x6,
                      double x7, double x8, double x9){
//    x6-x9 are dummy values
//    SWAT 2:2.2.24
//    λ is the latent heat of vaporization (MJ kg-1),
//    Eo is the potential evapotranspiration (mm d-1),
//    H0 is the extraterrestrial radiation (MJ m-2 d-1),
//    Tmx is the maximum air temperature for a given day (􏰊C),
//    Tmn is the minimum air temperature for a given day (􏰊C), and
//    Tavg is the mean air temperature for a given day (􏰊C).
    return 0.023 * Rad * sqrt(Tmax - Tmin) * (Tavg - 17.8) / lambda;  /* SWAT 2:2.2.24 */
}
double PET_Priestley_Taylor(double alpha, double Delta, double Gamma,
                            double Hnet, double G, double lambda,
                            double x7, double x8, double x9){
//    x7-x9 are dummy values
//    SWAT 2:2.2.23
//    λ is the latent heat of vaporization (MJ kg-1),
//    Eo is the potential evapotranspiration (mm d-1), 􏰔
//    pet is a coefficient,
//    Δ is the slope of the saturation vapor pressure-temperature curve,
//    de/dT (kPa  ̊C-1),
//    γ is the psychrometric constant (kPa  ̊C-1),
//    Hnet is the net radiation (MJ m-2 d-1), and
//    G is the heat flux density to the ground (MJ m-2 d-1).
    return alpha * Delta / (Delta + Gamma) * (Hnet - G) / lambda;   /* SWAT 2:2.2.23 */
}
double PET_Penman_Monteith(double Rad, double rho,
                       double ed, double Delta, double r_a, double r_s,
                       double Gamma, double lambda){
    /* http://www.fao.org/docrep/X0490E/x0490e06.htm#penman%20monteith%20equation
     Rn net radiation at the crop surface [MJ m-2 s-1],
     G soil heat flux density [MJ m-2 day-1],
     es saturation vapour pressure [kPa],
     ea actual vapour pressure [kPa],
     es - ea saturation vapour pressure deficit [kPa],
     ∆ slope vapour pressure curve [kPa °C-1],
     γ psychrometric constant [kPa °C-1].
     r_s;    (bulk) surface resistances. [s m-1]
     r_a;    aerodynamic resistances. [s m-1]
    Gamma;  psychrometric constant
    Delta;  the slope of the saturation vapour pressure temperature relationship
     */
    double E_rad, E_air, r_sa;
    E_rad = Delta * Rad ;
    E_air = rho * Cp * ed / r_a;
    r_sa = r_s / r_a;
    double    ETp = ( E_rad + E_air) / (Delta + Gamma * (1 + r_sa)); // eq 4.2.27 [ ]
    ETp = ETp / lambda; // [ kg/s = 0.001 m3/m2/s = mm/s ]
    ETp = ETp * 0.001;  /* [kg s-1] to [m s-1] */
    return ETp;
}
double PET_PM_openwater(double Delta, double Gamma, double lambda,
                        double Rad, double ed, double U2){
    // eq 4.2.29 and eq 4.2.30 [ ]
    double    ETp = (Delta * Rad * SecADay + Gamma * 6.43 * (1.0 + 0.536 * U2) * ed) / (Delta + Gamma); // eq 4.2.30 [ mm/day ]
    ETp = ETp / lambda; // [ mm / day ]
    ETp = ETp * 0.001 / SecADay;  /* [mm d-1] -> [m s-1] */
    return ETp;
}

//double PET_Penman_Monteith(double RH,
//                       double T, double Vel, double P,
//                       double Rn, double rl, double windH,
//                       double Gamma, double lai, double R_ref){
//    /* http://www.fao.org/docrep/X0490E/x0490e06.htm#penman%20monteith%20equation
//     Rn net radiation at the crop surface [MJ m-2 min-1],
//     G soil heat flux density [MJ m-2 min-1],
//     T mean daily air temperature at 2 m height [°C],
//     u2 wind speed at 2 m height [m s-1],
//     es saturation vapour pressure [kPa],
//     ea actual vapour pressure [kPa],
//     es - ea saturation vapour pressure deficit [kPa],
//     ∆ slope vapour pressure curve [kPa °C-1],
//     γ psychrometric constant [kPa °C-1].
//    */
//    double Delta, A;
//    double  G = 0.;  // Soil heat flux. Assume G = 0;
//    //    double    Gamma;  // psychrometric constant
//    //    double    Delta;  // the slope of the saturation vapour pressure temperature relationship
//    double    r_s;    // (bulk) surface resistances. [min m-1]
//    double    r_a;    // aerodynamic resistances. [min m-1]
//    double    ETp;    //
//    double    ea, es, ed; //represents the vapour pressure deficit of the air
//    double  rho = 1.225;    // Air density at sea level at 15C.
//
//    es = VaporPressure_Sat(T);              // eq 4.2.2 [kpa]
//    ea = es * RH;   // [kpa]
//    ed = es - ea ;  // [kpa]
//    Delta = SlopeSatVaporPressure(T, es);   // eq 4.2.3 [kPa C-1]
//    rho = AirDensity(P, T);                 // eq 4.2.4 [kg m-3]
//    r_a = AerodynamicResistance(Vel, rl, windH, 2.); // eq 4.2.25  [min m-1]
//    if(lai > 0.){
//        r_s = BulkSurfaceResistance(R_ref, lai);
//    }else{
//        r_s = 0.;
//    }
//    A = Rn - G; // eq 4.2.16 [MJ m-2 day-1]
//    ETp = (Delta * (Rn - G) + rho * Cp * ed / r_a)
//    / (Delta + Gamma * (1 + r_s / r_a)); // eq 4.2.27 [ ]
//    ETp = ETp / LAMBDA * 0.001 ; // mm/min =>> m/min;   kg/m2/min = mm/min
//    return ETp;
//}

double PlantCoeff(double Rs_ref, double Rmin, double LAI, double T, double r_a,
                  double Rn, double es, double ea,
                  double beta_s, double Gamma, double Delta){
    double Rmax = 5000.0 / 60. ; //s/m to min/m;  Rmin in min/m --convert in readin.
    double f_r, alpha_r, gamma_s, r_s, P_c, eta_s;
    f_r = 1.1 * 1.5 * Rn / (Rs_ref * LAI);
    f_r = f_r < 0 ? 0 : f_r;
    alpha_r = (1 + f_r) / (f_r + Rmin / Rmax);
    alpha_r = min(alpha_r, 10000.);
    eta_s = 1. - 0.0016 * (24.85 - T) * (24.85 - T);
    eta_s = max(0.0001, eta_s);
    gamma_s = 1. / (1. + 0.00025 * (es - ea));
    gamma_s = max(0.01, gamma_s);
    if(beta_s > 0){
        r_s = Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s);
        r_s = min(Rmax, r_s);
    }else{
        r_s = Rmax;
    }
    P_c = (1. + Delta / Gamma) / (1. + r_s / r_a + Delta / Gamma);
    P_c = min(max(0., P_c), 1.);
    return P_c;
}
double SoilMoistureStress(double ThetaS, double ThetaR, double SatRatio){
    double fc, beta_s;
    fc = ThetaS * FieldCapacityRatio; /* Assumed field capacity 75% */
    beta_s = (SatRatio * (ThetaS - ThetaR) - ThetaR) / (fc - ThetaR);
    beta_s = min(max(0., beta_s), 1.);
    beta_s = 0.5 * (1- cos(PI * beta_s));
    /* R test code:
     b=seq(0, 1, 0.01)
     plot(b, 0.5 * (1-cos(pi * b)))
     */
    return beta_s;
}

