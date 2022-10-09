//  Equations.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef Equations_hpp
#define Equations_hpp

#include <stdio.h>
#include <math.h>
#include "Macros.hpp"
#include "functions.hpp"


//double returnVal(double rArea, double rPerem, double eqWid, double ap_Bool);
//double CS_AreaOrPerem(int rivOrder, double rivDepth, double rivCoeff, double a_pBool);
double avgY_sf(double z1, double y1, double z2, double y2, double threshold);
double avgY_gw(double z1, double y1, double z2, double y2, double threshold);
double effKV(double ksatFunc, double gradY, double macKV, double KV, double areaF);
double effKVnew(double ksatFunc, double macKV, double KV, double areaF, double y0);

double effKH(double tmpY, double aqDepth, double MacD, double MacKsatH, double areaF, double ksatH);
double satKfun(double elemSatn, double beta);
double sat2psi(double elemSatn, double alpha, double beta, double mpsi);
double fun_recharge(double effk_us, double kgw, double Deficit, double ygw, double hus,double yus);
double effKRech(double ksatFunc,  double macKV, double KV, double areaF);
double ManningEquation(double Area, double rough, double y, double S);
double OverlandManning(double avg_y, double grad_y, double avg_sf, double A, double n);
double FieldCapacity (double Alpha, double Beta, double deficit );
double GreenAmpt(double k, double ti, double ts, double phi, double hf, double prcp, double sy);
inline
double sat2psi(double elemSatn, double alpha, double n){
    return -(pow(pow(elemSatn, n / (1 - n)) - 1, 1 / n) / alpha);
}

inline
double pow23(double x){
    double t = cbrt(x);
    return t * t;
}
inline
double sqpow2(double x, double y){
    return sqrt(x * x + y * y);
}
inline
double meanHarmonic(double k1, double k2, double d1, double d2){
//    return (d1 + d2) / ( d1 / k1 + d2 / k2);
    return (k1 * k2) * (d1 + d2) / ( d1 * k2 + d2 * k1);
}
inline
double meanArithmetic(double k1, double k2, double d1, double d2){
    return (k1 * d1 + k2 * d2) / (d1 + d2);
}
inline
double ManningEquation(double Area, double rough, double R, double S){
    //    double Q = 0.;
    if (S > 0) {
        return sqrt(S) * Area * pow23(R) / rough;
        //return sqrt(S) * Area * pow(Area / Perem, 2. / 3.) / rough;
    } else {
        return -1.0 * sqrt(-S) * Area * pow23(R) / rough;
    }
    //    return Q;
}

inline
double TemperatureOnElevation(double t, double Zi, double Zt){
    if( ifequal(Zi, NA_VALUE) || ifequal(Zt, NA_VALUE) ){
        return t;
    }else{
        return t + (Zt - Zi) * dTdZ;
    }
}
#endif /* Equations_hpp */
