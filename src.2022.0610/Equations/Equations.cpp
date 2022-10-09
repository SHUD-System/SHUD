//
//  Equations.cpp
//  Created by Lele Shu on 8/13/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Equations.hpp"
double avgY_sf(double z1, double y1, double z2, double y2, double threshold){
    double h1 = z1 + y1, h2 = z2 + y2;
//    double  dh = h1 - h2;
    if (h1 > h2) {
        if (y1 > threshold) {
//            if(y1 > dh){
//                return dh;
//            }else{
//                return y1;
//            }
//            if(y1 > y2 * 4.){ // VERY SLOW
//                return 0.5 * (y1 + y2);
//            }else{
//                return y1;
//            }
            //return ((yinabr > yi) ? 0 : 1.0 * yi);
            /* Note the if-else TRUE case can be possible only for
             * Kinematic case */
            return y1;
        } else {
            return 0.;
        }
    } else {
        if (y2 > threshold) {
//            if(y2 > -dh){
//                return -dh;
//            }else{
//                return y2;
//            }
//            if(y2 > y1 * 4. ){ // VERY SLOW
//                return 0.5 * (y1 + y2);
//            }else{
//                return y2;
//            }
            //return 0.5 * (yi + yinabr);
            //return ((yi > yinabr) ? 0 : 1.0 * yinabr);
            /* Note the if-else TRUE case can be possible only for
             * Kinematic case */
            return y2;
        } else {
            return 0.;
        }
    }
}
double avgY_gw(double z1, double y1, double z2, double y2, double threshold){
//    double h1 = z1 + y1, h2 = z2 + y2;
    y1 = max(y1, 0.);
    y2 = max(y2, 0.);
    return (y1 + y2) * .5;
//    if (h1 > h2) {
//        if (y1 > threshold) {
//            return y1;
//        } else {
//            return 0.;
//        }
//    } else {
//        if (y2 > threshold) {
//            return y2;
//        } else {
//            return 0.;
//        }
//    }
}
double effKV(double ksatFunc, double gradY, double macKV, double KV, double areaF)
{
    if (ksatFunc >= 0.98) {
        return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
    } else {
        if (fabs(gradY) * ksatFunc * KV <= 1 * KV * ksatFunc) {
            return KV * ksatFunc;
        } else {
            if (fabs(gradY) * ksatFunc * KV < (macKV * areaF + KV * (1 - areaF) * ksatFunc)) {
                return (macKV * areaF * ksatFunc + KV * (1 - areaF) * ksatFunc);
            } else {
                return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
            }
        }
    }
}
double effKRech(double ksatFunc,  double macKV, double KV, double areaF)
{
    double keff = 0.0;
    if(areaF <= 0. ){
        keff = KV * ksatFunc;
    }else{
        keff = (KV * (1 - areaF) + macKV * areaF) * ksatFunc;
    }
    return keff;
}
double effKVnew(double ksatFunc,  double macKV, double KV, double areaF, double y0)
{
    double keff = 0.0;
    
    if(areaF <= 0. ){
        keff = KV * ksatFunc;
    }else{
        keff = KV * (1 - areaF) + macKV * areaF;
    }
//    if(areaF <= 0. ){
//        keff = KV * ksatFunc;
//    }else if (y0 > EPSILON){
//        /*  Ponding water exists */
//        keff = KV * (1 - areaF) + macKV * areaF;
//    }else{
//        keff = (KV * (1 - areaF) + macKV * areaF) * ksatFunc;
//    }
    return keff;
}
double effKH(double Ygw, double aqDepth, double MacD, double Kmac, double AF, double Kmx){
    double effk = 0;
    if (MacD <= ZERO || Ygw < aqDepth - MacD) {
        effk = Kmx;
    } else {
        if (Ygw > aqDepth) {
            effk = (Kmac * MacD * AF +
                    Kmx * (aqDepth - MacD * AF)) / aqDepth;
        } else {
            effk = (Kmac * (Ygw - (aqDepth - MacD)) * AF +
                    Kmx * (aqDepth - MacD + (Ygw - (aqDepth - MacD)) * (1 - AF))) / Ygw;
        }
    }
    if( effk < 0. || effk > 1e9 ){
        fprintf(stderr, "Wrong effKH for ground water: %f [ m/s ]\n", effk / 60.);
        myexit(ERRDATAIN);
    }
    return effk;
}

double satKfun(double elemSatn, double n){
    double temp = -1. + pow(1. - pow(elemSatn, n / (n - 1.)), (n - 1.) / n);
    double ret = sqrt(elemSatn) * temp * temp;    
//    CheckNANi(ret, 0, "satKfun():ret");
    return ret;
}
/*************Critical update.*******************
double OverlandManning(double avg_y, double grad_y, double avg_sf, double A, double n){
    double Q;
    //flux[loci][locj] = crossA * pow23(avg_y) * grad_y / (sqrt(fabs(avg_sf)) * avg_rough);
    return Q;
} Critical update. The avg_sf of A->B differs from B->A. Which cause mass-balance issue.
*************/
double fun_recharge(double effk_us, double kgw, double Deficit, double ygw, double hus, double yus)
{
    double effk;
    double q = 0.;
    if(effk_us < ZERO || kgw < ZERO){
        return 0.;
    }
    effk = meanHarmonic(kgw, effk_us, ygw, yus);
    q = effk * (1. + hus / (Deficit * 0.5));
    return q;
}
double FieldCapacity (double Alpha, double Beta, double deficit ){
    //    sfc=((alpha *deficit)^(-beta) +1) ^(-1/beta)/alpha
    double    sfc;
    sfc = pow( pow(Alpha * deficit, -Beta) +1, -1/Beta) / Alpha;
    return sfc;
}
double GreenAmpt(double k, double ti, double ts, double phi, double hf, double h0, double Sy){
    double dTheta = ts - ti;
//    double dh = h0 + hf - phi;
    double q = 0.;
//    double D0;
    if(h0 <= 0 ){
        /* No rainfall, no ponding */
        return 0.;
    }
    hf = max(EPSILON, hf);
    if(dTheta <= 0.){
        q = 0.;
    }else{
        q = k * ( (h0 + hf - phi) * dTheta / hf  );
    }
    if( q >= 0. ){
        if( q > h0 ){
            q = h0;
        }else{
            q = q;
        }
    }else{
        q = q;
    }
    return q;
}

