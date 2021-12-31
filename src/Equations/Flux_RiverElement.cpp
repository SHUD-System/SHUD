//
//  Flux_Subsurface.cpp
//  SHUD
//
//  Created by Lele Shu on 1/18/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Flux_RiverElement.hpp"

double flux_R2E_GW(double yr, double zr,
                   double ye, double ze,
                   double Kele, double Kriv,
                   double L, double D_riv){
    double dh, A,  g, K, he, hr;
    double Q = 0.0;
//    K = (Kele + Kriv) * 0.5;  //debug ???
    if(Kele < ZERO || Kriv < ZERO){
        return 0.;
    }else{
//        K = meanHarmonic(Kele, Kriv, 1., D_riv);
        K= meanArithmetic(Kele, Kriv, 1., 1.);
    }
    he = ye + ze; //head of Ele GW
    hr = yr + zr; //head of river
    dh = hr - he;
    if( dh > ZERO){
        /* River to Element*/
        if( he > zr ){ // Element gw higher than river bed
            A = (yr + (he - zr) ) * .5 * L;
        }else{
            A = yr * L;
        }
        if(yr < EPSILON){
            Q = 0.;
        }else{
            g = dh / D_riv;
            Q = A * K * g;
        }
    }else if (dh < -ZERO){
        /* Element to River */
        if( ye > ZERO){
            A = (yr + (he - zr) ) * .5 * L;  /* Mean of River stage and he-zr  */
            g = dh / D_riv;
            Q = A * K * g;
        }else{
            Q = 0.;
        }
    }else{
        Q = 0.;
    }
    
//    CheckNANi(Q, 1, "River Flux Sub(Functopm:fun_Seg_sub)");
    return Q;
}
double flux_R2E_SF(double yr, double zr, double ye, double ze, double rough, double Len, double dist, double threshold){
    double dh=0., A=0., Q=0., P=0., he=0., hr=0.;
    he = ye + ze;
    hr = yr + zr;
    dh = hr - he;
    if( dh <= 0.){ // flow from Element to River
        if(ye <= threshold || dh == 0.){
            Q = 0.;
        }else{
//            A = Len * ye;
//            P = Len ;
//            Q = -1.0 * ManningEquation(A, rough,  A / P, ye / dist);
            P= cbrt(ye);
            Q = -1. / rough * Len * P* P* P* P* P *sqrt(ye/dist);
        }
    }else{
        A = Len * (hr-ze + ye) * 0.5;
        P = Len + (hr-ze + ye);
        Q = ManningEquation(A , rough,  A/P, ye / dist);
    }
    return Q;
}
