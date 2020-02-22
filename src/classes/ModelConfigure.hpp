//  ModelConfigure.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef ModelConfigure_hpp
#define ModelConfigure_hpp
#include "Macros.hpp"
#include "River.hpp"
#include "Node.hpp"
#include <string.h>

class calib_soil{
public:
    double infKsatV = 1.0;
    double macKsatV = 1.0;
    double infD = 1.0;
    double Alpha = 1.0;
    double Beta = 1.0;
    double hAreaF = 1.0;
};
class calib_geol{
public:
    double KsatH = 1.0;
    double KsatV = 1.0;
    double macKsatH = 1.0;
    double macD = 1.0;
    double ThetaS = 1.0;
    double ThetaR = 1.0;
    double vAreaF = 1.0;
};
class calib_landcover{
public:
    double VegFrac = 1.0;
    double Albedo = 1.0;
    double Rough = 1.0;
    double SoilDgd = 1.0;
    double RzD = 1.0;
    double ImpAF = 1.0;
//    double cLAI = 1.0; /* TS LAI */
    double cISmax = 1.0;
};


class Soil_Layer {
public:
    int     index;    /* index */
    double  infKsatV;    /* vertical saturated soil conductivity [m/day]*/
    double  ThetaS;    /* soil porosity [m3/m3] */
    double  ThetaR;    /* soil moisture residual [m3/m3] */
    double  ThetaFC;    /* Field Capacity [m3/m3] */
    double  Alpha;    /* soil curve parameter 1 [1/m]*/
    double  Beta;    /* soil curve parameter 2 [-] */
    double  hAreaF;    /* macroporous area fraction on horizontal section [m2/m2]*/
    double  macKsatV;    /* macroporous saturated vertical conductivity [m/day]*/
    double  infD;    /* depth from ground surface accross which head is calculated during infiltration [m]*/
    void    applyCalib(calib_soil *);
    void    checkValue();
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};

class Geol_Layer {
public:
    int     index;    /* index */
    double  KsatH;    /* horizontal saturated geology conductivity [m/day] */
    double  KsatV;    /* vertical saturated geology conductivity [m/day] */
    double  Sy;     /* Specific Yield  = ThetaS - ThetaR */
    double  geo_ThetaS;    /* geology porosity [m3/m3] */
    double  geo_ThetaR;    /* residual porosity [m3/m3] */
    double  geo_vAreaF;    /* macroporous area fraction on vertical section [m2/m2] */
    double  macKsatH;   /* macroporous saturated horizontal conductivity [m/day] */
    double  macD;       /* Depth of Macropore layer */
    void    applyCalib(calib_geol *);
    void    checkValue();
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};

class Landcover {
public:
    int index;    /* index */
//    double LAImax = NA_VALUE;    /* max LAI [m2/m2]*/  // debug: MUST remove this one, dummy.
    double VegFrac = NA_VALUE;/* Vegetation fraction [m2/m2] */
    double Albedo = NA_VALUE;    /* Albedo [1]*/
    double Rs_ref = NA_VALUE;
    double Rmin = NA_VALUE;    /* Minimum stomatal resistance */
    double Rough = NA_VALUE;    /* surface roughness factor */
    double RzD = NA_VALUE;    /* rootZone Depth [m] */
    double SoilDgrd = NA_VALUE; /* Soil degradation ratio [-] */
    double ImpAF = NA_VALUE; /* Impervious Area Fraction [m2/m2]*/
    void    applyCalib(calib_landcover *);
    void    checkValue();
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};


class globalCal{
    //public river_calib,
    //public calib_geol,
    //public calib_soil,
    //public calib_landcover {
public:
    calib_river criv;
    calib_geol cgeol;
    calib_soil csoil;
    calib_landcover clandc;
    
    double cAqD = 0; // +
    double cTemp = 0; // +
    double c_ic_gw = 0;
    double c_ic_riv = 0;
    
    double cETP = 1;
    double cPrep = 1;
    double cE_ic = 1;
    double cE_trans = 1;
    double cE_Evapo = 1;
    double cISmax = 1; /* */
    double cLAItsd = 1;  /* LAI TSD */
    double cMF = 1;  /* MF TSD */
    
    void copy(globalCal *p);
    void copy(const char **varname, int nvar,  double *x, int nx);
    void write(const char *fn);
    void read(const char *fn);
    void push(const char *var, double x);
    double getValue(const char *var);
    globalCal();
} ;


#endif /* ModelConfigure.hpp */

