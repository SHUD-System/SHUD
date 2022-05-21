//  River.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef River_hpp
#define River_hpp

#include <stdio.h>
#include "Macros.hpp"
#include "Node.hpp"
#include "functions.hpp"
#include "Equations.hpp"

class calib_river{
public:
    double rivRough = 1.0;
    double rivBankSlope = 1.0;
    double rivCwr = 1.0;
    double rivKsatH = 1.0;
    double rivDepth = 1.0;
    double rivWidth = 1.0;
    double rivSINU = 1.0;
    double rivBedThick = 1.0;
};
class river_para {
public:
    int index;
    double depth;       /* depth [m]*/
    double bankslope;   /* Bank slope [m/m] */
    double BottomWidth; /* Bottom Width [m] */
    double Sinuosity;   /* Sinuosity of river [m/m] */
    double rivRough;       /* Manning's roughness coeff [day/ m^(1/3)] */
    double Cwr;         /* Weir Discharge Coefficient [-]*/
    double KsatH;       /* Conductivity of river banks [m/day] */
    double BedThick;    /* Thickness of river bed sediment --- for ground water exchange [m] */
    
    
    void InitValue(double *x);
    void applyCalib(calib_river *);
    river_para();
    river_para( const river_para &obj);  // copy constructor
    void printInfo(FILE *fp);
    void printHeader(FILE *fp);
};

class _River : public river_para{
public:
    int index;
    
    /* Geometry */
//    int FromNode = NA_VALUE;    /* Upstream Node no. [-] */
//    int ToNode = NA_VALUE;    /* Dnstream Node no. [-] */
    int down = NA_VALUE;    /* down stream reach [-] */
    int type = NA_VALUE;    /* shape type [-] */
    int BC = NA_VALUE;    /* BC type [-]
                           BC > 0 = Neumann BC; Fix river stage
                           BC < 0 = Dirichlet BC; Fix Flux in/out*/
    int toLake = NA_VALUE, frLake = NA_VALUE;
    double yBC = 0.;
    double qBC = 0.;
    
    int reservoir = NA_VALUE; /* [-] */
//    double x = NA_VALUE;    /* Centroid x of river reach [m] */
//    double y = NA_VALUE; /* Centroid y of river reach [m] */
    double zbed = NA_VALUE;    /* bed elevation [m] */
    double zbank = 0.0;    /* bank elevation [m] */
    double Length = NA_VALUE;    /* Riv reach Length [m] */
    double BedSlope = NA_VALUE; /* Slope of river bed [m/m] */
    /* relation to its downstream */
    double Dist2DownStream = NA_VALUE; /* */
    double avgRough = NA_VALUE; /* Average of Mannin's n between this and downstream [day/ m^(1/3)] */
    
    /* Gemotry update with Ystage */
    double u_Ystage = NA_VALUE; /* River stage [m] */
    double u_CSperem = NA_VALUE; /* Cross-section peremeter [m] */
    double u_CSarea = NA_VALUE; /* Cross-section area [m2] */
    double u_eqWidth = NA_VALUE; /* Equivalent width of river [m] */
    double u_TopArea = NA_VALUE; /* Equivalent top area of river [m2] */
    double u_topWidth = NA_VALUE; /* top width of river [m] */
    /*methods */
    _River();
//    void initialRiver(river_para *, _Node *);
    void initialRiver(river_para *);
//    void applyNode(_Node *Node);
    void applyParameter(river_para *paras);
    void updateFrDownstream(_River *nabr);
    void updateRiver(double newY);
    void printInfo(FILE *fp);
    void printHeader(FILE *fp);
private:
    void checkGeomtry(void);
};

class RiverSegement{
public:
    int index;
    int iRiv = NA_VALUE;
    int iEle = NA_VALUE;
    double length = NA_VALUE;
    double eqDistance = NA_VALUE;
    double Cwr = NA_VALUE; /* Weir Discharge Coefficient */ // should be remove. debug.
    double KsatH = NA_VALUE; /* Conductivity of river banks */
};


double fun_CrossArea(double y, double w0, double s);
double fun_CrossPerem(double y, double w0, double s);
double fun_EqWidth(double y, double w0, double s);
double fun_TopWidth(double y, double w0, double s);




inline double fun_CrossArea(double y, double w0, double s){
    return y * (w0 + y * s);
}
inline double fun_CrossPerem(double y, double w0, double s){
    return 2.0 * sqpow2(y, y * s) + w0;
}
inline double fun_EqWidth(double y, double w0, double s){
    double w1 = fun_TopWidth(y, w0, s);
    return 0.5 * (w1 + w0);
}
inline double fun_TopWidth(double y, double w0, double s){
    return y * s * 2.0 + w0;
}

#endif /* River_hpp */
