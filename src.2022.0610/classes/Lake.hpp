/*******************************************************************************
 * File        : Lake.hpp                                                      *
 * Version     : Jul, 2018                                      *
 * Developer   :     Lele Shu (lele.shu@gmail.com)                  *
 *-----------------------------------------------------------------------------*
 *******************************************************************************/
#ifndef Lake_hpp
#define Lake_hpp

#include <stdio.h>
#include <iostream>
#include "Element.hpp"
#include "TabularData.hpp"
#include "functions.hpp"

class Orifice_Outlet{
public:
    // for details: http://content.alterra.wur.nl/Internet/webdocs/ilri-publicaOes/ publicaOes/Pub20/pub20-h8.0.pdf
    double  C0 = 0.61 ; //orifice coefficient (default 0.61);
    double  Area; // the orifice area;
    double  Hin;
    double  Hout;
    double  Q;
    Orifice_Outlet();
    void init(double c0, double a);
    double flow(double h0, double h1);
};
class Weir_Outlet{
public:
    // http://epg.modot.org/files/b/bc/749_Broad-Crested_Weir_Coefficients.pdf)
    double Cwr = 2.64;
    double Width;
    double h;
    double Q;
    Weir_Outlet();
    void init(double cw, double wd);
    double flow(double hx);
};
class Open_Outlet{
public:
    double Cwr = 0.6;
    double Width;
    double Hlake;
    double HRiver;
    double Q;
    Open_Outlet();
    void init(double cw, double wd);
    double flow(double hl, double hr);
};
class lake_calib{
public:
    double cksat = 1.;
};
class RiverIn{
public:
    
};
class RiverOut{
public:
};

class LakeBathymetry{
public:
    int     nvalue = NA_VALUE;
    int     *index; /* */
    double  *yi; /* Lake Stage */
    double  *ai; /* Top area */
    LakeBathymetry();
    ~LakeBathymetry();
    void InitValue(int n);
    void read(FILE *fp);
    double toparea(double y);
};

class _Lake{
public:
    double x = NA_VALUE;
    double y = NA_VALUE;
    double zmin = NA_VALUE;
    double area = NA_VALUE;
    
    int BC = 0; 
    int SS = 0;
    int NumEleLake = NA_VALUE;
    int NumEleBank = NA_VALUE;
    int NumRivIn = NA_VALUE;
    int NumRivOut = NA_VALUE;
    int *iEleLake;
    int *iEleBank;
    int *iRivIn;
    int *iRivOut;
    int *RivIn;
    int *RivOut;
    LakeBathymetry bathymetry;
//    RiverIn *RivIn;
//    RiverOut *RivOut;
    
    double *QEleSurf;
    double *QEleGW;
    double *QRivIn;
    double *QRivOut;
    
    double yStage;
    
    /* */
    double u_depth;
    double u_toparea;
    double u_volume;
    
    /* Methods */
    _Lake();
    ~_Lake();
    void Initialize();
    void update();
    void readLake(FILE *fp);
    void readBathymetry(FILE *fp);
};


#endif /* Lake_hpp */
