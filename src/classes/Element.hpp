//  Element.hpp
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef Element_hpp
#define Element_hpp

#include <stdio.h>
#include "Macros.hpp"
#include "Node.hpp"
#include "ModelConfigure.hpp"
#include "Equations.hpp"
#include "is_sm_et.hpp"

class Element_Calib{
    calib_soil csoil;
    calib_geol cgeol;
    calib_landcover clandc;
};
class Triangle{
public:
    int node[3];/* anti-clock-wise */
    int nabr[3];/* neighbor i shares edge i (0: on boundary) */
    int lakenabr[3] = {0, 0, 0}; /* index to LAKE neighbor.  */
    int nabrToMe[3] = {-1, -1, -1}; /* The index j of my nabor */
    double edge[3];/* edge i is from node i to node i+1 */
    double area = NA_VALUE;    /* area of element */
    double slope[3]; /* Slope from centroid to edge */
    double Dist2Edge[3];
    double x = NA_VALUE;    /* x of centroid */
    double y = NA_VALUE;    /* y of centroid */
    double z_bottom = NA_VALUE;    /* Aquifer Bottom Elevation of the triangle centroid */
    double z_surf = NA_VALUE;    /* Surface Elevation of the triangle centroid */
    double zcentroid = NA_VALUE;
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};
/* Definition of Global Variable Types */
class AttriuteIndex{
public:
    int iSoil = NA_VALUE;   /* soil type */
    int iGeol = NA_VALUE;   /* geology type */
    int iLC = NA_VALUE;     /* Land Cover type */
    int IC  = NA_VALUE;     /* initial condition type */
    int iForc = NA_VALUE;   /* precipitation (forcing) type */
    int iMF = NA_VALUE;     /* meltFactor */
    int iSS = 0; /* Index of TS Source/Sink on LAND SURFACE only;
                 0=No SS;
                 BC>0 = SS for LANDSURFACE; such as irrigation. ;
                 BC<0 = SS for Ground water; such as pumping. */
    int iBC = 0; /* Index of TS Boundary Conditions;
                 0 = No BC;
                 BC>0 = Neumann / Fix head / ;
                 BC<0 = Dirichlet /Fix Fluxes */
    int iLake = NA_VALUE;
    int ilakebank = 0;
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};

class _Element : public Triangle,
public Soil_Layer,
public Geol_Layer,
public Landcover ,
public AttriuteIndex
{    /* Data model for a triangular element */
public:
    int index;    /* Element No. */
    int RivID = 0;
    int RivSegID = 0;
    double windH = NA_VALUE;    /* wind measurement height */
    /* for calculation of dh/ds */
//    double surfH[3];    /* Total head in neighboring cells */
//    double surfX[3];    /* Center X location of neighboring cells */
//    double surfY[3];   /* Center Y location of neighboring cells */
//    double dhBYdx = NA_VALUE;    /* Head gradient in x dirn. */
//    double dhBYdy = NA_VALUE;    /* Head gradient in y dirn. */
//    double Avg_Sf = NA_VALUE;    /* Head gradient in normal */
    double Dist2Nabor[3];
    double FixPressure = NA_VALUE;  /* Pressure [Pa]*/
//    double FixGamma = NA_VALUE;     /* Psychrometric Constant [kPa °C-1] */
    double AquiferDepth = NA_VALUE; // Zmax - Zmin
    double WetlandLevel = NA_VALUE; // Aquiferdepth - infD
    double RootReachLevel = NA_VALUE; //Aquiferdepth - RzD
    double MacporeLevel = NA_VALUE; //Aquiferdepth - macD
    double avgRough[3];
    double depression = 0.0002; //Depression value. No overland flow before filling the depression. Default = 0.2 mm.
    double yBC = 0; // Ground head in m
    double QBC = 0; // Flux in m3/day
    double QSS = 0; // Flux in m3/day
    
    int iupdGW[3] = {0,0,0}; /* whether the Groundwater Flux value on this J is update */
    int iupdSF[3] = {0,0,0};  /* whether the Surface Flux value on this J is update */
    /* Value must be updated each loop */
    double u_qi; /* infiltration from surface to unsat zone */
    double u_qex; /* exfiltration from groundwater to surface */
    double u_qr; /* recharge into ground water */
    /* GW flow*/
    double u_effKH; /* Horizontal effective flow, for gw*/
    double u_satn; /* Saturation ratio */
    double u_wf = 0.;
    double u_deficit; /* deficit. aquiferdepth - Ygw */
    double u_theta; /* Soil moisture content [m3/m3] */ 
private:
    /* Infiltration */
    double u_phius; /* pressure head of the unsat zone from zmin*/
    double u_Ginfi; /* Gradient to infiltration*/
    double u_satKr; /* ratio to effective Infiltration K */
    double u_effkInfi; /* Effective K for infiltration */
    double Kmax;
    
    /*==== Methods ===================*/
public:
    void InitElement();
    void copyGeol(Geol_Layer *);
    void copySoil(Soil_Layer *);
    void copyLandc(Landcover *);
    void applyGeometry(_Node *Node);
    void applyNabor(_Node *Node, _Element *Ele);
    void updateElement(double Ysurf, double Yunsat, double Ygw);
    void updateLakeElement();
//    void updateWF(double dh, double dt);
    void Flux_Infiltration(double Ysurf, double Yunsat, double Ygw, double netprcp);
    double Flux_Recharge(double Yunsat, double Ygw);
    void printHeader(FILE *fp);
    void printInfo(FILE *fp);
};
#endif /* Element_hpp */

