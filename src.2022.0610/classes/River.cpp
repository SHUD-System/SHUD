//  River.cpp
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "River.hpp"

river_para::river_para(){
    
}
river_para::river_para(const river_para &obj){
    index = obj.index;
    depth = obj.depth;
    bankslope = obj.bankslope;
    BottomWidth = obj.BottomWidth;
    Sinuosity = obj.Sinuosity;
    rivRough = obj.rivRough;
    Cwr = obj.Cwr;
    KsatH = obj.KsatH;
    BedThick = obj.BedThick;
}
void river_para::InitValue(double *x){
    index       = (int) x[0];
    depth       = (double) x[1];
    bankslope   = (double) x[2];
    BottomWidth = (double) x[3];
    Sinuosity   = (double) x[4];
    rivRough    = (double) x[5] / 60.; /* [s m^(-1/3)]  to [min m^(-1/3)]*/
    Cwr         = (double) x[6];
    KsatH       = (double) x[7] / 1440.; /* [m/d] to [m/min] */
    BedThick    = (double) x[8];
//    printf("debug%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", index, depth, bankslope, BottomWidth, Sinuosity, rivRough, Cwr, KsatH);
}
void river_para::applyCalib(calib_river *g){
    depth += g->rivDepth;
    BottomWidth += g->rivWidth;
    bankslope += g->rivBankSlope;
    Sinuosity *= g->rivSINU;
    rivRough *= g->rivRough;
    Cwr *= g->rivCwr;
    KsatH *= g->rivKsatH;
    BedThick *= g->rivBedThick;
}

_River::_River(){
    
}
void _River::updateRiver(double newY){
    u_Ystage = newY;
    u_topWidth = fun_TopWidth(u_Ystage, BottomWidth, bankslope);
    u_CSarea = fun_CrossArea(u_Ystage, BottomWidth, bankslope);
    u_CSperem = fun_CrossPerem(u_Ystage, BottomWidth, bankslope);
    u_eqWidth = fun_EqWidth(u_Ystage, BottomWidth, bankslope);
    u_TopArea = u_eqWidth * Length;
    
    u_topWidth = fixMaxValue(u_topWidth, 0.);
    u_CSarea = fixMaxValue(u_CSarea, 0.);
    u_CSperem = fixMaxValue(u_CSperem, 0.);
    u_eqWidth = fixMaxValue(u_eqWidth, 0.);
    u_TopArea = fixMaxValue(u_TopArea, 0.);    
}

void _River::applyParameter(river_para *paras){
    depth = paras[type - 1].depth;
    bankslope = paras[type - 1].bankslope;
    BottomWidth = paras[type - 1].BottomWidth;
    Sinuosity = paras[type - 1].Sinuosity;
    rivRough = paras[type - 1].rivRough;
    Cwr = paras[type - 1].Cwr;
    KsatH = paras[type - 1].KsatH;
    BedThick = paras[type - 1].BedThick;
}
void _River::updateFrDownstream(_River *DownRiv){
    int idown = down - 1;
    if( idown >= 0){
        avgRough = 0.5 *(rivRough + DownRiv[idown].rivRough);
        Dist2DownStream = 0.5 * (Length + DownRiv[idown].Length);
    }else{
        // outlets
        avgRough = rivRough;
        Dist2DownStream = Length;
    }
}

void _River::initialRiver(river_para *rpara){
//    BedSlope = DZ / Length;
    applyParameter(rpara);
//    applyNode(rivnode);
}
void _River::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "index");
//    fprintf(fp, "%s\t", "FromNode");
//    fprintf(fp, "%s\t", "ToNode");
    fprintf(fp, "%s\t", "down");
    fprintf(fp, "%s\t", "type");
    fprintf(fp, "%s\t", "reservoir");
    fprintf(fp, "%s\t", "BC");
//    fprintf(fp, "%s\t", "x");
//    fprintf(fp, "%s\t", "y");
    fprintf(fp, "%s\t", "zbed");
    fprintf(fp, "%s\t", "zbank");
    fprintf(fp, "%s\t", "Length");
    fprintf(fp, "%s\t", "BedSlope");
    fprintf(fp, "%s\t", "avgRough");
    river_para::printHeader(fp);
    fprintf(fp, "\n");
}
void _River::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
//    fprintf(fp, "%d\t", FromNode);
//    fprintf(fp, "%d\t", ToNode);
    fprintf(fp, "%d\t", down);
    fprintf(fp, "%d\t", type);
    fprintf(fp, "%d\t", reservoir);
    fprintf(fp, "%d\t", BC);
//    fprintf(fp, "%g\t", x);
//    fprintf(fp, "%g\t", y);
    fprintf(fp, "%g\t", zbed);
    fprintf(fp, "%g\t", zbank);
    fprintf(fp, "%g\t", Length);
    fprintf(fp, "%g\t", BedSlope);
    fprintf(fp, "%g\t", avgRough);
    river_para::printInfo(fp);
    fprintf(fp, "\n");
}

void river_para::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "index");
    fprintf(fp, "%s\t", "depth");
    fprintf(fp, "%s\t", "bankslope");
    fprintf(fp, "%s\t", "BottomWidth");
    fprintf(fp, "%s\t", "Sinuosity");
    fprintf(fp, "%s\t", "Rough");
    fprintf(fp, "%s\t", "Cwr");
    fprintf(fp, "%s\t", "KsatH");
//    fprintf(fp, "\n");
}
void river_para::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
    fprintf(fp, "%g\t", depth);
    fprintf(fp, "%g\t", bankslope);
    fprintf(fp, "%g\t", BottomWidth);
    fprintf(fp, "%g\t", Sinuosity);
    fprintf(fp, "%g\t", rivRough);
    fprintf(fp, "%g\t", Cwr);
    fprintf(fp, "%g\t", KsatH);
}
