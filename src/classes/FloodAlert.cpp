//  FloodAlert.cpp
//
//  Created by Lele Shu on 9/7/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "FloodAlert.hpp"

void FloodAlert::InitAlert(int nr, int nrp){
    nriv = nr;
    pstage = new double*[nriv];
    pflux = new double*[nriv];
    itype = new int[nriv];
    
//    rivStage = new double[nriv];
//    rivBank = new double[nriv];
//    Q_high = new double[nriv];
//    Q_low = new double[nriv];
//    Y_high = new double[nriv];
//    Y_low = new double[nriv];
    
    ntype = nrp;
    para  = new river_para[nrp];
}
FloodAlert::FloodAlert(){}
FloodAlert::~FloodAlert(){
    if(pstage != NULL) delete[] pstage;
    if(pflux != NULL) delete[] pflux;
    if(itype != NULL) delete[] itype;
    if(para != NULL) delete[] para;
//    delete[]  rivStage;
//    delete[]  rivBank;
//    delete[]  Q_high;
//    delete[]  Q_low;
//    delete[]  Y_high;
//    delete[]  Y_low;
    
    fclose(fid);
}
void FloodAlert::pushRiverType(int index, int type){
    itype[index] = type - 1;
}
void FloodAlert::InitPara(river_para *p){
    for(int i = 0; i < ntype; i++){
        para[i] = p[i];
//        para[i].BottomWidth = p[i].BottomWidth;
//        para[i].depth = p[i].depth;
//        para[i].index = p[i].index;
//        para[i].Sinuosity = p[i].BottomWidth;
//        para[i].Rough = p[i].Rough;
//        para[i].Cwr = p[i].Cwr;
//        para[i].KsatH = p[i].KsatH;
    }
    
}
void FloodAlert::InitPointer(double *y, double *q){
    for (int i = 0; i < nriv; i++){
        pstage[i] = &(y[i]);
        pflux[i] = &(q[i]) ;
    }
}
void FloodAlert::InitFile(const char *fn){
    fid = fopen(fn, "w");
    CheckFile(fid, fn);
    fprintf(fid, "%s\t%s\t%s\t%s\t%s\t%s\n",
            "time", "ID","Type",
            "Stage_m", "Bank_m",
            "Discharge_m3/day");
}

int FloodAlert::FloodWarning(double t){
    return FloodWarning(t, fid);
}
inline  int FloodAlert::FloodWarning(double t, FILE *fp){
    /*If there is flood, print warning. return 1 - flood, 0 - no flood*/
    int flag = 0;
    for(int i = 0; i < nriv; i++){
        if( *(pstage[i]) > para[ itype[i] ].depth ){
            flag = 1;
            fprintf(fp, "%.1f\t%d\t%d\t%f\t%f\t%f\n",
                    t, i + 1, itype[i] + 1,
                    *(pstage[i]), para[ itype[i] ].depth,
                    *(pflux[i]));
        }
    }
    return flag;
}

