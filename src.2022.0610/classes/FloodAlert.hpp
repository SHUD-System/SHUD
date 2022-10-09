//  FloodAlert.hpp
//
//  Created by Lele Shu on 9/7/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef FloodAlert_hpp
#define FloodAlert_hpp

#include <stdio.h>
#include "River.hpp"

class FloodAlert{
public:
    FloodAlert();
    ~FloodAlert();
    river_para *para = NULL;
    
    void InitAlert(int nr, int nrp);
    void InitFile(const char *fn);
    void InitPara(river_para *p);
    void pushRiverType(int index, int type);
    void InitPointer(double *y, double *q);
    int FloodWarning(double t, FILE *fp); /*If there is flood, print warning. return 1 - flood, 0 - no flood*/
    int FloodWarning(double t);
    
private:
    int nriv = 0;;   /* Number of rivers */
    int ntype = 0;  /* Type of rivers */
    int *itype; /* Index of rivertype for each river*/
    double **pstage = NULL; /* Pointers to river stage*/
    double **pflux = NULL; /* Pointers to river flux*/
    
//    double *rivStage;   /*  */
//    double *rivBank;    /* */
//    
//    double *Q_high; /* */
//    double *Q_low;  /* */
//    double *Y_high; /* */
//    double *Y_low;  /* */
    FILE *fid;
    
};
#endif /* FloodAlert_hpp */
