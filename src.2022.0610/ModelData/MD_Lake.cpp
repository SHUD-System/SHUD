//
//  MD_Lake.cpp
//  shud
//
//  Created by Lele Shu on 2022/4/13.
//  Copyright © 2022 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"


int Model_Data::LakeUniqueID(){
    int *x = new int[NumEle];
    int nl = 0;
    for(int i = 0; i < NumEle; i++){
//        printf("Lake %d\n", Ele[i].iLake);
        if(Ele[i].iLake > 0){
            if( !checkexist(x, nl, Ele[i].iLake) ){
                x[nl] = Ele[i].iLake;
                nl++;
            }else{
                /* void */
            }
        }
    }
    for(int i = 0; i < NumLake; i++){
        printf("Lake %d\n", x[i]);
    }
    delete[] x;
    return nl;
}
void Model_Data::LakeInitialize(){
    int inabr;
    int ii, jj;
    
    lake        = new _Lake[NumLake];
    yLakeStg    = new double[NumLake];
    y2LakeArea  = new double[NumLake];
    QLakeSurf   = new double[NumLake];  //30
    QLakeSub    = new double[NumLake];
    QLakeRivIn  = new double[NumLake];
    QLakeRivOut = new double[NumLake];
    qLakePrcp   = new double[NumLake];
    qLakeEvap   = new double[NumLake];  //34
    
    for(int i = 0; i < NumRiv; i++){
        Riv[i].frLake = -1;
        if(Riv[i].down <= -4){
            lakeon = 1;
            Riv[i].toLake = (-3 - Riv[i].down) - 1;  /* 0 ~ inf */
        }else{
            Riv[i].toLake = NA_VALUE;
        }
    }
    for(int i = 0; i < NumLake; i++){
        lake[i].area    = 0.;
        lake[i].NumEleLake  = 0;
        lake[i].NumEleBank  = 0;
        lake[i].NumRivIn    = 0;
        lake[i].NumRivOut   = 0;
        /* ELEMENT TO LAKE */
        for (int j = 0; j < NumEle ; j++){
            if(Ele[j].iLake == i + 1){ /* Lake element */
                lake[i].NumEleLake++;
                lake[i].area += Ele[i].area;
            }else{
                for(int k = 0; k < 3; k++){
                    inabr = Ele[j].nabr[k] - 1;
                    if(inabr >= 0){
                        if(Ele[inabr].iLake == i + 1  && Ele[j].iLake <= 0 ){
                            /* Neighbor is in Lake */ /* BankElement is not in Lake */
                            lake[i].NumEleBank++;
                        }
                    }
                }
            }
        }
        /* RIVER TO LAKE */
        for (int j = 0; j < NumRiv ; j++){
            if(Riv[j].toLake == i){
//                xrivin[lake[i].NumRivIn] = j;
                lake[i].NumRivIn++;
            }
            if(Riv[j].frLake == i){
//                xrivout[lake[i].NumRivOut] = j;
                lake[i].NumRivOut++;
            }
        }
    }
    
    for(int i = 0; i < NumLake; i++){
        lake[i].iEleLake = new int[lake[i].NumEleLake];
        lake[i].iEleBank = new int[lake[i].NumEleBank];
        lake[i].iRivIn  = new int[lake[i].NumRivIn];
        lake[i].iRivOut = new int[lake[i].NumRivOut];
    }
    
    for(int i = 0; i < NumLake; i++){
        ii = 0;
        jj = 0;
        /* ELEMENT TO LAKE */
        for (int j = 0; j < NumEle ; j++){
            if(Ele[j].iLake == i + 1){ /* Lake element */
                lake[i].iEleLake[ii++] = j;
            }else{
                for(int k = 0; k < 3; k++){
                    inabr = Ele[j].nabr[k] - 1;
                    if(inabr >= 0){ /*  NO BOUNDARY ELEMENT*/
//                        printf("debug:iEle %d(%d),\t inabr.lake =%d \t j.ilake = %d\n", j+1, k+1, Ele[inabr].iLake, Ele[j].iLake);
                        if(Ele[inabr].iLake == i + 1 && Ele[j].iLake <= 0){ /* SELF is not lake element, nabor[k] is lake. */
                            lake[i].iEleBank[jj++] = j;
                        }
                    }
                }
            }
//            printf("\n");
        }
    }
    for(int i = 0; i < NumLake; i++){
        ii = 0;
        jj = 0;
        /* RIVER TO LAKE */
        for (int j = 0; j < NumRiv ; j++){
            if(Riv[j].toLake == i + 1){
                lake[i].iRivIn[ii++] = j;
            }
            if(Riv[j].frLake == i + 1){
                lake[i].iRivOut[jj++] = j;
            }
        }
    }
    
    for(int i = 0; i < NumEle; i++){
        if(Ele[i].iLake<=0){ /* Non-Lake Element */
            for(int j = 0; j < 3; j++){
                inabr = Ele[i].nabr[j] - 1;
                if(inabr >= 0){
                    if(Ele[inabr].iLake>0){ /* Neighbor is a Lake Elememt */
                        Ele[i].lakenabr[j] = Ele[inabr].iLake;
//                        printf("%d(%d) %d : %d\n", i+1, j+1, inabr+1, Ele[i].nabr[j]);
                    }
                }
            }
        }
    }
}
void Model_Data::lake_readBathy(const char *fn){
    TabularData tb;
    int ncol = 0;
    int nrow = 0;
    /*========== open *.bath file ==========*/
    FILE * fp =  fopen(fn, "r");
    CheckFile(fp, fn);
    for(int j = 0; j < NumLake; j++){
        nrow = tb.read(fp, &ncol);
        lake[j].bathymetry.InitValue(nrow);
        for (int i = 0; i < lake[j].bathymetry.nvalue; i++){
            lake[j].bathymetry.index[i] = (int) tb.x[i][0];
            lake[j].bathymetry.yi[i]    = (double) tb.x[i][1];
            lake[j].bathymetry.ai[i]    = (double) tb.x[i][2];
        }
        lake[j].zmin = lake[j].bathymetry.yi[0];
    }
//    for(int j = 0; j < NumLake; j++){
//        printf("Lake %d; \t Nvalue = %d\n", j+1, lake[j].bathymetry.nvalue);
//        for (int i = 0; i < lake[j].bathymetry.nvalue ; i++){
//            printf("%d : %f\t%f\n", lake[j].bathymetry.index[i], lake[j].bathymetry.yi[i], lake[j].bathymetry.ai[i]);
//        }
//    }
    fclose(fp);
}

void Model_Data::LakeTable(const char *fn){
    FILE *fp = fopen(fn, "w");
    fclose(fp);
}
