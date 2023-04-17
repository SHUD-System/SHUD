#ifndef functions_h
#define functions_h
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "Macros.hpp"
#include "funPlatform.hpp"
using namespace std;
/* Time Interpolation */
double timeInterp(double t, double t0, double t1, double x0, double x1);

/*==========Screen print function===============*/
double getSecond_cpu(void);
double getSecond_wall(void);

inline
int atInterval(double x, int intv){
    if ((unsigned long )x % intv == 0) {
        return 1;
    } else {
        return 0;
    }
}
/*==========Math function===============*/
double Eudist(double x1, double y1, double x2, double y2);
double min(double a, double b);
double max(double a, double b);
double min(double x[], int n);
double max(double x[], int n);

double mean(double *x, int n);
double stddeviation(double x[], int n);
//double polygonArea(double X[], double Y[], int n);/* Not ready yet.*/

/*==========misc function===============*/
void printDY(char* fn, double * dy, int n, double t); /* Print the dy value out. */
void printDY(char* fn); /* Print dy, create empty file. */
void printVector(FILE *fid, double * x, int xstart, int n, double t);
void printVectorBin(FILE *fid, double * x, int xstart, int n, double t);

void CheckNonZero(double x, int i, const char *s);
void CheckNonNegative(double x, int i, const char *s);
void CheckNonZero(int x, int i, const char *s);
void CheckNonZero(double x, int i, int j, const char *s);
void CheckNANi(double x, int i, const char *s);
void CheckNANij(double x, int i, const char *s);
void CheckNAvalue(double x, const char *s);
void CheckNA(int x, const char *s);
void creatFile(const char *fn);
void creatFile(const char *fn, const char *mode);
void compareVal(double x, double y);
void CheckFile(FILE * fp, const char *s);

int checkRange(int x, int xmin, int xmax, int i, const char *s);
int checkRange(double x, double xmin, double xmax, int i, const char *s);
bool checkexist(int *x, int n, int xkey);
void myexit(int flag);
void screeninfo(const char *s);

template<typename  T>
void screeninfo(const char *s, T x){
    char str[MAXLEN];
    sprintf(str, s, x);
    fprintf(stdout, "%s", str);
}
inline
void setValue (double *x, double val, int start, int n){
    for(int i = 0; i < n; i++){
        x[i + start] = val;
    }
}
inline
void setValue (double *x, int xstart, double *val, int vstart, int n){
    for(int i = 0; i < n; i++){
        x[i + xstart] = val[i + vstart];
    }
}
inline
void Global2Sub(int n1, int n2, int n3){
    setValue(uYsf,  0, globalY, n1 * 0, n1);
    setValue(uYus,  0, globalY, n1 * 1, n1);
    setValue(uYgw,  0, globalY, n1 * 2, n1);
    setValue(uYriv, 0, globalY, n1 * 3, n2);
    setValue(uYlake,0, globalY, n1 * 3 + n2, n3);
}


inline
void Sub2Global(double *x1, double *x2, double *x3, double *x4, double *x5,
                int n1, int n2, int n3){
    setValue(globalY, n1 * 0, x1, 0, n1);
    setValue(globalY, n1 * 1, x2, 0, n1);
    setValue(globalY, n1 * 2, x3, 0, n1);
    setValue(globalY, n1 * 3, x4, 0, n2);
    setValue(globalY, n1 * 3 + n2, x5, 0, n3);
}

inline
void Sub2Global(int n1, int n2, int n3){
    Sub2Global(uYsf, uYus, uYgw, uYriv, uYlake, n1, n2, n3);
}
double Quadratic(double s, double w, double dA);
double fun_dAtodY(double dA, double w0, double s);


inline double Eudist(double x1, double y1, double x2, double y2){
    //Euclidean distance calculator
    double dx = x2 - x1;
    double dy = y2 - y1;
    double d = sqrt(dx * dx + dy * dy);
    return d;
}
inline double min(double a, double b){
    return (a > b ? b : a);
}
inline
double max(double a, double b){
    return (a < b ? b : a);
}

inline double Quadratic(double s, double w, double dA){
    /* dA can be positive or negative. */
    double ret = 0., cc;
    s = fabs(s);
    cc = w * w + 4 * s * dA;
    if( cc < ZERO){
        if( cc < -1. * ZERO){
            fprintf(stderr, "Error in Quadratic\n");
        }
        ret = -1. * w / ( 2. * s );
    }else{
        ret = (-w + sqrt(cc)) / (2 * s);
    }
    return ret;
}

inline double fun_dAtodY(double dA, double w_top, double s){
    /* Delta_A_cs of flux to Dy in river stage @t=t+1*/
    /* w_top = topwidth at t */
    double dy = 0.;
    if(dA == 0.) return 0.;
    
    if( fabs(s) < EPS_SLOPE ){
        dy = dA / w_top;
    }else{
        dy = Quadratic(s, w_top,  dA);
    }
    return dy;
}

inline
int ifequal(double x, double y){
    if( fabs(x - y) < ZERO){
        return 1;
    }else{
        return 0;
    }
}

inline
double dhdx(double *x, double *y, double *h){
    return -1. * (y[2] * (h[1] - h[0]) +
                  y[1] * (h[0] - h[2]) +
                  y[0] * (h[2] - h[1])) /
    (x[2] * (y[1] - y[0]) +
     x[1] * (y[0] - y[2]) +
     x[0] * (y[2] - y[1]));
}
inline
double dhdy(double *x, double *y, double *h){
    return -1. * (x[2] * (h[1] - h[0]) +
                  x[1] * (h[0] - h[2]) +
                  x[0] * (h[2] - h[1])) /
    (y[2] * (x[1] - x[0]) +
     y[1] * (x[0] - x[2]) +
     y[0] * (x[2] - x[1]));
}

inline double fixMaxValue(double x, double defVal){
    if(x < defVal){
        return defVal;
    }else{
        return x;
    }
}

inline double FrozenFraction(double T, double high, double low){
    double x;
    if(T > high){
        return 0; /* All mobile water*/
    }else if(T < low){
        return 1; /* All Frozen water*/
    }else{
        x = (high - T) / (high - low);
        return min(1.0, max(x, 0.0));
    }
}
void PointPerpdicularOnLine(double *xx, double *yy,
                              double x, double y,
                              double x1, double y1,
                            double x2, double y2);
double ZOnLine(double x1, double y1, double z1,
             double x2, double y2, double z2,
               double x3, double y3);
inline
void checkExchangeValue(double **x, int i, int j, int inabr, int jnabr){
    if(inabr >=0 && jnabr >=0){
        if( (x[i][j] + x[inabr][jnabr]) <= ZERO ){
//            printf("1:(%d, %d) - (%d, %d)\t %E  + %E = %E\n",
//                   i, j,
//                   inabr+1, jnabr+1,
//                   x[i][j], x[inabr][jnabr], x[i][j] + x[inabr][jnabr]);
        }else{
            printf("ERROR:(%d, %d) - (%d, %d)\t %E  + %E = %E\n",
                   i, j,
                   inabr+1, jnabr+1,
                   x[i][j], x[inabr][jnabr], x[i][j] + x[inabr][jnabr]);
            myexit(ERRCONSIS);
        }
    }
}
#endif /* functions_h */



