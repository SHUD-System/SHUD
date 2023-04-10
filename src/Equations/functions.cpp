//
//functions.cpp
//
//Created by Lele Shu on 6 / 23 / 18.
// Copyright Â ©2018 Lele Shu.All rights reserved.
//

#include "functions.hpp"

void myexit(int flag){
    switch (flag) {
        case ERRFileIO:
            fprintf(stderr, "\nEXIT with error code %d(FILEIO)\n", flag);
            break;
        case ERRNAN:
            fprintf(stderr, "\nEXIT with error code %d(NAN/INF VALUE)\n", flag);
            break;
        case ERRCVODE:
            fprintf(stderr, "\nEXIT with error code %d(CVODE)\n", flag);
            break;
        case ERRCONSIS:
            fprintf(stderr, "\nEXIT with error code %d(Data Consistency)\n", flag);
            break;
        case ERRDATAIN:
            fprintf(stderr, "\nEXIT with error code %d(Data validation)\n", flag);
            break;
        case ERRSUCCESS:
            fprintf(stderr, "\nSuccess.\n");
            break;
        default:
            fprintf(stderr, "\nEXIT with error code %d(Undefined error)\n", flag);
            break;
    }
    fprintf(stderr, "\n\n\n");
    exit(flag);
}
bool checkexist(int *x, int n, int xkey){
    for(int i = 0; i < n; i++){
        if(x[i] == xkey){
            return 1;
        }
    }
    return 0;
}
int checkRange(int x, int xmin, int xmax, int i, const char *s){
    if(x < xmin || x > xmax){
        fprintf(stderr, "\nInvalid value %d for %s(%d) is out of range [%d, %d]\n", x, s, i+1,  xmin, xmax);
        return 1;
    }else{
        return 0;
    }
}
int checkRange(double x, double xmin, double xmax, int i, const char *s){
    if(x < xmin || x > xmax){
        fprintf(stderr, "\nInvalid value %f for %s(%d) is out of range [%g, %g]\n", x, s, i+1,  xmin, xmax);
        return 1;
    }else{
        return 0;
    }
}
double timeInterp(double t, double t0, double t1, double x0, double x1){
    double dt = t1 - t0;
    if (fabs(dt) < ZERO) {
        return x1;
    } else {
        return ((t1 - t) * x0 + (t - t0) * x1) / dt;
    }
}

double min(double *x, int n){
    double  ret;
    ret = x[0];
    for(int i = 0; i < n; i++){
        if(ret > x[i]){
            ret = x[i];
        }
    }
    return ret;
}
double max(double *x, int n){
    double  ret;
    ret = x[0];
    for(int i = 0; i < n; i++){
        if(ret < x[i]){
            ret = x[i];
        }
    }
    return ret;
}
void CheckNANi(double x, int i, const char *s)
{
    if (isnan(x) || isinf(x)) {
        printf("\nERROR: NAN error for %s %d\n", s, i + 1);
        myexit(ERRNAN);
    }
}
void CheckNANij(double x, int i, const char *s)
{
    if (isnan(x) || isinf(x)) {
        printf("\nERROR: NAN error for %s %d\n", s, i + 1);
        myexit(ERRNAN);
    }
}
double getSecond_wall(void){
#ifdef _OPENMP_ON
    static double t0 = 0.;
    double t1;
    double sec;
    t1 = omp_get_wtime();
    sec = (double)(t1 - t0);
    t0 = t1;
#else
    double sec;
//    clock_t t1;
//    static clock_t t0 = 0;
//    t1 = clock();
//    sec = (double)(t1 - t0) / CLOCKS_PER_SEC;
    static double t0 = get_wall_time();
    double t1;
    t1 = get_wall_time();
    sec = (double)(t1 - t0);
    t0 = t1;
#endif
    return sec;
}
double getSecond_cpu(void){
#ifdef _OPENMP_ON
    static double t0 = 0.;
    double t1;
    double sec;
    t1 = omp_get_wtime();
    sec = (double)(t1 - t0);
    t0 = t1;
#else
    double sec;
//    clock_t t1;
//    static clock_t t0 = 0;
//    t1 = clock();
//    sec = (double)(t1 - t0) / CLOCKS_PER_SEC;
    static double t0 = get_cpu_time();
    double t1;
    t1 = get_cpu_time();
    sec = (double)(t1 - t0);
    t0 = t1;
#endif
    return sec;
}
void CheckNonNegative(double x, int i, const char *s)
{
    if (x < 0.0 || isnan(x) || isinf(x) || fabs(x - NA_VALUE) < ZERO) {
        printf("ERROR: Negative Value %e for %s of Element %d is not allowed. Please check again.\n", x, s, i + 1);
        myexit(ERRNAN);
    }
}
void CheckNonZero(double x, int i, const char *s)
{
    if (x <= 0.0 || isnan(x) || isinf(x) || fabs(x - NA_VALUE) < ZERO) {
        printf("ERROR: Value %e for %s of Element %d is not allowed. Please check again.\n", x, s, i + 1);
        myexit(ERRNAN);
    }
}
void CheckNonZero(int x, int i, const char *s)
{
    if (x <= 0) {
        printf("ERROR: Value %d for %s of Element %d is not allowed. Please check again.\n", x, s, i + 1);
        myexit(ERRNAN);
    }
}
void compareVal(double x, double y)
{
    double d = fabs(x - y);
    if (d != 0.) {
        if (d > ZERO) {
            //if (d / fabs(x) > 1.e-4 || d / fabs(y) > 1.e-4) {
            printf("\nValues are different\n\n");
            printf("%f \t %f \t %e\n\n", x, y, d);
            printf("\n");
        } else {
            d = d;
            //printf("\nValues are different, due to truncation error.\n\n");
            //printf("%f \t %f, error=%e\n\n", x, y, d);
        }
    }
}

double mean(double x[], int n){
    double y = x[0];
    for (int i = 1; i < n; i++) {
        y += x[i];
    }
    return y / n;
}

double stddeviation(double x[], int n){
    double ret = 0.;
    double xmean;
    xmean = mean(x, n);
    for (int i = 0; i < n; ++i) {
        ret += (x[i] - xmean) * (x[i] - xmean);
    }
    return sqrt(ret / n);
}
void CheckFile(FILE * fp, const char *s){
    if (fp == NULL) {
        fprintf(stderr, "\n  Fatal Error: \n %s is in use or does not exist!\n", s);
        myexit(ERRFileIO);
    }
}

void CheckNA(double x, const char *s){
    if (fabs(x - NA_VALUE) < ZERO) {
        fprintf(stderr, "\n  Fatal Error: %s is NA_value \n", s);
        myexit(ERRDATAIN);
    }
}
void CheckNA(int x, const char *s){
    if (x == NA_VALUE || fabs(x) > 1e100 ) {
        fprintf(stderr, "\n  Fatal Error: %s is NA_value \n", s);
        myexit(ERRDATAIN);
    }
}
void creatFile(const char *fn){
    creatFile(fn, "w");
}
void creatFile(const char *fn, const char *mode){
    FILE *fp = fopen(fn, mode);
    CheckFile(fp, fn);
    fclose(fp);
}
void screeninfo(const char *s){
    screeninfo(s,"");
}

void printDY(char* fn, double * dy, int n, double t){
    /* Print the dy value out. */
    FILE *file_debug = fopen(fn, "ab");
    printVectorBin(file_debug, dy, 0, n, t);
    fclose(file_debug);
}
void printDY(char* fn){
    /* Print dy, create empty file. */
    FILE *file_debug = fopen(fn, "wb");
    fclose(file_debug);
}
void printVectorBin(FILE *fid, double * x, int xstart, int n, double t){
    fwrite (&t, sizeof (double), 1, fid);
    fwrite (x, sizeof (double), n, fid);
    fflush (fid);
}
void printVector(FILE *fid, double * x, int xstart, int n, double t){
    fprintf(fid, "%f\t", t);;
    for(int i = 0; i < n; i++){
        fprintf(fid, "%.3e\t", x[i + xstart]);
    }
    fprintf(fid, "\n");;
}


void PointPerpdicularOnLine(double *xx, double *yy,
                            double x, double y,
                            double x1, double y1,
                            double x2, double y2){
    // https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;
    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = -1.;
    
    if (len_sq != 0){ //in case of 0 length line
        param = dot / len_sq;
    }else{
        param = -1.;
    }
    
    if (param < 0.) {
        *xx = x1;
        *yy = y1;
    }else if (param > 1.) {
        *xx = x2;
        *yy = y2;
    }else{
        *xx = x1 + param * C;
        *yy = y1 + param * D;
    }
}
double ZOnLine(double x1, double y1, double z1,
               double x2, double y2, double z2,
               double x3, double y3){
    double D = Eudist(x1, y1, x2, y2);
    double dx = Eudist(x1, y1, x3, y3);
    double dz = z2 - z1;
    return z1 + dz / D  * dx;
}

