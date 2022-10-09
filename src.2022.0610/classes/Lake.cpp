
#include "Lake.hpp"
Orifice_Outlet::Orifice_Outlet(){}
void Orifice_Outlet::init(double c0, double a){
    C0 = c0;
    Area = a;
}
double Orifice_Outlet::flow(double h0, double h1){
    Hin = h0;
    Hout = h1;
    Q = C0 * Area * sqrt(2 * GRAV * (Hin - Hout) )  * 60.;
    return Q;
}

Weir_Outlet::Weir_Outlet(){}
void Weir_Outlet::init(double cw, double wd){
    Cwr = cw;
    Width = wd;
}
double Weir_Outlet::flow(double hx){
    h = hx;
    Q = Cwr * Width  * sqrt(h) * h;
    return Q;
}

Open_Outlet::Open_Outlet(){};
void Open_Outlet::init(double cw, double wd){
    Cwr = cw;
    Width = wd;
}
double Open_Outlet::flow(double hl, double hr){
    Hlake = hl;
    HRiver = hr;
    double Q = 0.;
    return Q;
}
LakeBathymetry::LakeBathymetry(){}

void LakeBathymetry::read(FILE *fp){
//    TabularData tb;
//    nvalue = tb.read(fp);
//    index = new int[nvalue];
//    yi = new double[nvalue];
//    ai = new double[nvalue];
//    for (int i = 0; i < nvalue; i++){
//        index[i]    = (int) tb.x[i][0];
//        yi[i]     = (double) tb.x[i][1];
//        ai[i]     = (double) tb.x[i][2];
//        printf("%d\t%f\t%f", index[i], yi[i], ai[i]);
//    }
}
LakeBathymetry::~LakeBathymetry(){
    if(nvalue > 0){
        delete index;
        delete yi;
        delete ai;
    }
}
double LakeBathymetry::toparea(double y){
    double ta = ai[0];
    double dy, da;
    if(y <= yi[0]){
        ta = ai[0];
    }else{
        for(int i = 1; i < nvalue; i++){
            if(y < yi[i]){
                da = (ai[i] - ta);
                dy = yi[i] - y;
                ta = da / dy * (y - yi[i-1]) + ta;
                break;
            }else{
                ta = ai[i];
            }
        }
        
    }
    return ta;
}
/*********************************************/
/*   ==========LAKE===========               */
/*********************************************/
_Lake::_Lake(){
}
_Lake::~_Lake(){
    if(NumEleBank > 0){
        delete iEleBank;
        delete QEleSurf;
        delete QEleGW;
    }
    if(NumRivIn > 0){
        delete QRivIn;
    }
    if(NumRivOut > 0){
        delete RivOut;
    }
}
void _Lake::Initialize(){
    QEleGW = new double[NumEleBank];
    QEleSurf = new double[NumEleBank];
    
    QRivIn = new double[NumRivIn];
    QRivOut = new double[NumRivOut];
}
void _Lake::update(){
    u_toparea = bathymetry.toparea(yStage + zmin);
//    CheckNA(u_toparea, "u_toparea in _Lake::update()");
}
void _Lake::readLake(FILE *fp){
    char str[MAXLEN];
    fscanf(fp, "%s %d", str, &NumEleBank);
    iEleBank = new int[NumEleBank];
    for(int i = 0; i < NumEleBank; i++){
        fscanf(fp, "%d", &(iEleBank[i]) );
    }
    
    fscanf(fp, "%s %d", str, &NumRivIn);
    RivIn = new int[NumRivIn];
    for(int i = 0; i < NumRivIn; i++){
        fscanf(fp, "%d", &(RivIn[i]) );
    }
    
    fscanf(fp, "%s %d", str, &NumRivOut);
    RivOut = new int[NumRivOut];
    for(int i = 0; i < NumRivOut; i++){
        fscanf(fp, "%d", &(RivOut[i]) );
    }
}
void _Lake::readBathymetry(FILE *fp){
//    bathymetry.read(fp);
}

void LakeBathymetry::InitValue(int n){
    nvalue = n;
    index = new int[nvalue];
    yi = new double[nvalue];
    ai = new double[nvalue];
    for(int i = 0; i < nvalue; i++){
        index[i] = 0;
        yi[i] = NA_VALUE;
        ai[i] = NA_VALUE;
    }
}
