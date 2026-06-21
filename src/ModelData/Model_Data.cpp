#include "Model_Data.hpp"
#include "is_sm_et.hpp"
#include <cassert> /* S5d.1 (#178) — DEBUG asserts in initialize_hot() */

Model_Data::Model_Data(){
}

Model_Data::Model_Data(FileIn *f_in, FileOut *f_out){
    pf_in = f_in;
    pf_out = f_out;
}
Model_Data::~Model_Data(){
    FreeData();
}
void Model_Data::TimeSpent(){
    /* S1d.2 (openMP #48) — `omp_get_wtime` is the OpenMP wall-clock
     * API; available iff `-fopenmp` was passed (which auto-defines
     * `_OPENMP`). Migrated from the retired legacy macro.
     * Independent of the SHUD-level feature switches. */
#ifdef _OPENMP
    double toc = omp_get_wtime();
    double dt = toc - tic;
    screeninfo("\n\tNumber of calls of f function:\t %ld \n", nFCall);
    printf("\n\tTime used by model:\t %.3f seconds.\n", dt);
    screeninfo("\n\nThe successful end. \n\n");

#else
    clock_t toc = (double)clock();
    double dt = (toc - tic) / CLOCKS_PER_SEC;
    screeninfo("\n\tNumber of calls of f function:\t %ld \n", nFCall);
//    screeninfo("\n\tNumber of calls of f function:\t %ld \n", nFCall2);
    screeninfo("\n\tTime used by model:\t %.3f seconds.\n", dt);
    screeninfo("\n\nThe successful end. \n\n");
#endif

}
void Model_Data::modelSummary(int end){
    char str[MAXLEN];
    screeninfo("\n========================================================\n");
    screeninfo("Summary:\n");
    screeninfo("\tProject name:\t %s\n ", pf_in->projectname);
    screeninfo("\tInput path:\t %s\n", pf_in->inpath);
    screeninfo("\tOutput path:\t %s\n", pf_in->outpath);
    screeninfo("\tCalibration file:\t %s\n", pf_in->file_calib);
    screeninfo("\tParameter file:\t %s\n", pf_in->file_para);
    screeninfo("\tModel starts at: %.2f day\n", CS.StartTime / 1440);
    screeninfo("\tModel ends at: %.2f day\n", CS.EndTime / 1440);
    screeninfo("\tModel time step(max): %.2f minutes\n", CS.MaxStep);
//    screeninfo("\tModel time step(normal): %.2f minutes\n", CS.SolverStep);
    screeninfo("\tModel total number of steps(minimum): %d \n", CS.NumSteps);
    sprintf(str,"\tSize of model: \tNcell = %d \tNriver = %d\t NSeg = %d", NumEle, NumRiv, NumSegmt);
    screeninfo(str);
    /* S1d.2 (openMP #48) — same migration as TimeSpent() above. */
#ifdef _OPENMP
    screeninfo("\n\n\tOpenMP enable. No of threads = %d\n", CS.num_threads);
    screeninfo("\n========================================================\n");
    if (end) {
        TimeSpent();
    } else {
        tic = omp_get_wtime();
        nFCall = 0;
        screeninfo("\nModel Starting ... \n\n");
    }
#else
    screeninfo("\n\tOpenMP disable");
    screeninfo("\n========================================================\n");
    if (end) {
        TimeSpent();
    } else {
        tic = (double)clock();
        nFCall = 0;
//        nFCall2 = 0;
        screeninfo("\nModel Starting ... \n\n");
    }
#endif
}

void Model_Data::malloc_Y(){
    NumY = 3 * NumEle + 1 * NumRiv + 1 * NumLake;
    globalY = new double[NumY];
}

void Model_Data::malloc_EleRiv(){
    
    /* allocate memory storage to flux terms */
    QeleSurf    = new double *[NumEle];
    QeleSub     = new double *[NumEle];
    QeleSurfTot = new double[NumEle];
    QeleSubTot  = new double[NumEle];
    QoutSurf    = new double[NumEle]; // 5
    
    Qe2r_Surf = new double[NumEle]; //5.1
    Qe2r_Sub  = new double[NumEle]; // 5.2
    /* S3b (PR-9): per-edge slots used only when lakeon. Always allocate
     * (NumEle-sized; cheap) so non-lake builds need no conditional
     * cleanup. Touched only inside the lake branches of fun_Ele_surface
     * / fun_Ele_sub and gathered in PassValue_legacy. */
    QeleSurf_lake = new double[NumEle * 3];
    QeleSub_lake  = new double[NumEle * 3];
    qEleEvapo_lake = new double[NumEle]; // S3b.4 (PR-9)
    qElePrep_lake  = new double[NumEle]; // S3b.4 (PR-9)

    qEleE_IC      = new double[NumEle];
    qEleEvapo      = new double[NumEle];
    qEleTrans      = new double[NumEle];
    qElePrep    = new double[NumEle];
    qEleTF      = new double[NumEle];
    qEleETP     = new double[NumEle];
    qPotEvap     = new double[NumEle];
    qPotTran    = new double[NumEle];
    qEleETA     = new double[NumEle];
    qEs     = new double[NumEle];
    qEu     = new double[NumEle];
    qEg     = new double[NumEle];
    qTu     = new double[NumEle];
    qTg     = new double[NumEle];
    
    qEleETloss  = new double[NumEle]; //10
    iBeta     = new double[NumEle]; // 10.1
    
    qEleNetPrep = new double[NumEle];
    qEleInfil   = new double[NumEle];
    qEleExfil   = new double[NumEle];
    qEleRecharge = new double[NumEle]; //13
    
    yEleIS      = new double[NumEle];
    yEleISmax   = new double[NumEle];
    yEleISsnowmax = new double[NumEle]; //16
    
    yEleSnow    = new double[NumEle];
    yEleSnowGrnd = new double[NumEle];
    yEleSnowCanopy = new double[NumEle];
    yEleGW      = new double[NumEle];    //20
    
    yEleSurf    = new double[NumEle];
    yEleUnsat   = new double[NumEle];
    yEleWetFront = new double[NumEle];  //23
    
    fu_Surf     = new double[NumEle];
    fu_Sub      = new double[NumEle];
    AccT_surf   = new _AccTemp[NumEle];
    AccT_sub    = new _AccTemp[NumEle];
    
    yRivStg     = new double[NumRiv];
    QrivSurf    = new double[NumRiv];
    QrivSub     = new double[NumRiv];   //26
    QrivDown    = new double[NumRiv];
    QrivUp      = new double[NumRiv]; //28
    QsegSurf    = new double[NumSegmt];
    QsegSub     = new double[NumSegmt];
    
    uYsf = new double[NumEle];
    uYus = new double[NumEle];
    uYgw = new double[NumEle];
//    uYele = new double[NumY1];  // 35
    
    if(NumRiv > 0){
        uYriv = new double[NumRiv];  // 35.1
    }
    
    for (int i = 0; i < NumEle; i++) {
        QeleSurf[i] = new double[3];
        QeleSub[i] = new double[3];
    }
    
    t_prcp  = new double[NumEle];  //
    t_temp  = new double[NumEle];  //
    t_rh    = new double[NumEle];  //
    t_sph    = new double[NumEle];  //
    t_wind  = new double[NumEle];  //
    t_rn    = new double[NumEle];  //
    t_vp    = new double[NumEle];  //
    t_lai   = new double[NumEle];  //
    t_mf    = new double[NumEle];  //
//    t_hc    = new double[NumEle];  //

    /* S5d.1 (#178) — ElementHotData SoA allocation. Sized NumEle (or
     * NumEle*3 for flat-3 arrays). Layout mirrors docs/s5d_hot_fields.yaml.
     * NumEle is read at this call site; if NumEle changes after this
     * call, hot must be reallocated (no such code path exists today).
     * Free in symmetric order at MD_readin.cpp Model_Data::FreeData(). */
    hot.nabr_flat       = new int[NumEle * 3];
    hot.lakenabr_flat   = new int[NumEle * 3];
    hot.edge_flat       = new double[NumEle * 3];
    hot.area            = new double[NumEle];
    hot.z_bottom        = new double[NumEle];
    hot.z_surf          = new double[NumEle];
    hot.iSoil           = new int[NumEle];
    hot.iLC             = new int[NumEle];
    hot.iMF             = new int[NumEle];
    hot.iForc           = new int[NumEle];
    hot.iLake           = new int[NumEle];
    hot.iBC             = new int[NumEle];
    hot.iSS             = new int[NumEle];
    hot.Dist2Nabor_flat = new double[NumEle * 3];
    hot.Dist2Edge_flat  = new double[NumEle * 3];
    hot.avgRough_flat   = new double[NumEle * 3];
    hot.FixPressure     = new double[NumEle];
    hot.WetlandLevel    = new double[NumEle];
    hot.RootReachLevel  = new double[NumEle];
    hot.depression      = new double[NumEle];
    hot.QBC             = new double[NumEle];
    hot.QSS             = new double[NumEle];
    hot.windH           = new double[NumEle];
    hot.u_qi            = new double[NumEle];
    hot.u_qex           = new double[NumEle];
    hot.u_effKH         = new double[NumEle];
    hot.u_satn          = new double[NumEle];
    hot.Sy              = new double[NumEle];
    hot.VegFrac         = new double[NumEle];
    hot.Albedo          = new double[NumEle];
    hot.Rough           = new double[NumEle];
    hot.ImpAF           = new double[NumEle];
}

void Model_Data::initialize_hot() {
    /* S5d.1 (#178) — populate ElementHotData SoA from _Element AoS.
     * Bitwise contract: every SoA value matches the AoS source EXACTLY
     * (assignment-only; no rounding or cast loss). Called from
     * Model_Data::initialize() AFTER element AoS load is complete and
     * BEFORE any RHS dispatch. Dynamic fields (u_qi, u_qex, u_effKH,
     * u_satn) are seeded here from current Ele[i] values; subsequent
     * writes via Ele[i].updateElement / Flux_Infiltration / etc. are
     * propagated by sync_hot_dynamic(i) at each call site. */
    for (int i = 0; i < NumEle; ++i) {
        for (int j = 0; j < 3; ++j) {
            hot.nabr_flat[3*i + j]       = Ele[i].nabr[j];
            hot.lakenabr_flat[3*i + j]   = Ele[i].lakenabr[j];
            hot.edge_flat[3*i + j]       = Ele[i].edge[j];
            hot.Dist2Nabor_flat[3*i + j] = Ele[i].Dist2Nabor[j];
            hot.Dist2Edge_flat[3*i + j]  = Ele[i].Dist2Edge[j];
            hot.avgRough_flat[3*i + j]   = Ele[i].avgRough[j];
        }
        hot.area[i]           = Ele[i].area;
        hot.z_bottom[i]       = Ele[i].z_bottom;
        hot.z_surf[i]         = Ele[i].z_surf;
        hot.iSoil[i]          = Ele[i].iSoil;
        hot.iLC[i]            = Ele[i].iLC;
        hot.iMF[i]            = Ele[i].iMF;
        hot.iForc[i]          = Ele[i].iForc;
        hot.iLake[i]          = Ele[i].iLake;
        hot.iBC[i]            = Ele[i].iBC;
        hot.iSS[i]            = Ele[i].iSS;
        hot.FixPressure[i]    = Ele[i].FixPressure;
        hot.WetlandLevel[i]   = Ele[i].WetlandLevel;
        hot.RootReachLevel[i] = Ele[i].RootReachLevel;
        hot.depression[i]     = Ele[i].depression;
        hot.QBC[i]            = Ele[i].QBC;
        hot.QSS[i]            = Ele[i].QSS;
        hot.windH[i]          = Ele[i].windH;
        hot.u_qi[i]           = Ele[i].u_qi;
        hot.u_qex[i]          = Ele[i].u_qex;
        hot.u_effKH[i]        = Ele[i].u_effKH;
        hot.u_satn[i]         = Ele[i].u_satn;
        hot.Sy[i]             = Ele[i].Sy;
        hot.VegFrac[i]        = Ele[i].VegFrac;
        hot.Albedo[i]         = Ele[i].Albedo;
        hot.Rough[i]          = Ele[i].Rough;
        hot.ImpAF[i]          = Ele[i].ImpAF;

#ifdef DEBUG
        /* S5d.1 (#178) — sample assertion to catch SoA-vs-AoS drift on
         * DEBUG builds. Spec: scenario "DEBUG 一致性 assertion 通过".
         * Sampling = full sweep across all elements (cheap; DEBUG only). */
        assert(hot.area[i]    == Ele[i].area);
        assert(hot.u_effKH[i] == Ele[i].u_effKH);
        assert(hot.iLake[i]   == Ele[i].iLake);
        assert(hot.VegFrac[i] == Ele[i].VegFrac);
        assert(hot.Sy[i]      == Ele[i].Sy);
        for (int j = 0; j < 3; ++j) {
            assert(hot.nabr_flat[3*i+j] == Ele[i].nabr[j]);
            assert(hot.edge_flat[3*i+j] == Ele[i].edge[j]);
        }
#endif
    }
}

void Model_Data::copyCalib(){
    for (int i = 0; i < NumSoil; i++) {
        Soil[i].applyCalib(&(gc.csoil));
    }
    for (int i = 0; i < NumGeol; i++) {
        Geol[i].applyCalib(&(gc.cgeol));
    }
    for (int i = 0; i < NumLC; i++) {
        LandC[i].applyCalib(&(gc.clandc));
    }
    for (int i = 0; i < NumRivType; i++) {
        Riv_Type[i].applyCalib(&(gc.criv));
    }
}
void Model_Data::InitFloodAlert(const char *fn){
    flood->InitAlert(NumRiv, NumRivType);
    flood->InitPointer(yRivStg, QrivDown);
    flood->InitPara(Riv_Type);
    for(int i =  0; i < NumRiv; i++){
        flood->pushRiverType(i, Riv[i].type);
    }
    flood->InitFile(fn);
}
double Model_Data::updateArea(){
    WatershedArea = 0.;
    for(int i = 0; i < NumEle; i++){
        WatershedArea += Ele[i].area;
    }
//    for(int i = 0; i < NumLake; i++){
//        WatershedArea += Lake[i].area;
//    }
    return WatershedArea;
}
double Model_Data::getArea(){
    return WatershedArea;
}
void Model_Data::rmSinks(){
    double zmin_nb, zmax_current;
    int inabr;
    for(int i = 0; i < NumEle; i++){
        zmax_current =  Ele[i].z_surf;
        zmin_nb = 1.0e200;
        for(int j = 0; j < 3; j++){
            inabr = Ele[i].nabr[j] - 1;
            if(inabr >= 0){ /* Nabr exists */
                zmin_nb = min(zmin_nb, Ele[inabr].z_surf);
            }
        }
        if( zmin_nb > zmax_current){
            if(Ele[i].RivID <= 0) {
                fprintf(stderr, "Warning: remove sink on %d, from %.2f to %.2f. dz = %.2f\n", i+1, zmax_current, zmin_nb, zmin_nb - zmax_current);
                Ele[i].z_surf = zmin_nb;
                Ele[i].z_bottom = zmin_nb - Ele[i].AquiferDepth;
                
            }else{
                /* Void*/
            }
        }
        
    }
    
    for (int i = 0; i < NumEle; i++) {
        Ele[i].InitElement();
    }
}

void Model_Data::debugData(const char *outdir){
    char fn[MAXLEN];
    char str[MAXLEN];
    sprintf(str, "%s/%s", outdir, "Debug_Table");
    sprintf(file_debug, "%s/%s", outdir, "DY.dat");
    if(NumEle > 0){
        sprintf(fn, "%s%s", str, "_Element.csv");
        ElementTable(fn);
    }
    if(NumRiv > 0){
        sprintf(fn, "%s%s", str, "_River.csv");
        RiverTable(fn);
    }
    if(NumLake > 0){
        sprintf(fn, "%s%s", str, "_Lake.csv");
        LakeTable(fn);
    }
}
void Model_Data::ElementTable(const char *fn){
    FILE *fp = fopen(fn, "w");
    Ele[0].printHeader(fp);
    for(int i = 0; i < NumEle; i++){
        Ele[i].printInfo(fp);
    }
    fclose(fp);
}
void Model_Data::RiverTable(const char *fn){
    FILE *fp = fopen(fn, "w");
    Riv[0].printHeader(fp);
    for(int i = 0; i < NumRiv; i++){
        Riv[i].printInfo(fp);
    }
    fclose(fp);
}
int Model_Data::ScreenPrintu(double t, unsigned long it){
    int flag = 0;
#ifdef DEBUG
    printf("%.0f min ~ %.4f day\t %.2f%% \n", t, t / 1440., (double)it / CS.NumSteps * 100 );
    flag = 1;
#else
    static double tnext = t;
    static unsigned long ncall1 = 0, ncall2 = 0, ncall3 = 0, ncall4 = 0, ncall5 = 0;
    if (t >= tnext) {
        printf("%6.2f d \t %5.2f%% \t %6.2f s \t %6ld %6ld %6ld %6ld %6ld\n",
               t / 1440, 100.0 * it / CS.NumSteps, getSecond_wall(),
               nFCall1 - ncall1, nFCall2 - ncall2, nFCall3 - ncall3, nFCall4 - ncall4, nFCall5 - ncall5
               );
        tnext += CS.screenIntv;
        ncall1 = nFCall1;
        ncall2 = nFCall2;
        ncall3 = nFCall3;
        ncall4 = nFCall4;
        ncall5 = nFCall5;
        flag = 1;
    }
#endif
    return flag;
}
int Model_Data::ScreenPrint(double t, unsigned long it){
    int flag = 0;
#ifdef DEBUG
    printf("%.0f min ~ %.4f day\t %.2f%% \n", t, t / 1440., (double)it / CS.NumSteps * 100 );
    flag = 1;
#else
    static double tnext = t;
    static unsigned long ncall = 0;
    double sec_cpu, sec_wall, Perctage;
    if (t >= tnext) {
        sec_cpu = getSecond_cpu();
        sec_wall = getSecond_wall();
        Perctage = 100.0 * it / CS.NumSteps;
        printf("%.2f day \t %.2f%% \t %.2f s \t %.2f s \t %ld \n", tnext / 1440, Perctage, sec_cpu, sec_wall, nFCall - ncall);
        pf_out->writeTime(t, Perctage, sec_cpu, sec_wall, nFCall - ncall);
        tnext += CS.screenIntv;
        ncall = nFCall;
        flag = 1;
    }
#endif
    return flag;
}

