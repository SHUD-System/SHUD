//  Model_Control.cpp
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Model_Control.hpp"
void PrintOutDt::defaultmode(){
    int dt = 1440;
    /* Element storage */
    dt_ye_gw = dt;
    dt_ye_surf = dt;
    dt_ye_snow = dt;
    dt_ye_ic = 0;
    dt_ye_unsat = dt;
    
    /* Element Fluxes */
    dt_qe_prcp = dt; // default output PRCP.
    dt_qe_infil = dt;
    dt_qe_et = dt;
    
    dt_qe_rech = dt;
    dt_qe_etp = dt;
    dt_qe_eta = dt;
    
    /* Element volume Fluxes */
    dt_Qe_sub = 0;
    dt_Qe_surf = 0;
    
    /* River Stage */
    dt_yr_stage = dt;
    /* River volume Fluxes */
    dt_Qr_up = dt;
    dt_Qr_down = dt;
    dt_Qr_sub = dt;
    dt_Qr_surf = dt;
    
    /* Lake  stage */
    dt_lake     = dt;
}
void PrintOutDt::calibmode(int dt ){
    /* Element storage */
    dt_ye_gw = dt;
    dt_ye_surf = dt;
    dt_ye_snow = 0;
    dt_ye_ic = 0;
    dt_ye_unsat = 0;
    
    /* Element Fluxes */
    dt_qe_prcp = dt; // default output PRCP.
    dt_qe_infil = 0;
    dt_qe_et = dt;
    
    dt_qe_rech = 0;
    dt_qe_etp = dt;
    dt_qe_eta = dt;
    
    /* Element volume Fluxes */
    dt_Qe_sub = 0;
    dt_Qe_surf = 0;
    
    /* River Stage */
    dt_yr_stage = dt;
    /* River volume Fluxes */
    dt_Qr_down = dt;
    dt_Qr_sub = dt;
    dt_Qr_surf = dt;
    
    /* Lake  stage */
    dt_lake     = dt;
}
Control_Data::Control_Data(){
}
Control_Data::~Control_Data(){
//    delete Tout;
}
void Control_Data::ExportResults(double t){
    for (int i = 0; i < NumPrint; i++){
        PCtrl[i].PrintData(dt, t);
    }
}

void Control_Data::updateSimPeriod(){
    updateSimPeriod(DayStart, DayEnd);
}
void Control_Data::updateSimPeriod(double day0, double day1){
    DayStart = day0;
    DayEnd = day1;
    StartTime = DayStart  * 1440;
    EndTime = (DayEnd) * 1440;
    NumSteps = (unsigned long) (EndTime - StartTime) / SolverStep;
}


void Control_Data::read(const char *fn){
    char    str[MAXLEN];
    char    optstr[MAXLEN];
#ifdef _OPENMP_ON
    num_threads = omp_get_max_threads();; /*Default number of threads for OpenMP*/
#else
    num_threads = 0;
#endif
    
    FILE *fp;
    fp = fopen(fn, "r");
    CheckFile(fp, fn);
    /* Read through parameter file to find parameters */
    double val;
    while (fgets(str, MAXLEN, fp)) {
        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0' || str[0] == ' ')
        {
            continue;
        }
        sscanf (str, "%s %lf", optstr, &val);
        /* Get Model Parameters */
        if (strcasecmp ("ASCII_OUTPUT", optstr) == 0)
            Ascii=  (int) val;
        else if (strcasecmp ("BINARY_OUTPUT", optstr) == 0)
            Binary = (int)  val;
        else if (strcasecmp ("SpinupDay", optstr) == 0)
            Spinup =  (int) val;
        else if (strcasecmp ("SCR_INTV", optstr) == 0)
            screenIntv = (int)  val;
        else if (strcasecmp ("VERBOSE", optstr) == 0)
            Verbose = (int)  val;
        else if (strcasecmp ("CloseBoundary", optstr) == 0)
            CloseBoundary = (int)  val;
        else if (strcasecmp ("INIT_MODE", optstr) == 0)
            init_type = (int)  val;
        else if (strcasecmp ("NUM_OPENMP", optstr) == 0)
            num_threads = (int)  val;
        else if (strcasecmp ("ABSTOL", optstr) == 0)
            abstol =  val;
        else if (strcasecmp ("RELTOL", optstr) == 0)
            reltol =  val;
        else if (strcasecmp ("INIT_SOLVER_STEP", optstr) == 0)
            InitStep =  val;
        else if (strcasecmp ("MAX_SOLVER_STEP", optstr) == 0)
            MaxStep =  val;
        else if (strcasecmp ("Update_IC_STEP", optstr) == 0)
            UpdateICStep =  val;
        else if (strcasecmp ("ET_Mode", optstr) == 0)
            ET_Mode =  val;
        else if (strcasecmp ("ET_STEP", optstr) == 0 || strcasecmp ("LSM_STEP", optstr) == 0)
            ETStep =  val;
        else if (strcasecmp ("START", optstr) == 0)
            DayStart =  val;/* Convert days to Minutes */
        else if (strcasecmp ("END", optstr) == 0)
            DayEnd =  val;    /* Convert days to Minutes */
        else if (strcasecmp ("Exfiltration", optstr) == 0)
            exfiltration =  val;
        else if (strcasecmp ("cryosphere", optstr) == 0)
            cryosphere =  val;
//        else if (strcasecmp ("STEPSIZE_FACTOR", optstr) == 0)
//            a =  val;
//        else if (strcasecmp ("MODEL_STEPSIZE", optstr) == 0)
//            b =  val;
        
        /* Element print out control */
        // Y ele
        else if (strcasecmp ("dt_ye_ic", optstr) == 0)
            dt_ye_ic =  val;
        else if (strcasecmp ("dt_ye_SNOW", optstr) == 0)
            dt_ye_snow =  val;
        else if (strcasecmp ("dt_ye_SURF", optstr) == 0)
            dt_ye_surf =  val;
        else if (strcasecmp ("dt_ye_UNSAT", optstr) == 0)
            dt_ye_unsat =  val;
        else if (strcasecmp ("dt_ye_GW", optstr) == 0)
            dt_ye_gw =  val;
        // q Ele
        else if (strcasecmp ("dt_qe_PRCP", optstr) == 0)
            dt_qe_prcp =  val;
        else if (strcasecmp ("dt_qe_ET", optstr) == 0){
            dt_qe_et  =  val;
            dt_qe_etp =  val;
            dt_qe_eta =  val;
        }
        else if (strcasecmp ("dt_qe_rech", optstr) == 0)
            dt_qe_rech =  val;
        else if (strcasecmp ("dt_qe_infil", optstr) == 0)
            dt_qe_infil =  val;
        // Q Ele
        else if (strcasecmp ("dt_Qe_sub", optstr) == 0)
            dt_Qe_sub =  val;
        else if (strcasecmp ("dt_Qe_subx", optstr) == 0)
            dt_Qe_subx =  val;
        else if (strcasecmp ("dt_Qe_surf", optstr) == 0)
            dt_Qe_surf =  val;
        else if (strcasecmp ("dt_Qe_surfx", optstr) == 0)
            dt_Qe_surfx =  val;
        else if (strcasecmp ("dt_Qe_rsub", optstr) == 0)
            dt_Qe_rsub =  val;
        else if (strcasecmp ("dt_Qe_rsurf", optstr) == 0)
            dt_Qe_rsurf =  val;
        
        /* River print out control */
        // y Riv
        else if (strcasecmp ("dt_yr_stage", optstr) == 0)
            dt_yr_stage =  val;
        // Q riv
        else if (strcasecmp ("dt_Qr_Surf", optstr) == 0)
            dt_Qr_surf=  val;
        else if (strcasecmp ("dt_Qr_Sub", optstr) == 0)
            dt_Qr_sub =  val;
        else if (strcasecmp ("dt_Qr_down", optstr) == 0)
            dt_Qr_down =  val;
        else if (strcasecmp ("dt_Qr_up", optstr) == 0)
            dt_Qr_up =  val;
        /* Lake print out control */
        // Y lake
        else if (strcasecmp ("dt_lake", optstr) == 0)
            dt_lake =  val;
        /* Unrecognized Parameter Flag */
        else{
            printf
            ("\n  Parameter:%s cannot be recognized. Please see User's Manual for more details!\n",
             optstr);
//            printf("Continue? (y/N)\n");
//            char cc = getchar();
//            if(cc =='Y' || cc=='y' ){
//                /* Void */
//            }else{
//                myexit(ERRFileIO);
//            }
        }
    }
    SolverStep = MaxStep; // Improve it in future for overcasting.
    fclose (fp);
    updateSimPeriod();
}
void Control_Data::write(const char *fn){
    
}
void Control_Data::getValue(const char *varname){
    
}
Print_Ctrl::Print_Ctrl(){}
void Print_Ctrl::setHeader(const char *s){
    strcpy(header, s);
}
void Print_Ctrl::open_file(int a, int b){
    Ascii = a;
    Binary = b;
    double tmp;
    if(strlen(filename) < 1){
        fprintf(stderr, "WARNING: filename (%s)is empty.\n;", filename);
    }
    sprintf(filea, "%s.csv", filename);
    sprintf(fileb, "%s.dat", filename);
    if (Binary){
        fid_bin = fopen (fileb, "wb");
        CheckFile(fid_bin, fileb);
        fwrite(header, sizeof(char), 1024, fid_bin);
        tmp = (double) StartTime;
        fwrite( &tmp, sizeof(tmp), 1, fid_bin);
        tmp = (double) NumVar;
        fwrite( &tmp, sizeof(tmp), 1, fid_bin);
        fwrite( icol, sizeof(double), NumVar, fid_bin);
        if(global_fflush_mode){
            fflush(fid_bin);
        }
    }
    if (Ascii){
        fid_asc = fopen (filea, "w");
        CheckFile(fid_asc, filea);
        fprintf(fid_asc, "%d\t %d\t %ld\n", 0, NumVar, StartTime);
        fprintf(fid_asc, "%s", "Time_min");
        for(int i = 0; i < NumVar; i++){
            fprintf(fid_asc, " \tX%d", i + 1);
        }
        fprintf(fid_asc, "\n");
        if(global_fflush_mode){
            fflush(fid_asc);
        }
    }
}
void Print_Ctrl::Init(long st, int n, const char *s, int dt, double *x, int iFlux){
    StartTime = st;
    NumVar  = n;
    PrintVar = new double*[NumVar];
    buffer  = new double[NumVar];
    icol    = new double[NumVar];
    strcpy(filename, s);
    if(strlen(filename) < 1){
        fprintf(stderr, "WARNING: filename (%s)is empty.\n;", filename);
    }
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    Interval = dt;
    for(int i=0; i<NumVar; i++){
        icol[i] = (double) (i + 1);
        PrintVar[i] = &x[i];
//        *(PrintVar[i]) = 0.0;
        buffer[i] = 0.0;
    }
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1.;
    }
}
void Print_Ctrl::InitIJ(long st, int n, const char *s, int dt, double **x, int j, int iFlux){
    StartTime = st;
    NumVar  = n;
    PrintVar = new double*[NumVar];
    buffer  = new double[NumVar];
    icol    = new double[NumVar];
    strcpy(filename, s);
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    Interval = dt;
    for(int i=0; i<NumVar; i++){
        icol[i] = (double) (i + 1);
        PrintVar[i] = &(x[i][j]);
//        *(PrintVar[i]) = 0.0;
        buffer[i] = 0.0;
    }
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1;
    }
}

void Print_Ctrl::Init(long st, int n, const char *s, int dt, double *x, int iFlux, int *flag_IO){
    StartTime = st;
    strcpy(filename, s);
    if(strlen(filename) < 1){
        fprintf(stderr, "WARNING: filename (%s)is empty.\n;", filename);
    }
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    
    Interval = dt;
    NumVar = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            NumVar++;
        }
    }
    if(NumVar <= 0){
        fprintf(stderr, "WARNING: Empty columns in %s.\n;", filename);
    }
    buffer = new double[NumVar];
    PrintVar = new double*[NumVar];
    icol    = new double[NumVar];
    int k = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            PrintVar[k] = &x[i];
            icol[k] = (double) (i + 1);
            k++;
        }
    }
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1.;
    }
}

void Print_Ctrl::InitIJ(long st, int n, const char *s, int dt, double **x, int j, int iFlux, int *flag_IO){
    StartTime = st;
    NumVar = n;
    PrintVar = new double*[NumVar];
    buffer = new double[NumVar];
    icol    = new double[NumVar];
    strcpy(filename, s);
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    Interval = dt;
    NumVar = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            NumVar++;
        }
    }
    if(NumVar <= 0){
        fprintf(stderr, "WARNING: Empty columns in %s.\n;", filename);
    }
    buffer = new double[NumVar];
    PrintVar = new double*[NumVar];
    int k = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            PrintVar[k] = &(x[i][j]);
            icol[k] = (double) (i + 1);
            k++;
        }
    }
    
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1;
    }
}
Print_Ctrl::~Print_Ctrl(){
    if(NumVar > 0){
        if(PrintVar != NULL ) delete[] PrintVar;
        if(buffer != NULL ) delete[] buffer;
    }
    close_file();
}
void Print_Ctrl::fun_printBINARY(double t, double dt){
    fwrite (&t, sizeof (double), 1, fid_bin);
    fwrite (buffer, sizeof (double), NumVar, fid_bin);
    if(global_fflush_mode){
        fflush (fid_bin); // DEBUG
    }
}
void Print_Ctrl::fun_printASCII(double t, double dt){
    fprintf(fid_asc, "%.1f\t", t);
    for (int i = 0; i < NumVar; i++){
        fprintf (fid_asc, "%e\t", buffer[i]);
    }
    fprintf (fid_asc, "\n");
    if(global_fflush_mode){
//    fflush (fid_asc); // DEBUG
    }
}
void Print_Ctrl::close_file(){
    if (Binary){
        if(fid_bin != NULL)
            fclose(fid_bin);
    }
    if (Ascii){
        if(fid_asc != NULL)
            fclose(fid_asc);
    }
    
}
void Print_Ctrl::PrintData(double dt, double t){
    long    t_long;
    NumUpdate++; /* Number of times to push data into the buffer*/
    for (int i = 0; i < NumVar; i++){
        buffer[i] += *(PrintVar[i]);
    }
    t_long = (long int)t;
    if ((t_long % Interval) == 0){
        for (int i = 0; i < NumVar; i++){
            buffer[i] *= tau / NumUpdate ; /* Get the mean in the time-interval*/
        }
        NumUpdate = 0;
        if(Ascii){
            fun_printASCII(t-Interval, dt);
        }
        if(Binary){
            fun_printBINARY(t-Interval, dt);
        }
        for (int i = 0; i < NumVar; i++){
            buffer[i] = 0.;  /* Reset the buffer */
        }
    }
}
