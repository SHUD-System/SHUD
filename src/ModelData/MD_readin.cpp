#include "Model_Data.hpp"

void Model_Data::read_calib(const char *fn){
    gc.read(fn);
}
void Model_Data::read_para(const char *fn){
    CS.read(fn);
}

void Model_Data::setIO_ele(int x){
    for(int i = 0; i < NumEle; i++){
        io_ele[i] = x;
    }
}
void Model_Data::setIO_riv(int x){
    for(int i = 0; i < NumRiv; i++){
        io_riv[i] = x;
    }
}
void Model_Data::setIO_lake(int x){
    for(int i = 0; i < NumLake; i++){
        io_lake[i] = x;
    }
}
void Model_Data::read_cfgout(const char *fn){
    TabularData tb;
    int idx, val;
    
    io_ele      = new int[NumEle];
    if(NumRiv > 0) io_riv      = new int[NumRiv];
    if(NumLake > 0) io_lake     = new int[NumLake];
    setIO_ele(1);
    setIO_riv(1);
    setIO_lake(1);
    
    FILE *fp;
    fp =  fopen(fn, "r");
    if (fp == NULL) {
        /* IF file is missing.
         ALL is 1*/
    }else{
        int nline;
        nline = tb.read(fp);
        if(tb.ncol != 2){
            printf("The column of element in cfg.output should be: \n");
            printf("%s\t%s\n", "index", "OFF/ON");
            printf("Actual column is:\n%s", tb.header);
            myexit(ERRFileIO);
        }
        setIO_ele(atoi(tb.header));
        for (int i = 0; i < nline; i++){
            idx     = (int) tb.x[i][0] - 1;
            val     = (int) tb.x[i][1];
            if(val > 0){
                io_ele[idx] = 1;
            }else{
                io_ele[idx] = 0;  /* output OFF for the element */
            }
            printf("%d: %d, %d\n", i+1, idx, io_ele[idx]);
        }
        
        if(NumRiv > 0 ){
            nline = tb.read(fp);
            if(tb.ncol != 2){
                printf("The column of river in cfg.output should be: \n");
                printf("%s\t%s\n", "index", "OFF/ON");
                printf("Actual column is:\n%s", tb.header);
                myexit(ERRFileIO);
            }
            setIO_riv(atoi(tb.header));
            for (int i = 0; i < nline; i++){
                idx     = (int) tb.x[i][0] - 1;
                val     = (int) tb.x[i][1];
                if(val > 0){
                    io_riv[idx] = 1;
                }else{
                    io_riv[idx] = 0;  /* output OFF for the River reach */
                }
                printf("%d: %d, %d\n", i+1, idx, io_ele[idx]);
            }
        }
            
        if(NumLake > 0 ){
            nline = tb.read(fp);
            if(tb.ncol != 2){
                printf("The column of lake in cfg.output should be: \n");
                printf("%s\t%s\n", "index", "OFF/ON");
                printf("Actual column is:\n%s", tb.header);
                myexit(ERRFileIO);
            }
            setIO_lake(atoi(tb.header));
            for (int i = 0; i < nline; i++){
                idx     = (int) tb.x[i][0] - 1;
                val     = (int) tb.x[i][1];
                if(val > 0){
                    io_lake[idx] = 1;
                }else{
                    io_lake[idx] = 0;  /* output OFF for the River reach */
                }
                printf("%d: %d, %d\n", i+1, idx, io_ele[idx]);
            }
        }
        fclose(fp);
    }
}
void Model_Data::read_rivseg(const char *fn){
    TabularData tb;
    NumSegmt = tb.read(fn);
    RivSeg = new RiverSegement[NumSegmt];
    printf("\tNumber of River segments: %d\n", NumSegmt);
    if(tb.ncol != 4){
        printf("The column of River Parameter should be: \n");
        printf("%s\t%s\t%s\t%s\n", "index", "iRiv", "iEle", "Length");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    for (int i=0; i<NumSegmt ; i++){
        RivSeg[i].index    = (int) tb.x[i][0];
        RivSeg[i].iRiv     = (int) tb.x[i][1];
        RivSeg[i].iEle     = (int) tb.x[i][2];
        RivSeg[i].length   = (double) tb.x[i][3];
    }
}
void Model_Data::read_riv(const char *fn){
    int iflag = 0;
    TabularData tb;
    FILE *fp;
    fp =  fopen(fn, "r");
    CheckFile(fp, fn);
    
    NumRiv = tb.read(fp);
    Riv = new _River[NumRiv];
    
    if(tb.ncol != 6){
        printf("The column of RiverAtt should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\n", "index", "Down", "type", "Bedslope", "Length", "BC");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    for (int i = 0; i < NumRiv; i++){
        Riv[i].index    = (int) tb.x[i][0];
        Riv[i].down     = (int) tb.x[i][1];
        Riv[i].type     = (int) tb.x[i][2];
        Riv[i].BedSlope = (double) tb.x[i][3];
        Riv[i].Length   = (double) tb.x[i][4];
        Riv[i].BC   = (double) tb.x[i][5];
        if(Riv[i].Length < ZERO){
            printf("Debug: %d\t%d\t%d\t%d\t%f\t%f\n", i, Riv[i].index, Riv[i].down, Riv[i].type, Riv[i].BedSlope, Riv[i].Length );
            printf("\tThe RIV %d is ZERO/NEGATIVE length\n", i + 1 );
            iflag = 1;
        }
        if(Riv[i].down < 1){
            printf("\tThe downstream of RIV %d is negtive (Outlet, %d)\n", i + 1, Riv[i].down);
        }
        if(Riv[i].BC > 0){
            printf("BC[%d] (Neumann Condition) applied to River %d\n", Riv[i].BC, i+1);
            irBC1 = 1;
        }
        if(Riv[i].BC < 0){
            printf("BC[%d] (Dirichlet Condition) applied to River %d\n", -Riv[i].BC, i+1);
            irBC2 = 1;
        }
//        printf("debug%d\t%d\t%d\t%d\t%f\t%f\n", i, Riv[i].index, Riv[i].down, Riv[i].type, Riv[i].BedSlope, Riv[i].Length );
    }
    if(iflag){
        myexit(ERRDATAIN);
    }
    NumRivType = tb.read(fp);
    Riv_Type = new river_para[NumRivType];
    
    if(tb.ncol != 9){
        printf("The column of River Parameter should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "index", "depth", "bankslope", "BottomWidth", "Sinuosity", "Rough", "Cwr", "KsatH", "BedThick");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    for (int i = 0; i < NumRivType; i++){
        Riv_Type[i].InitValue( tb.x[i] );
    }
//    
//    NumRivNode = tb.read(fp);
//    rivNode = new _Node[NumRivNode];
//    for (int i = 0; i < NumRivNode; i++){
//        rivNode[i].index    = (int) tb.x[i][0];
//        rivNode[i].x        = (double) tb.x[i][1];
//        rivNode[i].y        = (double) tb.x[i][2];
//        rivNode[i].zmax     = (double) tb.x[i][3];
//    }
    fclose(fp);
}

void Model_Data::read_mesh(const char *fn){
    TabularData tb;
    int ncol;
    /*========== open *.mesh file ==========*/
    FILE * fp =  fopen(fn, "r");
    CheckFile(fp, fn);
    
    NumEle = tb.read(fp, &ncol);
    if(tb.ncol != 8){
        printf("The column of Mesh in .mesh should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               "index", "Node1", "Node2", "Node3",
               "nabr1", "nabr2", "nabr3");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    Ele = new _Element[NumEle];
    for (int i=0; i<NumEle ; i++){
        Ele[i].index = (int) tb.x[i][0];
        Ele[i].node[0] = (int) tb.x[i][1];
        Ele[i].node[1] = (int) tb.x[i][2];
        Ele[i].node[2] = (int) tb.x[i][3];
        Ele[i].nabr[0] = (int) tb.x[i][4];
        Ele[i].nabr[1] = (int) tb.x[i][5];
        Ele[i].nabr[2] = (int) tb.x[i][6];
    }
    /* read in nodes information */
    NumNode = tb.read(fp);
    if(tb.ncol != 5){
        printf("The column of Points in .mesh should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t\n", "index", "x", "y", "AqD", "zmax");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    Node = new _Node[NumNode];
    for (int i=0; i<NumNode ; i++){
        Node[i].index   = (int) tb.x[i][0];
        Node[i].x       = tb.x[i][1];
        Node[i].y       = tb.x[i][2];
        Node[i].AqD     = tb.x[i][3];
        Node[i].zmax    = tb.x[i][4];
//        Node[i].Init();
    }
    fclose(fp);
}
void Model_Data::read_att(const char *fn){
    int nr;
    TabularData tb;
    nr = tb.read(fn);
    if( nr != NumEle){
        fprintf(stderr, "\nWARNING: number of rows (%d) in att DOSE NOT match the number of cell in .mesh file (%d).\n Press anykey to continue ...\n", nr, NumEle);
//        getchar();
    }
    if(tb.ncol != 9){
        printf("The column of .att should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               "index", "iSoil", "iGeol", "iLC", "iForc", "iMF", "iBC", "iSS", "Lake");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    for (int i=0; i<NumEle; i++){
        // index are ignored, since the read_mesh already take care of it.
        Ele[i].iSoil    = (int) tb.x[i][1];
        Ele[i].iGeol    = (int) tb.x[i][2];
        Ele[i].iLC      = (int) tb.x[i][3];
        Ele[i].iForc    = (int) tb.x[i][4];
        Ele[i].iMF      = (int) tb.x[i][5];
        Ele[i].iBC      = (int) tb.x[i][6];
        Ele[i].iSS      = (int) tb.x[i][7];
        Ele[i].iLake    = (int) tb.x[i][8];
        if( lakeon <=0 && Ele[i].iLake >0 ){
            lakeon = 1;
        }
//        printf("debug%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",i+1, Ele[i].iSoil, Ele[i].iGeol, Ele[i].iLC, Ele[i].iForc, Ele[i].iMF, Ele[i].iBC, Ele[i].iSS, Ele[i].iLake);
        if(Ele[i].iBC > 0){
            ieBC1 = 1;
        }
        if(Ele[i].iBC < 0){
            ieBC2 = 1;
        }
        if(Ele[i].iSS != 0){
            ieSS = 1;
        }
    }
}

void Model_Data::read_soil(const char *fn){
    TabularData tb;
    NumSoil = tb.read(fn);
    if(tb.ncol != 9){
        printf("The column of Soil file should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n", "index", "infKsatV", "ThetaS", "ThetaR", "infD", "Alpha","Beta","hAreaF", "macKsatV");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    Soil = new Soil_Layer[NumSoil];
    for (int i=0; i<NumSoil; i++){
        Soil[i].index       = (int) tb.x[i][0];
        Soil[i].infKsatV    = (double) tb.x[i][1] / 1440.;
        Soil[i].ThetaS      = (double) tb.x[i][2];
        Soil[i].ThetaR      = (double) tb.x[i][3];
        Soil[i].infD        = (double) tb.x[i][4];
        Soil[i].Alpha       = (double) tb.x[i][5];
        Soil[i].Beta        = (double) tb.x[i][6];
        Soil[i].hAreaF      = (double) tb.x[i][7];
        Soil[i].macKsatV    = (double) tb.x[i][8] / 1440.;
        //        printf("dbg%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
        //               (Soil[i].index),
        //               (Soil[i].infKsatV),
        //               (Soil[i].ThetaS),
        //               (Soil[i].ThetaR),
        //               (Soil[i].infD),
        //               (Soil[i].Alpha),
        //               (Soil[i].Beta),
        //               (Soil[i].hAreaF),
        //               (Soil[i].macKsatV)
        //               );
        
    }
}
void Model_Data::read_geol(const char *fn){
    TabularData tb;
    NumGeol = tb.read(fn);
    if(tb.ncol != 8){
        printf("The column of Geol file should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "index", "KsatH","KsatV", "geo_ThetaS", "geo_ThetaR", "geo_ThetaR", "macKsatH","macD");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    Geol = new Geol_Layer[NumGeol];
    for (int i=0; i<NumGeol; i++){
        Geol[i].index       = (int) tb.x[i][0];
        Geol[i].KsatH       = (double) tb.x[i][1] / 1440.;
        Geol[i].KsatV       = (double) tb.x[i][2] / 1440.;
        Geol[i].geo_ThetaS  = (double) tb.x[i][3];
        Geol[i].geo_ThetaR  = (double) tb.x[i][4];
        Geol[i].geo_vAreaF  = (double) tb.x[i][5];
        Geol[i].macKsatH    = (double) tb.x[i][6] / 1440.;
        Geol[i].macD        = (double) tb.x[i][7];
    }
}

void Model_Data::read_lc(const char *fn){
    TabularData tb;
    NumLC = tb.read(fn);
    if(tb.ncol < 7 | tb.ncol > 8){
        printf("The column of LandCover file should be: \n");
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "index",
               "Albedo", "VegFrac", "Rough","RzD", "SoilDgrd", "ImpAF");
        printf("Actual column is:\n%s", tb.header);
        myexit(ERRFileIO);
    }
    LandC = new Landcover[NumLC];
    for (int i=0; i<NumLC ; i++){
        LandC[i].index  = (int) tb.x[i][0];
        LandC[i].Albedo = (double) tb.x[i][1];
        LandC[i].VegFrac = (double) tb.x[i][2];
        LandC[i].Rough  = (double) tb.x[i][3] / 60.; /* [s m^(-1/3)]  to [min m^(-1/3)]*/
        LandC[i].RzD    = (double) tb.x[i][4];
        LandC[i].SoilDgrd = (double) tb.x[i][5];
        LandC[i].ImpAF  = (double) tb.x[i][6];
//        printf( "%d %lf %lf %lf %lf %lf %lf\n",
//               (LandC[i].index),
//               (LandC[i].Albedo),
//               (LandC[i].VegFrac),
//               (LandC[i].Rough),
//               (LandC[i].RzD),
//               (LandC[i].SoilDgrd),
//               (LandC[i].ImpAF)
//               );
    }
}
void Model_Data::read_forc_csv(const char *fn){
    FILE *fp = fopen(fn, "r");
    CheckFile(fp, fn);
    char path[MAXLEN];
    char str[MAXLEN];
    char shortname[MAXLEN];
    char longname[MAXLEN];
    int id;
    double lon, lat;
    path[0] = '\0';
    str[0] = '\0';
    fgets(str, MAXLEN, fp); // Dimension of the table
    sscanf(str, "%d %ld", &NumForc, &ForcStartTime);
    tsd_weather = new _TimeSeriesData[NumForc];
    for(int i=0; i < NumForc; i++){
        tsd_weather[i].initialize(Nforc + 1); /* Nforc= number of forcing variables. */
    }
    fgets(str, MAXLEN, fp);
    if( strlen(str) > 1){
        sscanf(str, "%s", path);
    }
    fgets(str, MAXLEN, fp);
#ifdef DEBUG
    printf("%s", str);
#endif
    for(int i=0; i < NumForc && fgets(str, MAXLEN, fp); i++)
    {
        sscanf(str, "%d %lf %lf %lf %lf %lf %s", &id, &lon, &lat, tsd_weather[i].xyz, tsd_weather[i].xyz+1, tsd_weather[i].xyz+2, shortname);
        sprintf(longname, "%s/%s", path, shortname);
#ifdef DEBUG
        printf("debug:%d %d, %lf, %lf, %lf, %lf, %lf, %s\n", i+1, id, lon, lat,
               tsd_weather[i].xyz[0], tsd_weather[i].xyz[1], tsd_weather[i].xyz[2], longname);
#endif
        tsd_weather[i].fn.assign(longname);
    }
    fclose(fp);
    printf("\tNumber of csv files: %d\n", NumForc);
    for(int i=0; i < NumForc; i++){
        printf("\t Reading %d/%d: \t%s\n", i+1, NumForc, tsd_weather[i].fn.c_str());
        tsd_weather[i].read_csv();
    }
}
void Model_Data::loadinput(){
    int nt=1;
    int flag = 1;
    if(flag) fprintf(stdout,"* \t Project name: %s\n", pf_in->projectname);
    if(flag) fprintf(stdout,"* \t Project input folder: %s\n", pf_in->inpath);
    if(flag) fprintf(stdout,"* \t Project output folder: %s\n", pf_in->outpath);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_para);
    read_para(pf_in->file_para);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_riv);
    read_riv(pf_in->file_riv);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_rivseg);
    read_rivseg(pf_in->file_rivseg);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_mesh);
    read_mesh(pf_in->file_mesh);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_att);
    read_att(pf_in->file_att);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_soil);
    read_soil(pf_in->file_soil);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_geol);
    read_geol(pf_in->file_geol);
    //    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_ibc);
    //    read_ibc(pf_in->file_ibc);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_lc);
    read_lc(pf_in->file_lc);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_forc);
    read_forc_csv(pf_in->file_forc);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_lai);
    read_lai(pf_in->file_lai);
//    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_rl);
//    read_rl(pf_in->file_rl);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_mf);
    read_mf(pf_in->file_mf);
    
    //    read_forc_binary(pf_in->file_forc);
    if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_calib);
    read_calib(pf_in->file_calib);
    if(ieBC1){
        if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_ebc1);
        read_bcEle1(pf_in->file_ebc1);
    }
    if(ieBC2){
        if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_ebc2);
        read_bcEle2(pf_in->file_ebc2);
    }
    if(irBC1){
        if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_rbc1);
        read_bcRiv1(pf_in->file_rbc1);
    }
    if(irBC2){
        if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_rbc2);
        read_bcRiv2(pf_in->file_rbc2);
    }
    if(ilBC1){
        if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_rbc2);
        read_bcLake1(pf_in->file_lbc1);
    }
    if(ilBC2){
        if(flag) fprintf(stdout,"%d \t Reading file: %s\n", nt++,pf_in->file_rbc2);
        read_bcLake2(pf_in->file_lbc2);
    }
    this->CS.num_threads = max(CS.num_threads,  pf_in->numthreads); /* Number of threads for OpenMP */
}
//void Model_Data::read_rl(const char *fn){
//    tsd_RL.fn = fn;
//    tsd_RL.readDimensions();
//    tsd_RL.read_csv();
//}
void Model_Data::read_lai(const char *fn){
    tsd_LAI.fn = fn;
    tsd_LAI.readDimensions();
    tsd_LAI.read_csv();
}
void Model_Data::read_mf(const char *fn){
    tsd_MF.fn = fn;
    tsd_MF.readDimensions();
    NumMeltF = tsd_MF.get_Ncol();
    tsd_MF.read_csv();
}
void Model_Data::read_ssEle(const char *fn){
    tsd_eleSS.fn = fn;
    tsd_eleSS.readDimensions();
    NumSSEle = tsd_eleSS.get_Ncol();
    tsd_eleSS.read_csv();
}
void Model_Data::read_bcEle1(const char *fn){
    tsd_eyBC.fn = fn;
    tsd_eyBC.readDimensions();
    NumBCEle1 = tsd_eyBC.get_Ncol();
    tsd_eyBC.read_csv();
}
void Model_Data::read_bcEle2(const char *fn){
    tsd_eqBC.fn = fn;
    tsd_eqBC.readDimensions();
    NumBCEle2 = tsd_eqBC.get_Ncol();
    tsd_eqBC.read_csv();
}
void Model_Data::read_bcRiv1(const char *fn){
    tsd_ryBC.fn = fn;
    tsd_ryBC.readDimensions();
    NumBCRiv1 = tsd_ryBC.get_Ncol();
    tsd_ryBC.read_csv();
}
void Model_Data::read_bcRiv2(const char *fn){
    tsd_rqBC.fn = fn;
    tsd_rqBC.readDimensions();
    NumBCRiv2 = tsd_rqBC.get_Ncol();
    tsd_rqBC.read_csv();
}
void Model_Data::read_bcLake1(const char *fn){
    tsd_lyBC.fn = fn;
    tsd_lyBC.readDimensions();
    NumBCLake1 = tsd_lqBC.get_Ncol();
    tsd_lyBC.read_csv();
}
void Model_Data::read_bcLake2(const char *fn){
    tsd_lqBC.fn = fn;
    tsd_lqBC.readDimensions();
    NumBCLake2 = tsd_lqBC.get_Ncol();
    tsd_lqBC.read_csv();
}
void Model_Data::FreeData(){
    
    for (int i = 0; i < NumEle; i++) {
        delete[] QeleSurf[i] ;
        delete[] QeleSub[i] ;
    }
    
    delete[]    io_ele;
    delete[]    io_riv;
    delete[]    io_lake;
    
    delete[]    QeleSurf;
    delete[]    QeleSub;
    delete[]    QeleSurfTot;
    delete[]    QeleSubTot;
    delete[]    QoutSurf; //5
    
    delete[]    Qe2r_Surf; // 5.1
    delete[]    Qe2r_Sub; // 5.2
    
    delete[]    qElePrep;
    delete[]    qEleTF;
    delete[]    qEleETP;
    delete[]    qEleETA;
    delete[]    qEleE_IC;
    delete[]    qEleTrans;
    delete[]    qEleEvapo;
    delete[]    qPotEvap;
    delete[]    qPotTran;
    
    delete[]    qEs;
    delete[]    qEu;
    delete[]    qEg;
    delete[]    qTu;
    delete[]    qTg;
    
    delete[]    qEleETloss; //10
    delete[]    iBeta; // 10.1
    
    delete[]    qEleNetPrep;
    delete[]    qEleInfil;
    delete[]    qEleExfil;
    delete[]    qEleRecharge; //13
    
    delete[]    yEleIS;
    delete[]    yEleISmax;
    delete[]    yEleISsnowmax; //16
    
    delete[]    yEleSnow;
    delete[]    yEleSnowGrnd;
    delete[]    yEleSnowCanopy;
    delete[]    yEleGW; //20
    
    delete[]    yEleSurf;
    delete[]    yEleUnsat;
    delete[]    yEleWetFront; //23
    
    delete[]    yRivStg;
    delete[]    QrivSurf;
    delete[]    QrivSub; // 26
    
    delete[]    QrivDown;
    delete[]    QrivUp;  // 28
    delete[]    QsegSurf;
    delete[]    QsegSub;
    
    
    delete[]    fu_Surf;
    delete[]    fu_Sub;
//    delete[]    AccT;
    
    if(NumLake > 0){
        delete[]    yLakeStg;
        delete[]    QLakeSurf; //30
        
        delete[]    QLakeSub;
        delete[]    QLakeRivIn;
        delete[]    QLakeRivOut;
        delete[]    qLakePrcp;
        delete[]    qLakeEvap; //34
    }
    delete[]    uYsf; //35.
    delete[]    uYus; //35.
    delete[]    uYgw; //35.
    delete[]    uYriv; //35.1
    delete[]    uYlake; //35.1
    delete[]    globalY;
    
    delete[]    t_prcp;
    delete[]    t_temp;
    delete[]    t_rh;
    delete[]    t_wind;
    delete[]    t_rn;
    delete[]    t_vp;
    delete[]    t_lai;
    delete[]    t_mf;
//    delete[]    t_hc;
    
    //read_rivseg()
    delete[] RivSeg;
    
    /* free river, read_riv */
    delete[] Riv;
    delete[] Riv_Type;
    delete[] rivNode;
    
    /* free mesh, read_mesh */
    delete[] Ele;
    delete[] Node;
    
    /* free soil, read_soil */
    delete[] Soil;
    /* free geol, read_geol */
    delete[] Geol;
    /* free lc, read_lc */
    delete[] LandC;
    
    /* free forcing data */
    delete[] tsd_weather;
    
    /* MD::initialize() */
    delete flood;
}

