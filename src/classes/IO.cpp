#include "IO.hpp"

void FileIn::saveProject(){
    FILE *fp;
    char fn[MAXLEN];
    sprintf(fn, "%s/%s.SHUD", outpath, projectname);
    fp = fopen(fn,"w");
    CheckFile(fp, fn);
    /*  */
    /*  */
    /* Configure */
    fprintf(fp, "PRJ \t %s\n", projectname);
    fprintf(fp, "INPATH \t %s\n", inpath);
    fprintf(fp, "OUTPATH \t %s\n", outpath);
    
    /* spatial data */
    fprintf(fp, "MESH \t %s\n", file_mesh);
    fprintf(fp, "ATT \t %s\n", file_att);
    fprintf(fp, "LAKE \t %s\n", file_lake);
    fprintf(fp, "RIV \t %s\n", file_riv);
    fprintf(fp, "RIVSEG \t %s\n", file_rivseg);
    /* calibration & config */
    fprintf(fp, "CALIB \t %s\n", file_calib);
    fprintf(fp, "PARA \t %s\n", file_para);
    fprintf(fp, "INIT \t %s\n", file_init);
    /* Parameters */
    fprintf(fp, "LC \t %s\n", file_lc);
    fprintf(fp, "SOIL \t %s\n", file_soil);
    fprintf(fp, "GEOL \t %s\n", file_geol);
    /* TSD */
    fprintf(fp, "FORC \t %s\n", file_forc);
    fprintf(fp, "LAI \t %s\n", file_lai);
    fprintf(fp, "MF \t %s\n", file_mf);
    fprintf(fp, "RL \t %s\n", file_rl);
    fprintf(fp, "LCM \t %s\n", file_lcm);
    /* bc */
    fprintf(fp, "EleBC1 \t %s\n", file_ebc1);
    fprintf(fp, "EleBC2 \t %s\n", file_ebc2);
    fprintf(fp, "RivBC1 \t %s\n", file_rbc1);
    fprintf(fp, "RivBC2 \t %s\n", file_rbc2);
    fprintf(fp, "LakeBC1 \t %s\n", file_lbc1);
    fprintf(fp, "LakeBC2 \t %s\n", file_lbc2);
    
    fclose(fp);
}

void FileIn:: setInFilePath(char * indir, char *  pjrname, int nthred){
    numthreads = nthred;
    setInFilePath(indir, pjrname);
}
void FileIn:: setInFilePath(char * indir, char *  pjrname){
    /*Model io*/
    sprintf(inpath, "%s", indir);
    sprintf(outpath, "output/%s.out", pjrname);
    sprintf(logfile, "%s.log", outpath);
    sprintf(projectname, "%s",  pjrname);
    
    /*Spatial Data*/
    sprintf(file_mesh, "%s/%s.%s", inpath, projectname, "sp.mesh");
    sprintf(file_att, "%s/%s.%s", inpath, projectname, "sp.att");
    sprintf(file_riv, "%s/%s.%s", inpath, projectname, "sp.riv");
    sprintf(file_rivseg, "%s/%s.%s", inpath, projectname, "sp.rivseg");
    sprintf(file_lake, "%s/%s.%s", inpath, projectname, "sp.lake");
    sprintf(file_lake_bathy, "%s/%s.%s", inpath, projectname, "lake.bathy");
    sprintf(file_lake_ic, "%s/%s.%s", inpath, projectname, "lake.ic");
    
    /* physical parameters */
    sprintf(file_lc, "%s/%s.%s", inpath, projectname, "para.lc");
    sprintf(file_soil, "%s/%s.%s", inpath, projectname, "para.soil");
    sprintf(file_geol, "%s/%s.%s", inpath, projectname, "para.geol");
    
    /* model configuration */
    sprintf(file_cfgout, "%s/%s.%s", inpath, projectname, "cfg.output");
    sprintf(file_para, "%s/%s.%s", inpath, projectname, "cfg.para");
    sprintf(file_calib, "%s/%s.%s", inpath, projectname, "cfg.calib");
    sprintf(file_init, "%s/%s.%s", inpath, projectname, "cfg.ic");
    
    /* Time-series data */
    sprintf(file_forc, "%s/%s.%s", inpath, projectname, "tsd.forc");
    sprintf(file_lai, "%s/%s.%s", inpath, projectname, "tsd.lai");
    sprintf(file_mf, "%s/%s.%s", inpath, projectname, "tsd.mf");
    sprintf(file_rl, "%s/%s.%s", inpath, projectname, "tsd.rl");
    sprintf(file_lcm, "%s/%s.%s", inpath, projectname, "tsd.lcm");
    sprintf(file_obs, "%s/%s.%s", inpath, projectname, "tsd.obs");
    
    sprintf(file_ebc1, "%s/%s.%s", inpath, projectname, "tsd.ebc1");
    sprintf(file_ebc2, "%s/%s.%s", inpath, projectname, "tsd.ebc2");
    sprintf(file_rbc1, "%s/%s.%s", inpath, projectname, "tsd.rbc1");
    sprintf(file_rbc2, "%s/%s.%s", inpath, projectname, "tsd.rbc2");
    sprintf(file_lbc1, "%s/%s.%s", inpath, projectname, "tsd.lbc1");
    sprintf(file_lbc2, "%s/%s.%s", inpath, projectname, "tsd.lbc2");
}

FileOut::FileOut(){
    setsuffix("");
}
void FileOut::setsuffix(const char *s){
    strcpy(suffix, s);
}
void FileOut::setOutFilePath(char *outdir, char *  prjname){
    sprintf(outpath, "%s", outdir);
    sprintf(projectname, "%s", prjname);
    updateFilePath();
}
void FileOut::updateFilePath(){
    createDir();
    //    printf("\nOutput path is \"%s\"\n", outpath);
    /****************************************************
     format： Prjname+Surffix.Component+Type+Value.Extension
     i.g. pj0.eleysurf.dat/ prject001.rivqdown.dat
     pj - project name. len = 1~n
     0  - suffix. len = 1~n
     ele- componencts, len = 3
            riv = River
            ele = Element
            lak = Lake
     y  - Type,
            y - State Variable[m],
            v- specific flux [m3/m2/day],
            q - flux variable[m3/day]; len =1
     surf- variable, len = 1~n;
     dat- File extension, len =3
            dat = bindary,
            csv = ASCII spreadsheet.
     *****************************************************/
    //rivers
    sprintf(riv_Q_down, "%s/%s%s.rivqdown", outpath, projectname, suffix);
    sprintf(riv_Q_up, "%s/%s%s.rivqup", outpath, projectname, suffix);
    sprintf(riv_Q_surf, "%s/%s%s.rivqsurf", outpath, projectname, suffix);
    sprintf(riv_Q_sub, "%s/%s%s.rivqsub", outpath, projectname, suffix);
    sprintf(riv_y_stage, "%s/%s%s.rivystage", outpath, projectname, suffix);
    //cells
    sprintf(ele_y_snow, "%s/%s%s.eleysnow", outpath, projectname, suffix);
    sprintf(ele_y_ic, "%s/%s%s.eleyic", outpath, projectname, suffix);
    sprintf(ele_y_surf, "%s/%s%s.eleysurf", outpath, projectname, suffix);
    sprintf(ele_y_unsat, "%s/%s%s.eleyunsat", outpath, projectname, suffix);
    sprintf(ele_y_gw, "%s/%s%s.eleygw", outpath, projectname, suffix);
    //cell-fluxes
    sprintf(ele_q_ET[0], "%s/%s%s.elevetic", outpath, projectname, suffix);
    sprintf(ele_q_ET[1], "%s/%s%s.elevettr", outpath, projectname, suffix);
    sprintf(ele_q_ET[2], "%s/%s%s.elevetev", outpath, projectname, suffix);
    sprintf(ele_q_ETP, "%s/%s%s.elevetp", outpath, projectname, suffix);
    sprintf(ele_q_ETA, "%s/%s%s.eleveta", outpath, projectname, suffix);
    sprintf(ele_q_prcp, "%s/%s%s.elevprcp", outpath, projectname, suffix);
    sprintf(ele_q_netprcp, "%s/%s%s.elevnetprcp", outpath, projectname, suffix);
    sprintf(ele_q_infil, "%s/%s%s.elevinfil", outpath, projectname, suffix);
    sprintf(ele_q_exfil, "%s/%s%s.elevexfil", outpath, projectname, suffix);
    sprintf(ele_q_rech, "%s/%s%s.elevrech", outpath, projectname, suffix);
    
    //    sprintf(ele_Q_surf, "%s/%s%s.eleqsurf", outpath, projname);
    //    sprintf(ele_Q_sub, "%s/%s%s.eleqsub", outpath, projname);
    sprintf(ele_Q_subTot, "%s/%s%s.eleqsub", outpath, projectname, suffix);
    sprintf(ele_Q_sub0, "%s/%s%s.eleqsub1", outpath, projectname, suffix);
    sprintf(ele_Q_sub1, "%s/%s%s.eleqsub2", outpath, projectname, suffix);
    sprintf(ele_Q_sub2, "%s/%s%s.eleqsub3", outpath, projectname, suffix);
    sprintf(ele_Q_surfTot, "%s/%s%s.eleqsurf", outpath, projectname, suffix);
    sprintf(ele_Q_surf0, "%s/%s%s.eleqsurf1", outpath, projectname, suffix);
    sprintf(ele_Q_surf1, "%s/%s%s.eleqsurf2", outpath, projectname, suffix);
    sprintf(ele_Q_surf2, "%s/%s%s.eleqsurf3", outpath, projectname, suffix);
    sprintf(ele_Q_rsurf, "%s/%s%s.eleqrsurf", outpath, projectname, suffix);
    sprintf(ele_Q_rsub, "%s/%s%s.eleqrsub", outpath, projectname, suffix);
    
    //cell_wb
    sprintf(ewb_q_in, "%s/%s%s.ewbqin", outpath, projectname, suffix);
    sprintf(ewb_q_out, "%s/%s%s.ewbqout", outpath, projectname, suffix);
    sprintf(ewb_q_io, "%s/%s%s.ewbqio", outpath, projectname, suffix);
    sprintf(ewb_dh, "%s/%s%s.ewbydh", outpath, projectname, suffix);
    
    sprintf(lake_Q_rivin, "%s/%s%s.lakqrivin", outpath, projectname, suffix);
    sprintf(lake_Q_rivout, "%s/%s%s.lakqrivout", outpath, projectname, suffix);
    sprintf(lake_Q_surf, "%s/%s%s.lakqsurf", outpath, projectname, suffix);
    sprintf(lake_Q_sub, "%s/%s%s.lakqsub", outpath, projectname, suffix);
    sprintf(lake_y_stage, "%s/%s%s.lakystage", outpath, projectname, suffix);
    sprintf(lake_a_area, "%s/%s%s.lakatop", outpath, projectname, suffix);
    sprintf(lake_q_evap, "%s/%s%s.lakvevap", outpath, projectname, suffix);
    sprintf(lake_q_prcp, "%s/%s%s.lakvprcp", outpath, projectname, suffix);
    
    sprintf(Init_update, "%s/%s%s.cfg.ic.update", outpath, projectname, suffix);
    sprintf(Init_bak, "%s/%s%s.cfg.ic.bak", outpath, projectname, suffix);
    
    sprintf(Calib_bak, "%s/%s%s.cfg.calib.bak", outpath, projectname, suffix);
    
    sprintf(floodout, "%s/%s%s.flood.csv", outpath, projectname, suffix);
    sprintf(obs_sim, "%s/%s%s.ovs.csv", outpath, projectname, suffix);
    
    sprintf(File_Time, "%s/%s%s.time.csv", outpath, projectname, suffix);
    fid_time = fopen(File_Time, "w");
    CheckFile(fid_time, "Time log file");
    fprintf(fid_time, "time_Minutes \t Time_Days \t Task_perc \t CPUTime_s \t WallTime_s \t Num_fcall \n");

}
void FileOut::createDir(){
    mkdir_p(outpath, 0777);
}
void FileIn::setOutpath(const char *fn){
    strcpy(outpath, fn);
}
void FileIn::setCalibFile(const char *fn){
    strcpy(file_calib, fn);
}
void FileIn::readProject(const char *fn){
    
    char    str[MAXLEN];
    char    optstr[MAXLEN];
    char    val[MAXLEN];
    FILE *fp;
    fp = fopen(fn, "r");
    CheckFile(fp, fn);
    /* Read through parameter file to find parameters */

    while (fgets(str, MAXLEN, fp)) {
        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0' || str[0] == ' ')
        {
            continue;
        }
        sscanf (str, "%s %s\n", optstr, val);
        /* Get Model Parameters */
        if (strcasecmp ("PRJ", optstr) == 0){
            strcpy(projectname, val);
            sprintf(inpath, "input/%s", projectname);
            sprintf(outpath, "output/%s.out", projectname);
            setInFilePath(inpath, projectname);
        }
        else if (strcasecmp ("INPATH", optstr) == 0)
            strcpy(inpath, val);
        else if (strcasecmp ("OUTPATH", optstr) == 0)
            strcpy(outpath, val);
        /* SPATIAL DATA */
        else if (strcasecmp ("MESH", optstr) == 0)
            strcpy(file_mesh, val);
        else if (strcasecmp ("ATT", optstr) == 0)
            strcpy(file_att, val);
        else if (strcasecmp ("RIV", optstr) == 0)
            strcpy(file_riv, val);
        else if (strcasecmp ("RIVSEG", optstr) == 0)
            strcpy(file_rivseg, val);
        else if (strcasecmp ("LAKE", optstr) == 0)
            strcpy(file_lake, val);
        
        /* CALIB PARA DATA */
        else if (strcasecmp ("CALIB", optstr) == 0)
            strcpy(file_calib, val);
        else if (strcasecmp ("PARA", optstr) == 0)
            strcpy(file_para, val);
        else if (strcasecmp ("INIT", optstr) == 0)
            strcpy(file_init, val);
        
        /* PARAMETERS DATA */
        else if (strcasecmp ("LC", optstr) == 0)
            strcpy(file_lc, val);
        else if (strcasecmp ("SOIL", optstr) == 0)
            strcpy(file_soil, val);
        else if (strcasecmp ("GEOL", optstr) == 0)
            strcpy(file_geol, val);
        
        /* TSD DATA */
        else if (strcasecmp ("FORC", optstr) == 0)
            strcpy(file_forc, val);
        else if (strcasecmp ("LAI", optstr) == 0)
            strcpy(file_lai, val);
        else if (strcasecmp ("RL", optstr) == 0)
            strcpy(file_rl, val);
        else if (strcasecmp ("MF", optstr) == 0)
            strcpy(file_mf, val);
        else if (strcasecmp ("LCM", optstr) == 0)
            strcpy(file_lcm, val);
        
        /* init/bnd condition DATA */
        else if (strcasecmp ("ELEBC1", optstr) == 0)
            strcpy(file_ebc1, val);
        else if (strcasecmp ("ELEBC2", optstr) == 0)
            strcpy(file_ebc2, val);
        else if (strcasecmp ("RivBC1", optstr) == 0)
            strcpy(file_rbc1, val);
        else if (strcasecmp ("RivBC2", optstr) == 0)
            strcpy(file_rbc2, val);
        else if (strcasecmp ("LakeBC1", optstr) == 0)
            strcpy(file_lbc1, val);
        else if (strcasecmp ("LakeBC2", optstr) == 0)
            strcpy(file_lbc2, val);
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
    fclose (fp);
}
void FileOut::setOutpath(const char *fn){
    strcpy(outpath, fn);
    updateFilePath();
}
void FileOut::copy(FileOut *p){
    strcpy(   outpath, p->outpath );
    strcpy(   projectname, p->projectname );
    strcpy(   riv_Q_down, p->riv_Q_down );
    strcpy(   riv_Q_surf, p->riv_Q_surf );
    strcpy(   riv_Q_sub, p->riv_Q_sub );
    strcpy(   riv_y_stage, p->riv_y_stage );
    strcpy(   ele_y_snow, p->ele_y_snow );
    strcpy(   ele_y_ic, p-> ele_y_ic);
    strcpy(   ele_y_surf, p->ele_y_surf );
    strcpy(   ele_y_unsat, p->ele_y_unsat );
    strcpy(   ele_y_gw, p->ele_y_gw );
    for(int i = 0; i < 3; i++){
        strcpy(   ele_q_ET[i], p->ele_q_ET[i] );
    }
    strcpy(   ele_q_ETP, p->ele_q_ETP );
    strcpy(   ele_q_ETA, p->ele_q_ETA );
    strcpy(   ele_q_prcp, p->ele_q_prcp );
    strcpy(   ele_q_netprcp, p->ele_q_netprcp );
    strcpy(   ele_Q_surfTot, p->ele_Q_surfTot );
    strcpy(   ele_Q_rsurf, p->ele_Q_rsurf );
    strcpy(   ele_Q_rsub, p->ele_Q_rsub );
    strcpy(   ele_Q_subTot, p->ele_Q_subTot );
    strcpy(   ele_q_infil, p->ele_q_infil );
    strcpy(   ele_q_exfil, p->ele_q_exfil );
    strcpy(   ele_q_rech, p->ele_q_rech );
    
    strcpy(   ewb_q_in, p->ewb_q_in );
    strcpy(   ewb_q_out, p->ewb_q_out );
    strcpy(   ewb_q_io, p->ewb_q_io );
    strcpy(   ewb_dh, p->ewb_dh );
    
    strcpy(   lake_Q_rivin, p->lake_Q_rivin );
    strcpy(   lake_Q_rivout, p->lake_Q_rivout );
    strcpy(   lake_Q_surf, p->lake_Q_surf );
    strcpy(   lake_Q_sub, p->lake_Q_sub );
    strcpy(   lake_y_stage, p->lake_y_stage );
    strcpy(   lake_q_evap, p->lake_q_evap );
    strcpy(   lake_q_prcp, p->lake_q_prcp );
    
    strcpy(   Init_update, p->Init_update );
    strcpy(   Init_bak, p->Init_bak );
    
    strcpy(   floodout, p->floodout );
    strcpy(   obs_sim, p->obs_sim );

}

void FileOut::writeTime(double t){
    long lt = (long)t;
    if( lt % 1440 == 0){
        t1 = clock();
        sec = (double)(t1 - t0) / CLOCKS_PER_SEC;
        t0 = t1;
        fprintf(fid_time, "%.2f \t %.2f \t %.2f\n", t, t / 1440., sec);
//        fflush(fid_time);
    }
}
void FileOut::writeTime(double t, double Percentage, double sec_cpu, double sec_wall, unsigned long NumStep){
    fprintf(fid_time, "%.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %ld \n",
            t, t / 1440., Percentage, sec_cpu, sec_wall, NumStep);
    if(global_fflush_mode){
        fflush(fid_time);
    }
}
