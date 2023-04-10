//  IO.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef IO_hpp
#define IO_hpp

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <time.h>
#include "Macros.hpp"
#include "functions.hpp"
#include "funPlatform.hpp"

class FileIn {
public:
    int numthreads = 0;
    /*Model io*/
    char projectname[MAXLEN];
    char outpath[MAXLEN];
    char inpath[MAXLEN];
    char logfile[MAXLEN];
    
    /*Spatial Data*/
    char file_mesh[MAXLEN];
    char file_att[MAXLEN];
    char file_riv[MAXLEN];
    char file_rivseg[MAXLEN];
    char file_lake[MAXLEN];
    char file_lake_ic[MAXLEN];
    char file_lake_bathy[MAXLEN];
    
    /* model configuration */
    char file_init[MAXLEN];
    char file_para[MAXLEN];
    char file_calib[MAXLEN];
    char file_cfgout[MAXLEN];
    
    /* physical parameters */
    char file_lc[MAXLEN];
    char file_soil[MAXLEN];
    char file_geol[MAXLEN];
    
    /* Time-series data */
    char file_ebc1[MAXLEN];
    char file_ebc2[MAXLEN];
    char file_rbc1[MAXLEN];
    char file_rbc2[MAXLEN];
    char file_lbc1[MAXLEN];
    char file_lbc2[MAXLEN];
    
    char file_forc[MAXLEN];
    char file_lai[MAXLEN];
    char file_rl[MAXLEN];
    char file_mf[MAXLEN];
    char file_lcm[MAXLEN];
    
    /* Time-series data, observation data for calibration*/
    char file_obs[MAXLEN];
    
    void setInFilePath(char * inpath, char *  projectname);
    void setInFilePath(char * indir, char *  pjrname, int nthred);
    void readProject(const char *fn); /* read project file from file fn*/
    void saveProject(); /* save project file to outdir */
    void setCalibFile(const char *fn);
    void setOutpath(const char *fn);
    void copy(FileIn *fin);
}               ;
class FileOut{
public:
    char suffix[MAXLEN];
    char outpath[MAXLEN];
    char projectname[MAXLEN];
    //rivers
    char riv_Q_down[MAXLEN];
    char riv_Q_up[MAXLEN];
    char riv_Q_surf[MAXLEN];
    char riv_Q_sub[MAXLEN];
    char riv_y_stage[MAXLEN];
    //cells
    char ele_y_snow[MAXLEN];
    char ele_y_ic[MAXLEN];
    char ele_y_surf[MAXLEN];
    char ele_y_unsat[MAXLEN];
    char ele_y_gw[MAXLEN];
    
    //cell-fluxes
    char ele_q_ET[3][MAXLEN];
    char ele_q_ETP[MAXLEN];
    char ele_q_ETA[MAXLEN];
    char ele_q_prcp[MAXLEN];
    char ele_q_netprcp[MAXLEN];
//    char ele_Q_surf[MAXLEN];
//    char ele_Q_sub[MAXLEN];
    char ele_Q_surfTot[MAXLEN];
    char ele_Q_subTot[MAXLEN];
    char ele_Q_sub0[MAXLEN];
    char ele_Q_sub1[MAXLEN];
    char ele_Q_sub2[MAXLEN];
    char ele_Q_surf0[MAXLEN];
    char ele_Q_surf1[MAXLEN];
    char ele_Q_surf2[MAXLEN];
    
    char ele_Q_rsurf[MAXLEN]; // Element to River via Surface
    char ele_Q_rsub[MAXLEN]; // Element to River via Subsurface
    char ele_q_infil[MAXLEN];
    char ele_q_exfil[MAXLEN]; // Exfiltration
    char ele_q_rech[MAXLEN];
    
    //cell_wb
    char ewb_q_in[MAXLEN];
    char ewb_q_out[MAXLEN];
    char ewb_q_io[MAXLEN];
    char ewb_dh[MAXLEN];
    
    //Lakes
    char lake_Q_rivin[MAXLEN];
    char lake_Q_rivout[MAXLEN];
    char lake_Q_surf[MAXLEN];
    char lake_Q_sub[MAXLEN];
    char lake_y_stage[MAXLEN];
    char lake_a_area[MAXLEN];
    char lake_q_evap[MAXLEN];
    char lake_q_prcp[MAXLEN];
    
    //Init_Condition
    char Init_update[MAXLEN];
    char Init_bak[MAXLEN];
    
    //calib backup
    char Calib_bak[MAXLEN];
    
    //Time steps
    FILE *fid_time;
    char File_Time[MAXLEN];
    
    // Flood Warnings
    char floodout[MAXLEN];
    
    char obs_sim[MAXLEN];
    
    FileOut();
    void setsuffix(const char *s);
    void setOutFilePath( char *outpath, char * prjname);
    void updateFilePath();
    void createDir();
    void setOutpath(const char *fn);
    void writeTime(double t);
    void writeTime(double t, double Percentage, double sec_cpu, double sec_wall, unsigned long NumStep);
    void copy(FileOut *fout);
private:
    clock_t t0 = clock();
    clock_t t1 = t0;
    double sec;
};

#endif /* IO_h */
