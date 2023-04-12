//  CommandIn.cpp
//
//  Created by Lele Shu on 9/29/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "CommandIn.hpp"
void CommandIn::SHUD_help(void ){
    printf ("\n\nUsage:\n");
    printf ("./shud [-0gv] [-p project_file] [-c Calib_file] [-o output] [-n Num_Threads] <project_name>\n\n");
    printf (" -0 Dummy simulation. Load input and write output, but no calculation.\n");
    printf (" -f fflush for each time interval. fflush export data frequently, but slow down performance on cluster.\n");
    printf (" -g Sequential coupled Surface-Unsaturated-Saturated-River mode.\n");
    printf (" -v Dummy simulation. Load input and write output, but no calculation.\n");
    printf (" -c Calibration file (.cfg.calib). \n");
    printf (" -o output folder. Default is output/projname.out\n");
    printf (" -p projectfile, which includes the path to input files and output path.\n");
    printf (" -n Number of threads to run with OpenMP. \n");
}

void CommandIn::parse(int argc, char **argv){
    if(argc<=1){
        SHUD_help();
        myexit(ERRSUCCESS);
    }
    while ((c = getopt (argc, argv, "0fgvc:e:n:o:p:")) != -1){
        switch (c){
            case '0':
                dummy_mode = 1;
                break;
            case 'f':
                global_fflush_mode = 1;
                break;
            case 'g':
                global_implicit_mode = 0;
                break;
            case 'v':
                global_verbose_mode = 1;
                break;
            case 'c':
                strcpy(calibfile, optarg);
                break;
            case 'e':
                strcpy(dir_cmaes, optarg);
                break;
            case 'o':
                strcpy(outpath, optarg);
                iout = 1;
                break;
            case 'n':
                n_lambda = atoi(optarg) ;
                break;
            case 'p':
                strcpy(prjfile, optarg);
                iprj = 1;
                break;
            case '?':
                if (optopt == 'p')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",
                             optopt);
                myexit(-1);
            default:
                break;
        }
    }
    if(iprj){
       // void
    }else{
        strcpy(prjname, argv[optind]);
    }
    
    
#ifdef DEBUG
    printf("\t\t\t * Debug mode enable.\n");
#endif
    
#ifdef _OPENMP_ON
    printf("\t\t * openMP enabled. Maximum Threads = %d\n", omp_get_max_threads());
#else
    printf("\t\t * openMP disabled.\n");
#endif
}
CommandIn::CommandIn(){
    strcpy(prjname, "");
    strcpy(outpath, "");
    strcpy(inpath, "");
    strcpy(prjfile, "");
    strcpy(calibfile, "");
    strcpy(dir_cmaes, "cmaes");
}
void CommandIn::setFileIO(FileIn *fin, FileOut *fout){
    if(iprj){
        fin->readProject(prjfile);
    }else{
        sprintf(inpath, "input/%s", prjname);
        fin->setInFilePath(inpath, prjname, n_lambda);
        if (iout){
            fin->setOutpath(outpath);
        }
    }
    if(calibfile[0] != '\0' ){
        fin->setCalibFile(calibfile);
    }
    fout->setOutFilePath(fin->outpath, fin->projectname);
    fin->saveProject();
}
int CommandIn::getNumberThreads(){
    if(n_lambda < 0 ){
        return 0;
    }else{
        return n_lambda;
    }
}
