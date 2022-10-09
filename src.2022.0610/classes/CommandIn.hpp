//  CommandIn.hpp
//
//  Created by Lele Shu on 9/29/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef CommandIn_hpp
#define CommandIn_hpp

#include <stdio.h>
#include "Macros.hpp"
#include "IO.hpp"

class CommandIn{
public:
    char    prjname[MAXLEN];
    char    outpath[MAXLEN];
    char    inpath[MAXLEN];
    char    prjfile[MAXLEN];
    char    calibfile[MAXLEN];
    char    dir_cmaes[MAXLEN];
    int     c;
    int     iprj = 0;       /*boolean flag, iprj = 1 read from prj file. iprj = 0, read based on prjname */
    int     iout = 0;   /*User defined output folder*/
    int     n_lambda = 2; /* Number of Threads in OpenMP or MPI */
            CommandIn();
    void    parse(int argc, char **argv);
    void    setFileIO(FileIn *fin, FileOut *fout);
    int     getNumberThreads();
private:
    void    SHUD_help();
};
#endif /* CommandIn_hpp */
