//
//  funPlatform.hpp
//  SHUD
//
//  Created by Lele Shu on 4/10/23.
//  Copyright © 2023 Lele Shu. All rights reserved.
//

#ifndef funPlatform_hpp
#define funPlatform_hpp

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>

#include "Macros.hpp"

//  Windows
#ifdef _WIN32
#include <Windows.h>

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
#endif

double get_wall_time();
double get_cpu_time();

void mkdir_p( char *dir, int mode);

#endif /* funPlatform_hpp */
