//  TabularData.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef TabularData_hpp
#define TabularData_hpp

#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "Macros.hpp"
#include "functions.hpp"
using namespace std;
class TabularData {
public:
    int nrow = 0;
    int ncol = 0;
    double **x;
    char header[MAXLEN];
    
    TabularData();
    ~TabularData();
    int read(const char *fn);
    int read(FILE *fp);
    int read(FILE *fp, int *NumCol);
private:
    void reset();
};
#endif /* TabularData_hpp */
