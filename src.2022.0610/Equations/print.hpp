//  print.hpp
//
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef print_hpp
#define print_hpp
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Macros.hpp"
#include "IO.hpp"
#include "functions.hpp"
#include "Model_Data.hpp"

void PrintDataNew (Print_Ctrl PCtrl, double tmpt, double dt, Control_Data *CS);
void SHUDlogo(void);

#endif /* print_hpp */
