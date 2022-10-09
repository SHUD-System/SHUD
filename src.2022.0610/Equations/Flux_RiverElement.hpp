//
//  Flux_Subsurface.hpp
//  SHUD
//
//  Created by Lele Shu on 1/18/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#ifndef Flux_Subsurface_hpp
#define Flux_Subsurface_hpp

#include <stdio.h>
#include "Macros.hpp"
#include "Equations.hpp"
double flux_R2E_SF(double yr, double zr,
                   double ye, double ze,
                   double rough, double Len, double dist, double threshold);
double flux_R2E_GW(double yr, double zr,
                   double ye, double ze,
                   double Kele, double Kriv, 
                   double L, double D_riv);

#endif /* Flux_Subsurface_hpp */
