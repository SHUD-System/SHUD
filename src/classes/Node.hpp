//  Node.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//
#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include "Macros.hpp"
class _Node {    /* Data model for a node */
public:
    int index;    /* Node no. */
    double x = NA_VALUE;    /* x coordinate */
    double y = NA_VALUE;    /* y coordinate */
    double zmin = NA_VALUE;    /* z bed rock elevation */
    double AqD = NA_VALUE;    /* Aquifer Depth */
    double zmax = NA_VALUE;    /* z surface elevation */
    
    _Node();
    void Init();
    void Init(double cAqD);
};

#endif /* Node_hpp */
