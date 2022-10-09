//
//  Node.cpp
//  SHUD
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Node.hpp"
_Node::_Node(){
    
}
void _Node::Init(){
    Init(0.0);
}
void _Node::Init(double cAqD){
    zmin = zmax - (AqD + cAqD);
    if(AqD < 1.){
        fprintf(stderr, "WARNING:: Aqd of Node(%d) = %f\n",index, AqD);
        fprintf(stderr, "Press anykey to continue ...\n");
//        getchar();
    }
}
