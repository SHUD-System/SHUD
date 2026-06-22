//
//  AccTemperature.hpp
//  shud
//
//  Created by Lele Shu on 3/29/21.
//  Copyright © 2021 Lele Shu. All rights reserved.
//

#ifndef AccTemperature_hpp
#define AccTemperature_hpp

#include <stdio.h>
#include <queue>

class _AccTemp{
private:
    double  Time_start = -9999.;
    double  T_AccDay = 0.;
    int     N_of_day;
public:
    std::queue<double> que;
    int     MaxLen = 28;
    double  ACC;
    void setLength(int x){
        MaxLen = x;
    }
private:
    void push(double x){
        que.push(x);        /* PUSH the new value */
        
        ACC += x;           /* Add the new value */
        if(que.size() > MaxLen){
            ACC -= que.front(); /* Remove the OLDEST value */
            que.pop();
        }else{
            /* void */
        }
        
    };
public:
    
    _AccTemp(){
        Time_start = -9999.;
        T_AccDay = 0.;
        N_of_day = 0;
    };
    
    void push(double x, double tnow){
        T_AccDay += x;
        N_of_day++;
        if( (tnow - Time_start) >= 1440. ){
            push(T_AccDay /  N_of_day);  /* Push the mean value */
            T_AccDay = 0.;
            N_of_day = 0;
            Time_start = tnow;
        }else{
            
        }
    }
    double getACC(){
        /* S6b.1 (#184): defensive divide-zero guard. In the current call
         * graph (MD_ET.cpp:153-156, push() immediately precedes getACC()),
         * Time_start = -9999. (L17 class init) guarantees the very first
         * push(x, tnow) call enqueues — for any tnow >= -8559, the
         * (tnow - Time_start) >= 1440 condition holds — so the empty-queue
         * path is unreachable and this fix is bitwise-neutral on all
         * current goldens (B1a-tag SHA256 PASS). Guards against future
         * call-graph changes that might invoke getACC() before any push.
         * 0.0 fallback per master plan §4.12 / §S2.15. */
        return que.empty() ? 0.0 : ACC / que.size();
    }
};
#endif /* AccTemperature_hpp */
