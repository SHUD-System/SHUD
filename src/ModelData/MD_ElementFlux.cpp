#include "Model_Data.hpp"
void Model_Data::fun_Ele_Recharge(int i, double t){
    qEleRecharge[i] = Ele[i].Flux_Recharge(uYus[i] , uYgw[i]);
//    qEleRecharge[i] = Ele[i].u_qr;
}

void Model_Data::fun_Ele_Infiltraion(int i, double t){
    Ele[i].Flux_Infiltration(uYsf[i] , uYus[i] , uYgw[i], qEleNetPrep[i]  );
    qEleInfil[i] = Ele[i].u_qi;
    qEleExfil[i] = Ele[i].u_qex * vExf;
}
void Model_Data::fun_Ele_surface(int i, double t){
    int j, inabr;
    double  Ymean, dh, s, CrossA, Q;
    double isf, nsf; // Available Y in Surface of this/nabor element    
    isf = uYsf[i] - qEleInfil[i] + qEleExfil[i];
    isf = isf < 0. ? 0. : isf;
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        if (inabr >= 0) {
            /***************************************************************************/
            /* Surface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            nsf = uYsf[inabr] - qEleInfil[inabr] + qEleExfil[inabr];
            nsf = nsf < 0. ? 0. : nsf;
            dh = (isf + Ele[i].zmax) - (nsf + Ele[inabr].zmax);
            Ymean = avgY_sf(Ele[i].zmax, isf,
                            Ele[inabr].zmax, nsf,
                            Ele[i].depression);
            if(Ymean <= 0.){
                Q = 0.;
            }else{
                s = dh / Ele[i].Dist2Nabor[j];
                CrossA = Ymean * Ele[i].edge[j];
                if(s > 0 && isf <=0){
                    Q = 0.;
                }else if(s < 0 && nsf <=0){
                    Q = 0.;
                }else{
                    Q= ManningEquation(CrossA, Ele[i].avgRough[j], Ymean, s);
                }
//                CheckNANi(QeleSurf[i][j], i, "QeleSurf[i][j]");
            } //end of ifelse Avg_Y_Surf < 0.
        } else {
            Q= 0;
        } // end of if
        QeleSurf[i][j] = Q;
    } // end of for loop
}// end of functions


void Model_Data::fun_Ele_sub(int i, double t){
    int j, inabr;
    double  Ymean, dh, Kmean, grad, Q;
    
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        if (inabr >= 0) {
            /***************************************************************************/
            /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            dh = (uYgw[i] + Ele[i].zmin) - (uYgw[inabr] + Ele[inabr].zmin);
            if(dh > 0. && uYgw[i] <= EPSILON){
                Q = 0.;
            }else if(dh < 0. && uYgw[inabr]<= EPSILON){
                Q = 0.;
            }else{
                Ymean = avgY_gw(Ele[i].zmin, uYgw[i], Ele[inabr].zmin, uYgw[inabr], 0.002);  
                grad = dh / Ele[i].Dist2Nabor[j];
                /* It should be weighted average. However, there is an ambiguity about distance used */
                Kmean = 0.5 * (Ele[i].u_effKH + Ele[inabr].u_effKH);
                Q = Kmean * grad * Ymean * Ele[i].edge[j];
//                CheckNA(Q, "Q in Model_Data::fun_Ele_sub");
            }
        } else {
            Q = 0;
        } // end of if
        QeleSub[i][j] = Q;
    } // end of for loop
}// end of functions
