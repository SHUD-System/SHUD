#include "Model_Data.hpp"
void Model_Data::fun_Ele_lakeVertical(int i, double t){
    /*  elements in LAKEs */
    qEleInfil[i] = 0.;
    qEleRecharge[i] = 0.;
    qEleExfil[i] = 0.;
    qEleTrans[i] = 0.;
    qEs[i] = 0.;
    qEu[i] = 0.;
    qEg[i] = 0.;
    qTu[i] = 0.;
    qTg[i] = 0.;
    qEleE_IC[i] = 0.;
    qEleTrans[i] = 0.;
    qEleEvapo[i] = qPotEvap[i];
    qEleETA[i] = qEleE_IC[i] + qEleEvapo[i] + qEleTrans[i];
}
void Model_Data::fun_Ele_lakeHorizon(int i, double t){
    for (int j = 0; j < 3; j++) {
        QeleSurf[i][j] = 0.;
        QeleSub[i][j] = 0.;
    }
}
void Model_Data::fun_Ele_Recharge(int i, double t){
    qEleRecharge[i] = Ele[i].Flux_Recharge(uYus[i] , uYgw[i]);
    qEleRecharge[i] *= fu_Sub[i];
//    CheckNANi(qEleRecharge[i], i, "Model_Data::fun_Ele_Recharge():qEleRecharge[i]");
}

void Model_Data::fun_Ele_Infiltraion(int i, double t){
    Ele[i].Flux_Infiltration(uYsf[i] , uYus[i] , uYgw[i], qEleNetPrep[i]  );
    qEleInfil[i] = Ele[i].u_qi * fu_Surf[i];
    qEleExfil[i] = Ele[i].u_qex * fu_Surf[i];
}
void Model_Data::fun_Ele_surface(int i, double t){
    int j, inabr, ilake;
    double  Ymean, dh, s, CrossA, Q, B;
    double isf, nsf; // Available Y in Surface of this/nabor element
//    isf = uYsf[i] - qEleInfil[i] + qEleExfil[i];
    isf = uYsf[i];
    isf = isf < 0. ? 0. : isf;
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        ilake = Ele[i].lakenabr[j] - 1;
        B = Ele[i].edge[j];
        if(ilake >= 0){  /* For Lake element */
            nsf = yLakeStg[ilake];
            nsf = nsf < 0. ? 0. : nsf;
            Q = WeirFlow_jtoi(lake[ilake].zmin, nsf,
                              Ele[i].z_surf, isf,
                              Ele[i].z_surf, 0.6, B, 0.01); /* func WeirFlow_jtoi is */
            QLakeSurf[ilake] += Q;  /* Positive of QLakeSurf = Element to Lake */
//            CheckNANi( QLakeSurf[ilake] , i, "QLakeSurf[ilake] in Model_Data::fun_Ele_surface");
        }else if (inabr >= 0) {
            /***************************************************************************/
            /* Surface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
//            nsf = uYsf[inabr] - qEleInfil[inabr] + qEleExfil[inabr];
            nsf = uYsf[inabr];
            nsf = nsf < 0. ? 0. : nsf;
            dh = (isf + Ele[i].z_surf) - (nsf + Ele[inabr].z_surf);
            Ymean = avgY_sf(Ele[i].z_surf, isf,
                            Ele[inabr].z_surf, nsf,
                            Ele[i].depression);
            Ymean = min(Ymean, MAXYSURF);/* HARD CODE.When Ymean > 0.5, the solver oscilates; namely, the program slows down dramatically. */
//            Ymean = min(Ymean, fabs(dh));
            if(Ymean <= 0.){
                Q = 0.;
            }else{
                s = dh / Ele[i].Dist2Nabor[j];
                CrossA = Ymean * B;
                if(s > 0 && isf <=0){
                    Q = 0.;
                }else if(s < 0 && nsf <=0){
                    Q = 0.;
                }else{
                    Q= ManningEquation(CrossA, Ele[i].avgRough[j], Ymean, s);
                }
//                CheckNANi(QeleSurf[i][j], i, "QeleSurf[i][j]");
            } //end of ifelse Avg_Y_Surf < 0.
        }else{
            Q = 0;
            if(CS.CloseBoundary){
                /* Void */
            }else{
                if(isf > Ele[i].depression){
                    s = isf / Ele[i].Dist2Edge[j] * 0.5;
                    if(s > 0.){
                        Q = sqrt(s) * cbrt(isf * isf * isf * isf * isf) * B / Ele[i].Rough;
                    }
                }
            }
        } // end of if
        QeleSurf[i][j] = Q;
//        CheckNANi(QeleSurf[i][j], i, "QeleSurf[i][j]");
    } // end of for loop
}// end of functions


void Model_Data::fun_Ele_sub(int i, double t){
    int j, inabr, ilake;
    double  Ymean, dh, Kmean, grad, Q;
    
    for (j = 0; j < 3; j++) {
        inabr = Ele[i].nabr[j] - 1;
        ilake = Ele[i].lakenabr[j] - 1;
        if(ilake >= 0){ /* For Lake element */
            dh = (uYgw[i] + Ele[i].z_bottom) - (yLakeStg[ilake] + lake[ilake].bathymetry.yi[0]);
            if(dh > 0. && uYgw[i] <= 0.02){ /* Depression condition */
                Q = 0.;
            }else if(dh < 0. && yLakeStg[ilake]<= 0.02){ /* Depression condition */
                Q = 0.;
            }else{
                Ymean = avgY_gw(Ele[i].z_bottom, uYgw[i], lake[ilake].bathymetry.yi[0], yLakeStg[ilake], 0.002);
                grad = dh / Ele[i].Dist2Nabor[j];
                /* It should be weighted average. However, there is an ambiguity about distance used */
                Kmean = 0.5 * (Ele[i].u_effKH + Ele[inabr].u_effKH);
                Q = Kmean * grad * Ymean * Ele[i].edge[j];
//                CheckNANi(Q, i, "Q in Model_Data::fun_Ele_sub");
            }
            QLakeSub[ilake] += Q; /* Positive of QLakeSub = Element to Lake */
        }else if (inabr >= 0) {
            /***************************************************************************/
            /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            dh = (uYgw[i] + Ele[i].z_bottom) - (uYgw[inabr] + Ele[inabr].z_bottom);
            if(dh > 0. && uYgw[i] <= 0.02){
                Q = 0.;
            }else if(dh < 0. && uYgw[inabr]<= 0.02){
                Q = 0.;
            }else{
                Ymean = avgY_gw(Ele[i].z_bottom, uYgw[i], Ele[inabr].z_bottom, uYgw[inabr], 0.002);  
                grad = dh / Ele[i].Dist2Nabor[j];
                /* It should be weighted average. However, there is an ambiguity about distance used */
                Kmean = 0.5 * (Ele[i].u_effKH + Ele[inabr].u_effKH);
                Q = Kmean * grad * Ymean * Ele[i].edge[j];
//                CheckNANi(Q, i, "Q in Model_Data::fun_Ele_sub");
            }
        }else {
            Q = 0;
            if(CS.CloseBoundary){
                /* Void */
            }else{
                if(uYgw[i] > Ele[i].depression * 10.){
                    grad = uYgw[i]  / Ele[i].Dist2Edge[j] * 0.5;
                    if(grad > 0.){
                        Q = Ele[i].u_effKH * grad;
                    }
                }
//                CheckNANi(Q, i, "Q in Model_Data::fun_Ele_sub");
            }
        } // end of if
        QeleSub[i][j] = Q * fu_Sub[i];
//        CheckNANi(QeleSub[i][j], i, "Q in Model_Data::fun_Ele_sub");
    } // end of for loop
}// end of functions
