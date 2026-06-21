#include "Model_Data.hpp"
#include <cassert>
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
    /* S5d.2-5a (#179) — write via QeleSurfAt / QeleSubAt accessors;
     * jagged 2-D `QeleSurf[i][j]` retired. Row-major flat: at(i,j)
     * ↔ `_flat[3*i + j]`. */
    for (int j = 0; j < 3; j++) {
        QeleSurfAt(i, j) = 0.;
        QeleSubAt(i, j) = 0.;
    }
}
void Model_Data::fun_Ele_Recharge(int i, double t){
    /* S5d.1 (#178) — Flux_Recharge writes Ele[i].u_qr (and reads many
     * private fields of _Element). It is a _Element member method, so
     * the call site keeps AoS dispatch; sync_hot_dynamic(i) is invoked
     * AFTER to refresh the u_* SoA mirror (u_qr is NOT in the SoA
     * roster — grep showed 0 RHS reads of u_qr in MD_ElementFlux/MD_f/
     * MD_ET — but updateElement / Flux_Infiltration / Flux_Recharge all
     * touch the writer set so the sync call is uniform). */
    qEleRecharge[i] = Ele[i].Flux_Recharge(uYus[i] , uYgw[i]);
    sync_hot_dynamic(i);
    qEleRecharge[i] *= fu_Sub[i];
//    CheckNANi(qEleRecharge[i], i, "Model_Data::fun_Ele_Recharge():qEleRecharge[i]");
}

void Model_Data::fun_Ele_Infiltraion(int i, double t){
    /* S5d.1 (#178) — Flux_Infiltration writes Ele[i].u_qi, u_qex (and
     * u_effkInfi privately). Call sync_hot_dynamic AFTER so subsequent
     * SoA reads of u_qi / u_qex in this RHS step see the just-written
     * AoS values. */
    Ele[i].Flux_Infiltration(uYsf[i] , uYus[i] , uYgw[i], qEleNetPrep[i]  );
    sync_hot_dynamic(i);
    qEleInfil[i] = hot.u_qi[i] * fu_Surf[i];
    qEleExfil[i] = hot.u_qex[i] * fu_Surf[i];
}
void Model_Data::fun_Ele_surface(int i, double t){
    /* S5d.1 (#178) — all Ele[i].<hot-field> reads rerouted to hot.<field>
     * SoA mirror. Cross-element reads (Ele[inabr].z_surf) likewise.
     * No AoS writes in this function; SoA mirror is up-to-date because
     * initialize_hot() + sync_hot_dynamic invariants hold. */
    int j, inabr, ilake;
    double  Ymean, dh, s, CrossA, Q, B;
    double isf, nsf; // Available Y in Surface of this/nabor element
//    isf = uYsf[i] - qEleInfil[i] + qEleExfil[i];
    isf = uYsf[i];
    isf = isf < 0. ? 0. : isf;
    for (j = 0; j < 3; j++) {
        inabr = hot.nabr_flat[3*i + j] - 1;
        ilake = hot.lakenabr_flat[3*i + j] - 1;
        B = hot.edge_flat[3*i + j];
        if(ilake >= 0){  /* For Lake element */
            nsf = yLakeStg[ilake];
            nsf = nsf < 0. ? 0. : nsf;
            Q = WeirFlow_jtoi(lake[ilake].zmin, nsf,
                              hot.z_surf[i], isf,
                              hot.z_surf[i], 0.6, B, 0.01); /* func WeirFlow_jtoi is */
            /* S3b.2 (PR-9): shared write `QLakeSurf[ilake] += Q` replaced
             * with deterministic per-edge slot. PassValue_legacy() will gather
             * QeleSurf_lake -> QLakeSurf. Will be replaced by
             * rhs_deterministic_gather() in S3c (PR-11). */
            QeleSurf_lake[i*3 + j] = Q;
//            CheckNANi( QLakeSurf[ilake] , i, "QLakeSurf[ilake] in Model_Data::fun_Ele_surface");
        }else if (inabr >= 0) {
            /***************************************************************************/
            /* Surface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
//            nsf = uYsf[inabr] - qEleInfil[inabr] + qEleExfil[inabr];
            nsf = uYsf[inabr];
            nsf = nsf < 0. ? 0. : nsf;
            dh = (isf + hot.z_surf[i]) - (nsf + hot.z_surf[inabr]);
            Ymean = avgY_sf(hot.z_surf[i], isf,
                            hot.z_surf[inabr], nsf,
                            hot.depression[i]);
            Ymean = min(Ymean, MAXYSURF);/* HARD CODE.When Ymean > 0.5, the solver oscilates; namely, the program slows down dramatically. */
//            Ymean = min(Ymean, fabs(dh));
            if(Ymean <= 0.){
                Q = 0.;
            }else{
                s = dh / hot.Dist2Nabor_flat[3*i + j];
                CrossA = Ymean * B;
                if(s > 0 && isf <=0){
                    Q = 0.;
                }else if(s < 0 && nsf <=0){
                    Q = 0.;
                }else{
                    Q= ManningEquation(CrossA, hot.avgRough_flat[3*i + j], Ymean, s);
                }
//                CheckNANi(QeleSurf[i][j], i, "QeleSurf[i][j]");
            } //end of ifelse Avg_Y_Surf < 0.
        }else{
            Q = 0;
            if(CS.CloseBoundary){
                /* Void */
            }else{
                if(isf > hot.depression[i]){
                    s = isf / hot.Dist2Edge_flat[3*i + j] * 0.5;
                    if(s > 0.){
                        Q = sqrt(s) * cbrt(isf * isf * isf * isf * isf) * B / hot.Rough[i];
                    }
                }
            }
        } // end of if
        /* S5d.2-5a (#179) — flat write via accessor. */
        QeleSurfAt(i, j) = Q;
//        CheckNANi(QeleSurfAt(i, j), i, "QeleSurfAt(i, j)");
    } // end of for loop
}// end of functions


void Model_Data::fun_Ele_sub(int i, double t){
    /* S5d.1 (#178) — all Ele[i].<hot-field> and Ele[inabr].<hot-field>
     * reads rerouted to hot.<field>[<idx>] SoA mirror. Identical
     * algorithm; reads only. */
    int j, inabr, ilake;
    double  Ymean, dh, Kmean, grad, Q;

    for (j = 0; j < 3; j++) {
        inabr = hot.nabr_flat[3*i + j] - 1;
        ilake = hot.lakenabr_flat[3*i + j] - 1;
        if(ilake >= 0){ /* For Lake element */
            assert(inabr >= 0);
            dh = (uYgw[i] + hot.z_bottom[i]) - (yLakeStg[ilake] + lake[ilake].bathymetry.yi[0]);
            if(dh > 0. && uYgw[i] <= 0.02){ /* Depression condition */
                Q = 0.;
            }else if(dh < 0. && yLakeStg[ilake]<= 0.02){ /* Depression condition */
                Q = 0.;
            }else{
                Ymean = avgY_gw(hot.z_bottom[i], uYgw[i], lake[ilake].bathymetry.yi[0], yLakeStg[ilake], 0.002);
                grad = dh / hot.Dist2Nabor_flat[3*i + j];
                /* It should be weighted average. However, there is an ambiguity about distance used */
                Kmean = 0.5 * (hot.u_effKH[i] + hot.u_effKH[inabr]);
                Q = Kmean * grad * Ymean * hot.edge_flat[3*i + j];
//                CheckNANi(Q, i, "Q in Model_Data::fun_Ele_sub");
            }
            /* S3b.3 (PR-9): shared write `QLakeSub[ilake] += Q` replaced
             * with deterministic per-edge slot. PassValue_legacy() will gather
             * QeleSub_lake -> QLakeSub. Will be replaced by
             * rhs_deterministic_gather() in S3c (PR-11). */
            QeleSub_lake[i*3 + j] = Q;
        }else if (inabr >= 0) {
            /***************************************************************************/
            /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
            /***************************************************************************/
            dh = (uYgw[i] + hot.z_bottom[i]) - (uYgw[inabr] + hot.z_bottom[inabr]);
            if(dh > 0. && uYgw[i] <= 0.02){
                Q = 0.;
            }else if(dh < 0. && uYgw[inabr]<= 0.02){
                Q = 0.;
            }else{
                Ymean = avgY_gw(hot.z_bottom[i], uYgw[i], hot.z_bottom[inabr], uYgw[inabr], 0.002);
                grad = dh / hot.Dist2Nabor_flat[3*i + j];
                /* It should be weighted average. However, there is an ambiguity about distance used */
                Kmean = 0.5 * (hot.u_effKH[i] + hot.u_effKH[inabr]);
                Q = Kmean * grad * Ymean * hot.edge_flat[3*i + j];
//                CheckNANi(Q, i, "Q in Model_Data::fun_Ele_sub");
            }
        }else {
            Q = 0;
            if(CS.CloseBoundary){
                /* Void */
            }else{
                if(uYgw[i] > hot.depression[i] * 10.){
                    grad = uYgw[i]  / hot.Dist2Edge_flat[3*i + j] * 0.5;
                    if(grad > 0.){
                        Q = hot.u_effKH[i] * grad;
                    }
                }
//                CheckNANi(Q, i, "Q in Model_Data::fun_Ele_sub");
            }
        } // end of if
        /* S5d.2-5a (#179) — flat write via accessor. */
        QeleSubAt(i, j) = Q * fu_Sub[i];
//        CheckNANi(QeleSubAt(i, j), i, "Q in Model_Data::fun_Ele_sub");
    } // end of for loop
}// end of functions
