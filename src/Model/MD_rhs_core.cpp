/* MD_rhs_core.cpp — S1a RHS core scaffolding (openMP issue #44).
 *
 * `Model_Data::rhs_update()` is a PURE CARRY-OVER of
 * `Model_Data::f_update` defined in
 * SHUD/src/ModelData/MD_update.cpp:63-153. The only structural
 * difference vs legacy is the function NAME (`rhs_update` instead of
 * `f_update`) — same `Model_Data::` member, same signature, same
 * `this`-relative member / global access (per spec
 * rhs-core-scaffolding Scenario "Source carry-over diff is
 * structural-only" allowing function-name + namespace + qualifier
 * differences only).
 *
 * All variable names, loop bounds, branch predicates, expression
 * trees, floating-point operation order, and global-state read /
 * write timing are preserved byte-for-byte.
 *
 * `SHUD_DUMP_RHS` hook tag string is `"f_update"` (NOT `"rhs_update"`)
 * — rename would break snapshot lookup against the 24 goldens shipped
 * in PR #53 + #54. See spec Scenario "Source carry-over diff is
 * structural-only" final bullet.
 *
 * `timeNow` is NOT written here — it is assigned at f.cpp::f() L22
 * before `rhs_core()` dispatch (single-source per spec Scenario
 * "`timeNow` not double-written").
 */
#include "MD_rhs_core.hpp"
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif

void Model_Data::rhs_update(double *Y, double *DY, double t){
    for (int i = 0; i < NumEle; i++) {
//        uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
//        uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
        for(int j = 0; j < 3; j++){
            QeleSub[i][j] = 0.;
            QeleSurf[i][j] = 0.;
            QeleSubTot[i] = 0.;
            QeleSurfTot[i] = 0.;
        }
        uYsf[i] = Y[iSF];
        uYus[i] = Y[iUS];
        if(Ele[i].iBC == 0){ // NO BC
//            uYgw[i] = max(0.0, Y[iGW]);
            uYgw[i] = Y[iGW];
            Ele[i].QBC = 0.;
        }else if(Ele[i].iBC > 0){ // BC fix head
            Ele[i].yBC = tsd_eyBC.getX(t, Ele[i].iBC);
            uYgw[i] = Ele[i].yBC;
            Ele[i].QBC = 0.;
        }else{ // BC fix flux to GW
            Ele[i].QBC = tsd_eqBC.getX(t, -Ele[i].iBC);
        }
        qEleExfil[i] = 0.;
        qEleInfil[i] = 0.;
        /***** SS and BC *****/
//        for(int j = 0; j<3;j++){
//            Ele[i].iupdSF[j] = 0;
//            Ele[i].iupdGW[j] = 0;
//        }
/********* Below are remove because the bass-balance issue. **********/
//        for (int j = 0; j < 3; j++) {
//            if(Ele[i].nabr[j] > 0){
//                Ele[i].surfH[j] = (Ele[Ele[i].nabr[j] - 1].zmax + uYsf[Ele[i].nabr[j] - 1]);
//            }else{
//                Ele[i].surfH[j] = (Ele[i].zmax + uYsf[i]);
//            }
//        }
//        Ele[i].dhBYdx = dhdx(Ele[i].surfX, Ele[i].surfY, Ele[i].surfH);
//        Ele[i].dhBYdy = dhdy(Ele[i].surfX, Ele[i].surfY, Ele[i].surfH);
//        Ele[i].Avg_Sf = sqpow2(Ele[i].dhBYdx, Ele[i].dhBYdy);
    }//end of for j=1:NumEle

    for (int i = 0; i < NumRiv; i++ ){
        uYriv[i] = Y[iRIV];
        /* qrivsurf and qrivsub are calculated in Element fluxes.
         qrivDown and qrivUp are calculated in River fluxes. */
        Riv[i].updateRiver(uYriv[i]);
        /***** SS and BC *****/
        Riv[i].qBC = 0.0;
        if(Riv[i].BC == 0){
            /* Void */
        }else if(Riv[i].BC < 0){ // Fixed Flux INTO river Reaches.
            Riv[i].qBC = tsd_rqBC.getX(t, -Riv[i].BC);
        }else if (Riv[i].BC > 0){ // Fixed Stage of river reach.
            Riv[i].yBC = tsd_ryBC.getX(t, Riv[i].BC);
            uYriv[i] = Riv[i].yBC;
        }
#ifdef DEBUG
        CheckNANi(uYriv[i], i, "uYriv in f_update.");
#endif
    }

    for (int i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    for (int i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }
    for (int i = 0; i < NumLake; i++) {
        yLakeStg[i] = Y[iLAKE];
        lake[i].yStage = yLakeStg[i];
        lake[i].update();
        y2LakeArea[i] = lake[i].u_toparea;
        QLakeSub[i] = 0.;
        QLakeSurf[i] = 0.;
        qLakeEvap[i] = 0.;
        qLakePrcp[i] = 0.;
        QLakeRivIn[i] = 0.;
        QLakeRivOut[i] = 0.;
    }
    for (int i = 0; i < NumY; i++){
        DY[i] = 0.;
    }
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_update", t, DY, NumY);
#endif
}

void Model_Data::rhs_core(double *Y, double *DY, double t){
    /* S1a mixed-mode: only `rhs_update` is migrated; flux and apply
     * still call into legacy `f_loop` / `f_applyDY`. S1b / S1c
     * replace these fallbacks with `rhs_flux` / `rhs_apply`. */
    rhs_update(Y, DY, t);
    f_loop(t);
    f_applyDY(DY, t);
}
