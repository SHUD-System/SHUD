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

/* S1b (openMP #45) — `Model_Data::rhs_flux` is a PURE CARRY-OVER of
 * `Model_Data::f_loop` defined in
 * `SHUD/src/ModelData/MD_f.cpp:11-74` (PR #43 post-state; range
 * includes both #ifdef SHUD_DUMP_RHS probe blocks + after-PassValue
 * no-op hook). The only structural difference vs legacy is the
 * function NAME (`rhs_flux` instead of `f_loop`) — same
 * `Model_Data::` member, same `(double t)` signature, same
 * `this`-relative member / global access (per spec
 * rhs-core-flux-extraction Scenario "Diff legacy f_loop vs new
 * rhs_flux" allowing only function-name + namespace qualifier
 * differences).
 *
 * Process order is the strict 6-step sequence preserved byte-for-byte:
 *   1. Element pass 1  (lake updateLakeElement / fun_Ele_lakeVertical
 *                       + qLakeEvap/qLakePrcp accum; non-lake f_etFlux
 *                       + updateElement + fun_Ele_Infiltraion +
 *                       fun_Ele_Recharge)
 *   2. Element pass 2  (lake fun_Ele_lakeHorizon; non-lake
 *                       fun_Ele_surface + fun_Ele_sub)
 *   3. Segment pass    (fun_Seg_surface + fun_Seg_sub)
 *   4. River pass      (Flux_RiverDown)
 *   5. Lake clamp pass (min/max on qLakeEvap)
 *   6. PassValue()     (zero-reset + segment re-accumulate)
 * Followed by the after-PassValue no-op SHUD_DUMP_RHS hook (legacy).
 *
 * Dump tag strings `"f_loop_before_passvalue"` (MD_f.cpp:67) and
 * `"f_loop"` (MD_f.cpp:72) are preserved verbatim so PR #54's 12
 * `_before_passvalue.bin` + 12 unsuffixed snapshot goldens (4 case
 * × 3 t_values) stay addressable. Renaming = breaks goldens.
 *
 * No sub-function split (rhs_element_vertical / rhs_segment_compute
 * etc. belong to S3, not S1b — per spec Scenario "No sub-function
 * split inside rhs_flux"). */
void Model_Data:: rhs_flux(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            Ele[i].updateLakeElement();
            fun_Ele_lakeVertical(i, t);
            qLakeEvap[Ele[i].iLake - 1] += qEleEvapo[i] / lake[Ele[i].iLake - 1].NumEleLake;
            qLakePrcp[Ele[i].iLake - 1] += qElePrep[i] / lake[Ele[i].iLake - 1].NumEleLake;
        }else{
            f_etFlux(i, t);
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
            fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
        }
    }
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            fun_Ele_lakeHorizon(i, t);
        }else{
            /*========surf/gw flow Function==============*/
            fun_Ele_surface(i, t);  // AFTER infiltration, do the lateral flux. ESP for overland flow.
            fun_Ele_sub(i, t);
        }
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    for (i = 0; i < NumLake; i++) {
        qLakeEvap[i] = min(qLakeEvap[i], qLakePrcp[i] + yLakeStg[i]);
        qLakeEvap[i] = max(0, qLakeEvap[i]);
    }
    /* #43 (S1-pre-B): before-PassValue probe. Dumps Qe2r_Surf
     * (length NumEle), a PassValue() write-set member that carries
     * the previous iteration's element-to-river surface flux state.
     * PassValue (see body at L182-205) zero-resets Qe2r_Surf[0..NumEle-1]
     * and then accumulates QsegSurf over NumSegmt segments; capturing
     * Qe2r_Surf HERE gives a deterministic snapshot of the value that
     * is about to be cleared + re-derived by PassValue. PR #54 round-1
     * fix F4 replaced the prior QeleSurfTot probe payload (which
     * f_update zero-resets so the snapshot was always all zeros and
     * thus useless as a before-vs-after PassValue diff) with this
     * write-set member. Site tag "f_loop_before_passvalue" is distinct
     * from the no-op "f_loop" hook below + the "f_update" hook in
     * MD_update.cpp:151; the writer SHUD_DUMP_FNAME_SUFFIX env
     * disambiguates output files (`snapshot_t<v>_before_passvalue.bin`
     * vs `snapshot_t<v>.bin`). SHUD_DUMP_RHS=0 builds emit zero code
     * (compile-switch neutrality contract). */
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop_before_passvalue", t, Qe2r_Surf, NumEle);
#endif
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop", t, NULL, 0);
#endif
}

void Model_Data::rhs_core(double *Y, double *DY, double t){
    /* S1b mixed-mode: `rhs_update` (S1a) + `rhs_flux` (S1b) are
     * migrated to the new path; `rhs_apply` is still the legacy
     * `f_applyDY` fallback (S1c replaces that final fallback). */
    rhs_update(Y, DY, t);
    rhs_flux(t);
    f_applyDY(DY, t);
}
