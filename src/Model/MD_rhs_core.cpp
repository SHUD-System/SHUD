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
#include "MD_adjacency.hpp"  /* S4 PR-10 (#154): 7 adjacency lists +
                              * `build_adjacency_lists()` declarations.
                              * Lists are built once from
                              * `Model_Data::initialize()` (MD_initialize.cpp);
                              * PR-10 BUILDS but does NOT YET USE them in
                              * rhs_core (PR-11 / S3c will replace PassValue's
                              * in-loop gather). The include here satisfies
                              * spec s4-adjacency-topology Scenario
                              * "MD_rhs_core.cpp 文件顶部 SHALL 包含
                              * #include MD_adjacency.hpp". */
#include <cstdlib>   /* std::abort -- S1d.1 OMP stub regression guard */
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
            /* S3b.4 (PR-9): shared writes
             *   qLakeEvap[Ele[i].iLake-1] += qEleEvapo[i] / NumEleLake
             *   qLakePrcp[Ele[i].iLake-1] += qElePrep[i]  / NumEleLake
             * extracted into deterministic per-element slots. The
             * division (by lake.NumEleLake) is now per-element. Gather
             * (after RivLoop, BEFORE the lake clamp below) sums to
             * per-lake qLakeEvap / qLakePrcp. */
            qEleEvapo_lake[i] = qEleEvapo[i] / lake[Ele[i].iLake - 1].NumEleLake;
            qElePrep_lake[i]  = qElePrep[i]  / lake[Ele[i].iLake - 1].NumEleLake;
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
    /* S3b.4 (PR-9): transitional gather per-element -> per-lake.
     * Must run BEFORE the lake clamp below (clamp reads qLakeEvap /
     * qLakePrcp). Cannot live in PassValue because PassValue() is
     * called AFTER the clamp. Will be replaced by
     * rhs_deterministic_gather() in S3c (PR-11). */
    if(lakeon){
        for (i = 0; i < NumLake; i++) {
            qLakeEvap[i] = 0.;
            qLakePrcp[i] = 0.;
        }
        for (i = 0; i < NumEle; i++) {
            if(Ele[i].iLake > 0){
                int ilake = Ele[i].iLake - 1;
                qLakeEvap[ilake] += qEleEvapo_lake[i];
                qLakePrcp[ilake] += qElePrep_lake[i];
            }
        }
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

/* S1c (openMP #46) — `Model_Data::rhs_apply` is a PURE CARRY-OVER of
 * `Model_Data::f_applyDY` defined in
 * `SHUD/src/ModelData/MD_f.cpp:76-182`. The only structural difference
 * vs legacy is the function NAME (`rhs_apply` instead of `f_applyDY`)
 * — same `Model_Data::` member, same `(double *DY, double t)`
 * signature, same `this`-relative member / global access (per spec
 * rhs-core-applydy-extraction Scenario "Diff vs B0 f_applyDY shows
 * only signature changes" allowing only function-name + namespace
 * qualifier differences).
 *
 * River DY uses serial 3-step formula (length + area clamp +
 * `fun_dAtodY()`); the simpler `_omp` "/ u_TopArea" 1-step formula
 * MUST NOT be used (per spec Scenario "River DY divisor is
 * Riv[i].Length not u_TopArea"). Lake DY divide-by-zero UB is
 * preserved verbatim — no `y2LakeArea[i] == 0` guard / epsilon /
 * ternary (per spec Scenario "Lake DY divide-by-zero UB preserved";
 * numerical-stability fixes deferred to S5 / S6).
 *
 * Commented-out debug printf blocks preserved verbatim.
 *
 * `SHUD_DUMP_RHS` hook tag string is `"f_applyDY"` (NOT `"rhs_apply"`)
 * — rename would break snapshot lookup against any goldens shipped in
 * PR #53 / #54 that key on the legacy tag. See spec Scenario "Diff
 * vs B0 f_applyDY shows only signature changes" final bullet. */
void Model_Data::rhs_apply(double *DY, double t){
    double area;
    int isf, ius, igw;
    for (int i = 0; i < NumEle; i++) {
        isf = iSF; ius = iUS; igw = iGW;
        area = Ele[i].area;
        QeleSurfTot[i] = Qe2r_Surf[i];
        QeleSubTot[i] = Qe2r_Sub[i];
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurf[i][j];
            QeleSubTot[i] += QeleSub[i][j];
            CheckNANij(QeleSurf[i][j], i, "QeleSurf[i][j]");
            CheckNANij(QeleSub[i][j], i, "QeleSub[i][j]");
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area - qEs[i];
        DY[ius] = qEleInfil[i] - qEleRecharge[i] - qEu[i] - qTu[i];
        DY[igw] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area - qEg[i] - qTg[i];
//        if(uYgw[i] > Ele[i].WetlandLevel ){ /* IF GW above surface, exfil from GW, OR, from UNSAT*/
//            DY[igw] += - qEleExfil[i];
//        }else{
//            if(uYus[i] > 0.){
//                DY[ius] += - qEleExfil[i];
//            }else{
//                DY[igw] += - qEleExfil[i];
//            }
//        }
        /* Boundary condition and Source/Sink */
        if(Ele[i].iBC == 0){
        }else if(Ele[i].iBC > 0){ // Fix head of GW.
            DY[igw] = 0;
        }else if(Ele[i].iBC < 0){ // Fix flux in GW
            DY[igw] += Ele[i].QBC / area;
        }

        if(Ele[i].iSS == 0){
        }else if(Ele[i].iSS > 0){ // SS in Landusrface
            DY[isf] += Ele[i].QSS / area;
        }else if(Ele[i].iSS < 0){ // SS in GW
            DY[igw] += Ele[i].QSS / area;
        }
        /* Convert with specific yield */
        DY[ius] /= Ele[i].Sy;
        DY[igw] /= Ele[i].Sy;

//        if(i+1==ID_ELE && DY[i]+uYsf[i] > 0.1 ){  // debug only
//            printf("%.3f, %d: %.2e, %.2e, %.2e | (%.2e, %.2e, %.2e, %.2e, %.2e )\n",
//                   t, i+1,
//                   uYsf[i], uYus[i], uYgw[i], qEleInfil[i], - qEleRecharge[i], - qEu[i],  - qTu[i], QeleSurfTot[i] / area);
//            printf("\n");
//        }
//        if(i +1== ID_ELE){  // debug only
//            printf("%.3f, %d: %f, %f, %f | (%f, %f, %f), %.2e, %.2e\n",
//                   t, i+1,
//                   DY[i], DY[ius], DY[igw], QeleSurf[i][0] / area, QeleSurf[i][1] / area, QeleSurf[i][2] / area,
//                   -QeleSurfTot[i] / area, qEleRecharge[i]);
//            printf("\n");
//        }
        if(Ele[i].iLake > 0){
            DY[i] = 0.;
            DY[ius] = 0.;
            DY[igw] = 0.;
        }
#ifdef DEBUG
        CheckNANi(DY[i], i, "DY[i] (Model_Data::f_applyDY)");
        CheckNANi(DY[ius], i, "DY[ius] (Model_Data::f_applyDY)");
        CheckNANi(DY[igw], i, "DY[igw] (Model_Data::f_applyDY)");
#endif
    }
    for (int i = 0; i < NumRiv; i++) {
        if(Riv[i].BC > 0){
//            Newmann condition.
            DY[iRIV] = 0.;
        }else{
            DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].Length; // dA on CS
            if(DY[iRIV] < -1. * Riv[i].u_CSarea){ /* The negative dA cannot larger then Availalbe Area. */
                DY[iRIV] = -1. * Riv[i].u_CSarea;
            }
            DY[iRIV] = fun_dAtodY(DY[iRIV], Riv[i].u_topWidth, Riv[i].bankslope);

//            if(i+1 == ID_RIV && fabs(DY[iRIV])> 0.001){
//                printf("%d, up:%f, down:%f, sub:%f, surf:%f, dy=%f\n", i+1,
//                       - QrivUp[i] / Riv[i].u_TopArea, - QrivDown[i] / Riv[i].u_TopArea,
//                       - QrivSub[i] / Riv[i].u_TopArea, - QrivSurf[i] / Riv[i].u_TopArea,
//                       DY[iRIV]);
//                i=i;
//            }
        }
#ifdef DEBUG
        CheckNANi(DY[i + 3 * NumEle], i, "DY[i] of river (Model_Data::f_applyDY)");
#endif
    }
    for(int i = 0; i < NumLake; i++){
//        DY[i + 3 * NumEle + NumRiv]
        DY[iLAKE] = qLakePrcp[i] - qLakeEvap[i]  +
                    (QLakeRivIn[i] - QLakeRivOut[i] + QLakeSub[i] + QLakeSurf[i] ) / y2LakeArea[i] ;
//        if(fabs(DY[iLAKE]) > 1.0e-4){
//            printf("%f: %g + %g\n",t, yLakeStg[i], DY[iLAKE]);
//            i=i;
//        }
#ifdef DEBUG
        CheckNANi(DY[iLAKE], i, "DY[i] of LAKE (Model_Data::f_applyDY)");
#endif
    }
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_applyDY", t, DY, 3 * NumEle + NumRiv + NumLake);
#endif
}

/* S1d.1 (openMP #47) — `rhs_core` four-arg ExecPolicy dispatch.
 *
 * Per spec exec-policy-enum + design.md D7 / D8:
 *   - `ExecPolicy::Serial` runs the full new-path chain
 *     (`rhs_update` -> `rhs_flux` -> `rhs_apply`) extracted in
 *     S1a/b/c — bitwise-identical to legacy `f_update/f_loop/f_applyDY`
 *     (verified vs B0-tag in tasks 4.6a/b).
 *   - `ExecPolicy::StrictOMP` and `ExecPolicy::ProductionOMP` are
 *     S1-phase compile-time stubs: each calls `std::abort()` so any
 *     runtime call SIGABRTs immediately. `assert(false)` is forbidden
 *     because `-DNDEBUG` (release / EXTRA_CXXFLAGS=-DNDEBUG smoke
 *     compile) strips assert to a no-op and would let execution fall
 *     through silently to the next statement — destroying the
 *     contract that an OMP-policy call cannot impersonate Serial.
 *   - SHUD_ENABLE_OPENMP_RHS=0 (default) `#ifdef`s the OMP cases
 *     out of the translation unit; the resulting binary contains no
 *     OMP-path symbols (verified in tasks 4.9 / 5.10c via `nm`).
 *   - `default:` branch also aborts to catch ABI drift / future
 *     enumerator additions that haven't been wired up here.
 *
 * No template specialization, no virtual dispatch — plain switch.
 * The switch is compile-time-known at every caller in this stage
 * (f.cpp always passes `ExecPolicy::Serial`), so the compiler
 * eliminates the branch in optimized builds.
 */
void Model_Data::rhs_core(double *Y, double *DY, double t, ExecPolicy policy){
    switch (policy) {
        case ExecPolicy::Serial:
            rhs_update(Y, DY, t);
            rhs_flux(t);
            rhs_apply(DY, t);
            break;
#ifdef SHUD_ENABLE_OPENMP_RHS
        case ExecPolicy::StrictOMP:
            /* S2+ scope. S1 stub: abort to prevent silent fall-through.
             * NOT assert(false) — -DNDEBUG strips it. */
            std::abort();
        case ExecPolicy::ProductionOMP:
            std::abort();
#endif
        default:
            /* Catches enumerator additions not yet wired + ABI drift. */
            std::abort();
    }
}
