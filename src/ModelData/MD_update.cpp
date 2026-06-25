#include "Model_Data.hpp"
#include <vector>  /* P1e PR-B0 (#323): recompute_for_output scratch DY */
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif

void Model_Data::f_updatei(double  *Y, double *DY, double t, int flag){
    switch (flag) {
        case 1:
            for (int i = 0; i < NumEle; i++) {
                uYsf[i] = (Y[i] >= 0.) ? Y[i] : 0.;
            }
            break;
        case 2:
            for (int i = 0; i < NumEle; i++) {
                uYus[i] = (Y[i] >= 0.) ? Y[i] : 0.;
            }
            break;
        case 3:
            for (int i = 0; i < NumEle; i++) {
                uYgw[i] = (Y[i] >= 0.) ? Y[i] : 0.;
                if(Ele[i].iBC == 0){ // NO BC
                    uYgw[i] = max(0.0, Y[i]);
                    Ele[i].QBC = 0.;
                }else if(Ele[i].iBC > 0){ // BC fix head
                    Ele[i].yBC = tsd_eyBC.getX(t, Ele[i].iBC);
                    uYgw[i] = Ele[i].yBC;
                    Ele[i].QBC = 0.;
                }else{ // BC fix flux to GW
                    Ele[i].QBC = tsd_eqBC.getX(t, -Ele[i].iBC);
                }
            }
            break;
        case 4:
            for (int i = 0; i < NumRiv; i++) {
                uYriv[i] = (Y[i] >= 0.) ? Y[i] : 0.;
//                uYriv[i] = Y[i];
//                QrivSurf[i] = 0.;
//                QrivSub[i] = 0.;
                QrivUp[i] = 0.;
                QrivDown[i] = 0.;
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
            }
            break;
        case 5:
            for (int i = 0; i < NumLake; i++) {
                uYlake[i] = (Y[i] >= 0.) ? Y[i] : 0.;
            }
            break;
        default:
            break;
    }
}
/* P1d.2.0 PR-C0 (#291): Model_Data f_update legacy carry-over deleted.
 * Live counterpart is Model_Data::rhs_update in MD_rhs_core.cpp
 * (~L58-149), which CVODE invokes via f.cpp:54 -> MD->rhs_core(...).
 * The "f_update" SHUD_DUMP_RHS tag string at MD_rhs_core.cpp:147 +
 * MD_rhs_dump.{h,cpp} default site name are preserved as the
 * golden-file dump contract. */
void Model_Data::summary (N_Vector udata){
    double  *Y;
    /* S1d.2 (openMP #48) — generic N_Vector data accessor (see f.cpp
     * + Macros.hpp comments). Replaces the prior backend-specific
     * NV_DATA_OMP / NV_DATA_S branch. */
    Y = N_VGetArrayPointer(udata);
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y[iSF];
        yEleUnsat[i] = Y[iUS];
        
        if(Ele[i].iBC > 0){
            yEleGW[i] = Ele[i].yBC;
        }else{
            yEleGW[i] = Y[iGW];
        }
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y[iRIV];
        //        uYriv[i] = Y[iRIV];
        if(Riv[i].BC > 0){
            yRivStg[i] = Riv[i].yBC;
        }else{
            yRivStg[i] = Y[iRIV];
        }
    }
}
/* P1e PR-B0 (#323): recompute flux/gather caches from Y(tout) before
 * ExportResults so PCtrl-aliased output buffers (QrivDown, QrivUp,
 * QrivSurf, QrivSub, QLakeRivIn, QLakeRivOut, QLakeSurf, QLakeSub,
 * qLakeEvap, qLakePrcp, Qe2r_Surf, Qe2r_Sub, qEle*, QeleSubTot,
 * QeleSurfTot, etc.) reflect deterministic Y(tout)-derived state, NOT
 * the side-effect cache left by the last internal-step f() at
 * t_internal != tout under CV_NORMAL mode (per
 * docs/p1e/p1e_rivqdown_cache_audit.md conclusion + spec
 * p1e-strict-omp-rhs L260-285 + design D5 option 1).
 *
 * Mechanics: re-run the full RHS chain `rhs_update -> rhs_flux ->
 * rhs_apply` exactly once at (Y, t)=(udata, t) using a local scratch
 * DY buffer. rhs_apply IS idempotent at fixed (Y, t) (verified by
 * Phase 4.5 verifier — MD_rhs_core.cpp L647-649 leading `=` resets +
 * inner j∈[0,3) `+=` accumulates per-i), so calling it is safe; the
 * DY_scratch mutation is harmless because the scratch buffer is
 * discarded on return.
 *
 * Why rhs_apply MUST be called (Phase 6 fix per outer #323): without
 * this call QeleSubTot[i] and QeleSurfTot[i] stay at the values
 * written by rhs_update (set to zero at MD_rhs_core.cpp:91 / :104) —
 * PrintData tau-averaging would then emit silent all-zero data into
 * any PCtrl-aliased *.eleQsubTot.dat / *.eleQsurfTot.dat output when
 * DT_QE_SUB > 0 or DT_QE_SURF > 0. The earlier "skip rhs_apply"
 * rationale was based on a non-idempotent-+= misread; Phase 4.5
 * verifier refuted that. See docs/p1e/p1e_pr_b0_rivqdown_recompute.md
 * §"rhs_apply Phase 6 fix rationale".
 *
 * Sibling-cache coverage: rhs_update + rhs_flux + rhs_apply is the
 * ground-truth RHS chain that originally populates every river / lake
 * / element cache exposed via PCtrl::Init() (MD_initialize.cpp
 * L307-391) + flood->InitPointer alias (Model_Data.cpp:427). Reusing
 * the existing chain (instead of a hand-rolled per-channel recompute)
 * guarantees all sibling caches are recomputed in the same pass, with
 * byte-equivalent iteration / floating-point operation order.
 *
 * Side effects extend beyond PCtrl-aliased output buffers to all
 * RHS-touched scratch state (Ele[i].{QBC, u_effKH, ...},
 * Riv[i].{u_Ystage, u_CSarea, ...}, lake[i].{u_toparea, ...},
 * hot.*[i]); these are deterministic functions of (Y, t) and are
 * overwritten on the next f() call, so the leakage is benign.
 *
 * Under SHUD_ENABLE_DIAGNOSTICS builds, the 5 shud_diag::ScopeTimer
 * instrumentations inside rhs_flux (MD_rhs_core.cpp L365 / L407 /
 * L423 / L433 / L496) will see double-counted bucket entries on this
 * extra call; subtract one outer-tick bucket per ScopeTimer in
 * diagnostics post-processing if exact attribution is needed.
 * Default builds (no DIAGNOSTICS macro) unaffected.
 *
 * Counter discipline: `nFCall` is NOT incremented. This call is a
 * tout-boundary cache refresh for output, not a CVODE-driven RHS
 * evaluation; nFCall semantics (per Model_Data.hpp L58 + S5c-C #175)
 * remain the count of solver-internal f() invocations.
 *
 * Profile / diagnostics: deliberately not wrapped in shud_profile or
 * shud_diag timers — recompute time is attributed to t_other (not
 * t_RHS_total) so existing profile decomposition stays interpretable.
 *
 * Determinism: `rhs_update + rhs_flux + rhs_apply` are deterministic
 * functions of (Y, t, time-series state, model config). Same (Y, t)
 * -> same caches.
 *
 * Scope: this call site is wired in only on the coupled-mode MainLoop
 * (SHUD(), shud.cpp:203). Uncoupled mode (SHUD_uncouple(),
 * shud.cpp:412-413 under `-g` CLI flag) runs 5 split CVode integrations
 * each with its own internal-step cache and does NOT receive this
 * helper call; that scope is deferred to a future issue. Workaround
 * for uncoupled users: do not use `-g` for runs that require bitwise
 * reproducibility — default coupled mode is unaffected.
 */
void Model_Data::recompute_for_output(N_Vector udata, double t){
    double *Y = N_VGetArrayPointer(udata);
    /* Scratch DY consumed by rhs_update (zeroes it L218-220) +
     * rhs_apply (writes derivative components L657-659/etc.); ignored
     * on rhs_flux input. NumY = 3*NumEle + NumRiv + NumLake, matches
     * the solver state vector layout. */
    std::vector<double> DY_scratch(NumY, 0.0);
    rhs_update(Y, DY_scratch.data(), t);
    rhs_flux(t);
    /* PR-B0 Phase 6 fix (cand-2): rhs_apply IS idempotent at fixed
     * (Y, t). Without this call, QeleSubTot/QeleSurfTot stay at the
     * rhs_update zero-out (MD_rhs_core.cpp:91 / :104) → PCtrl-aliased
     * *.eleQsubTot.dat / *.eleQsurfTot.dat would silently emit
     * all-zero data when DT_QE_SUB>0 or DT_QE_SURF>0. The DY_scratch
     * mutation (rhs_apply writes derivative components) is harmless —
     * scratch is discarded on return. */
    rhs_apply(DY_scratch.data(), t);
}
void Model_Data::summary (N_Vector u1, N_Vector u2, N_Vector u3, N_Vector u4, N_Vector u5){

    /* S1d.2 (openMP #48) — generic N_VGetArrayPointer collapse of the
     * prior backend-split blocks (NV_Ith_OMP / NV_Ith_S). Hoist the
     * five pointer fetches out of the inner loops so the per-element
     * indexing stays a single load-and-store; bitwise-equivalent to
     * `NV_Ith_S(v, i)` because nvector_serial's NV_Ith_S expands to
     * `((NV_DATA_S(v))[i])` (verified against sundials 6.0.0). */
    double *Y1 = N_VGetArrayPointer(u1);
    double *Y2 = N_VGetArrayPointer(u2);
    double *Y3 = N_VGetArrayPointer(u3);
    double *Y4 = N_VGetArrayPointer(u4);
    double *Y5 = N_VGetArrayPointer(u5);
    for (int i = 0; i < NumEle; i++){
        yEleSurf[i] = Y1[i];
        yEleUnsat[i] = Y2[i];
        yEleGW[i] = Y3[i];
        if(Ele[i].iBC > 0){
            yEleGW[i] = Ele[i].yBC;
        }
    }
    for (int i = 0; i < NumRiv; i++){
        yRivStg[i] = Y4[i];
        if(Riv[i].BC > 0){
            yRivStg[i] = Riv[i].yBC;
        }
    }
    for (int i = 0; i < NumLake; i++){
        yLakeStg[i] = Y5[i];
    }
    Sub2Global(yEleSurf, yEleUnsat, yEleGW, yRivStg, yLakeStg, NumEle, NumRiv, NumLake);
//    printVector(stdout, yEleSurf, 0, NumEle, 0);
//    printVector(stdout, yEleUnsat, 0, NumEle, 0);
//    printVector(stdout, yEleGW, 0, NumEle, 0);
//    printVector(stdout, yRivStg, 0, NumRiv, 0);
//
//    printVector(stdout, globalY, 0, NumEle, 0);
//    printVector(stdout, globalY, NumEle, NumEle, 0);
//    printVector(stdout, globalY, NumEle*2, NumEle, 0);
//    printVector(stdout, globalY, NumEle*3, NumRiv, 0);
}

int Model_Data::PrintInit (const char *fn, double t){
    unsigned long t_long = (long) t;
    if( t_long % CS.UpdateICStep ){
        return 0;
    }
    FILE           *fp;
    fp = fopen (fn, "w");
    CheckFile(fp, fn);
    /************* Element status **************/
    fprintf (fp, "%d\t %d \t%lf\n", NumEle, 6, t);
    fprintf (fp, "%s\t%s\t%s\t%s\t%s\t%s\n","Index",
             "Canopy", "Snow", "Surface", "Unsat", "GW");
    for (int i = 0; i < NumEle; i++){
        fprintf (fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", i+1, yEleIS[i], yEleSnow[i], yEleSurf[i], yEleUnsat[i], yEleGW[i]);
    }
    /************* River Reach status **************/
    fprintf (fp, "%d\t%d\n", NumRiv, 2);
    fprintf (fp, "%s\t%s\n", "Index", "Stage");
    for (int i = 0; i < NumRiv; i++){
        fprintf (fp, "%d\t%lf\n", i+1, yRivStg[i]);
    }
    /************* Lake status **************/
    if(NumLake > 0){
        fprintf (fp, "%d\t%d\n", NumLake, 2);
        fprintf (fp, "%s\t%s\n", "Index", "LakeStage");
        for (int i = 0; i < NumLake; i++){
            fprintf (fp, "%d\t%lf\n", i+1,yLakeStg[i]);
        }
    }
    fclose (fp);
    return 1;
}
