#include "Model_Data.hpp"
#include <vector>  /* recompute_for_output scratch DY */
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
/* The coupled-mode state update is Model_Data::rhs_update in
 * MD_rhs_core.cpp, which CVODE reaches through MD->rhs_core(...).
 * Its SHUD_DUMP_RHS site name is "f_update", which is also the
 * default site in MD_rhs_dump.{h,cpp}; snapshot files are selected by
 * that name, so it must not be renamed. */
void Model_Data::summary (N_Vector udata){
    double  *Y;
    /* Generic N_Vector data accessor (see f.cpp + Macros.hpp
     * comments): works for both the serial and the OpenMP N_Vector
     * backend, so no NV_DATA_OMP / NV_DATA_S branch is needed. */
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
/* Recompute the flux/gather caches from Y(tout) before ExportResults
 * so the output buffers registered with PCtrl (QrivDown, QrivUp,
 * QrivSurf, QrivSub, QLakeRivIn, QLakeRivOut, QLakeSurf, QLakeSub,
 * qLakeEvap, qLakePrcp, Qe2r_Surf, Qe2r_Sub, qEle*, QeleSubTot,
 * QeleSurfTot, etc.) hold values derived from Y(tout), NOT whatever
 * the last internal-step f() left behind at t_internal != tout under
 * CV_NORMAL mode.
 *
 * Mechanics: re-run the full RHS chain `rhs_update -> rhs_flux ->
 * rhs_apply` exactly once at (Y, t)=(udata, t) using a local scratch
 * DY buffer. rhs_apply is idempotent at fixed (Y, t) (each per-element
 * total is reset with `=` before its three edges are added with `+=`),
 * so calling it is safe; the DY_scratch mutation is harmless because
 * the scratch buffer is discarded on return.
 *
 * Why rhs_apply MUST be called: rhs_update zeroes QeleSubTot[i] and
 * QeleSurfTot[i] and only rhs_apply refills them. Without it the
 * *.eleQsubTot.dat / *.eleQsurfTot.dat outputs (enabled when
 * DT_QE_SUB > 0 or DT_QE_SURF > 0) would silently be all zero.
 *
 * Reusing the RHS chain (instead of a hand-written per-array
 * recompute) refreshes every river / lake / element cache exposed
 * through PCtrl::Init() and flood->InitPointer in one pass, with the
 * same iteration and floating-point operation order as the solver.
 *
 * Side effects extend beyond the output buffers to all RHS-touched
 * scratch state (Ele[i].{QBC, u_effKH, ...},
 * Riv[i].{u_Ystage, u_CSarea, ...}, lake[i].{u_toparea, ...},
 * hot.*[i]); these are deterministic functions of (Y, t) and are
 * overwritten on the next f() call, so this is benign.
 *
 * Under SHUD_ENABLE_DIAGNOSTICS builds, the shud_diag::ScopeTimer
 * buckets inside rhs_flux also count this extra call; subtract one
 * call per output step in post-processing if exact attribution is
 * needed. Default builds are unaffected.
 *
 * `nFCall` is NOT incremented: this is a cache refresh for output,
 * not a CVODE-driven RHS evaluation, and nFCall counts only
 * solver-internal f() invocations.
 *
 * Deliberately not wrapped in shud_profile or shud_diag timers —
 * recompute time is attributed to t_other (not t_RHS_total) so the
 * profile decomposition stays interpretable.
 *
 * Scope: called only from the coupled-mode main loop (SHUD() in
 * shud.cpp). Uncoupled mode (SHUD_uncouple(), `-g` CLI flag) runs 5
 * split CVode integrations, each with its own internal-step cache,
 * and does NOT call this helper; do not use `-g` for runs that
 * require bitwise-reproducible output.
 */
void Model_Data::recompute_for_output(N_Vector udata, double t){
    double *Y = N_VGetArrayPointer(udata);
    /* Scratch DY consumed by rhs_update (zeroes it) + rhs_apply
     * (writes derivative components); not used by rhs_flux.
     * NumY = 3*NumEle + NumRiv + NumLake, matches the solver state
     * vector layout. */
    std::vector<double> DY_scratch(NumY, 0.0);
    rhs_update(Y, DY_scratch.data(), t);
    rhs_flux(t);
    /* Required: refills QeleSubTot/QeleSurfTot, which rhs_update
     * zeroed (see the function comment above). */
    rhs_apply(DY_scratch.data(), t);
}
void Model_Data::summary (N_Vector u1, N_Vector u2, N_Vector u3, N_Vector u4, N_Vector u5){

    /* N_VGetArrayPointer works for both the serial and the OpenMP
     * N_Vector backend, so no NV_Ith_OMP / NV_Ith_S split is needed.
     * The five pointer fetches are hoisted out of the loops so the
     * per-element indexing stays a single load-and-store; this is
     * bitwise-equivalent to `NV_Ith_S(v, i)`, which expands to
     * `((NV_DATA_S(v))[i])` (sundials 6.0.0). */
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
