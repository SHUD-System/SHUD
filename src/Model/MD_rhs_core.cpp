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
                              * rhs_core (PR-11 / S3c will replace PassValue_legacy's
                              * in-loop gather). The include here satisfies
                              * spec s4-adjacency-topology Scenario
                              * "MD_rhs_core.cpp 文件顶部 SHALL 包含
                              * #include MD_adjacency.hpp". */
#include <cstdlib>   /* std::abort -- S1d.1 OMP stub regression guard */
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif
/* S5c-B (#174): RHS 7-bucket diagnostic timer. Header is empty under
 * default (SHUD_ENABLE_DIAGNOSTICS undefined) — zero new code in hot
 * path, zero floating-point change, B1a-tag bitwise contract intact. */
#include "MD_diagnostics.hpp"

#ifdef SHUD_ENABLE_DIAGNOSTICS
namespace shud_diag {
/* Definition of the 7-bucket nanosecond accumulators declared in
 * MD_diagnostics.hpp. Single TU storage — single-threaded driver makes
 * this race-free under B1a contract. */
long long g_rhs_timer_ns[RHS_BUCKET_COUNT] = {0, 0, 0, 0, 0, 0, 0};
}  // namespace shud_diag
#endif

/* P1d.2.1 (#277) — NUMA first-touch gate, defined in shud.cpp L49 and
 * set once by emit_numa_token() at SHUD() entry (L70-91). Read here to
 * gate the steady-state element first-touch warm-up below (mirrors the
 * extern usage in Model_Data.cpp L10 / MD_initialize.cpp L27). */
extern int g_numa_first_touch_enabled;

void Model_Data::rhs_update(double *Y, double *DY, double t){
    /* P1d.2.1 (#277) — steady-state element first-touch warm-up.
     * Mirrors allocation-time flat3 zero-write in Model_Data.cpp::
     * malloc_EleRiv L302-317; complements (does NOT replace) it by
     * re-touching the element-owned hot flux/scratch arrays on each
     * CVODE RHS evaluation so each owner thread page-faults its own
     * NUMA-local pages every step. Pure zero-write, no reduction /
     * accumulation / non-zero write (design D2 + R1). Every field
     * below is re-zeroed or fully overwritten LATER in this same
     * rhs_update body before any read:
     *   QeleSubAt(i,j)  -> re-zeroed L70   (`QeleSubAt(i,j)  = 0.`)
     *   QeleSurfAt(i,j) -> re-zeroed L71   (`QeleSurfAt(i,j) = 0.`)
     *   QeleSubTot[i]   -> re-zeroed L72   (`QeleSubTot[i]   = 0.`)
     *   QeleSurfTot[i]  -> re-zeroed L73   (`QeleSurfTot[i]  = 0.`)
     *   Qe2r_Surf[i]    -> re-zeroed L134  (`Qe2r_Surf[i]    = 0.`)
     *   Qe2r_Sub[i]     -> re-zeroed L135  (`Qe2r_Sub[i]     = 0.`)
     * Field set = element-owned subset only (no river / lake fields)
     * per OQ1 doc docs/p1d/p1d_first_touch_design.md. Gated by
     * g_numa_first_touch_enabled to match allocation-time gate
     * (preserves deterministic-vs-serial behavior when OMP_PROC_BIND
     * is unset; see shud.cpp L41-L77). */
    if (g_numa_first_touch_enabled) {
        int i;
#pragma omp parallel for schedule(static) default(none) shared(QeleSurfTot, QeleSubTot, Qe2r_Surf, Qe2r_Sub) private(i)
        for (i = 0; i < NumEle; i++) {
            for (int j = 0; j < 3; j++) {
                QeleSurfAt(i, j) = 0.0;
                QeleSubAt(i, j)  = 0.0;
            }
            QeleSurfTot[i] = 0.0;
            QeleSubTot[i]  = 0.0;
            Qe2r_Surf[i]   = 0.0;
            Qe2r_Sub[i]    = 0.0;
        }
    }

    for (int i = 0; i < NumEle; i++) {
//        uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
//        uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
        /* S5d.2-5a (#179) — flat zero via accessor. */
        for(int j = 0; j < 3; j++){
            QeleSubAt(i, j) = 0.;
            QeleSurfAt(i, j) = 0.;
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
    /* P1d.2.3 (#279) — steady-state lake first-touch warm-up.
     * Mirrors allocation-time first-touch + PR-C rhs_update element
     * first-touch (L65-98) + PR-D rhs_flux river first-touch (L327-358);
     * complements (does NOT replace) the allocation-time zero-write by
     * re-touching the lake-owned hot flux/scratch arrays on each CVODE
     * RHS evaluation so each owner thread page-faults its own NUMA-local
     * pages every step. Pure zero-write, no reduction / accumulation /
     * non-zero write (design D2 + R1). Every field below is re-zeroed
     * LATER in this same rhs_update lake block (the for-i NumLake loop
     * just below) before any read in rhs_flux / rhs_apply:
     *   QLakeSub[i]    -> re-zeroed L214 (`QLakeSub[i]    = 0.`)
     *   QLakeSurf[i]   -> re-zeroed L215 (`QLakeSurf[i]   = 0.`)
     *   qLakeEvap[i]   -> re-zeroed L216 (`qLakeEvap[i]   = 0.`)
     *   qLakePrcp[i]   -> re-zeroed L217 (`qLakePrcp[i]   = 0.`)
     *   QLakeRivIn[i]  -> re-zeroed L218 (`QLakeRivIn[i]  = 0.`)
     *   QLakeRivOut[i] -> re-zeroed L219 (`QLakeRivOut[i] = 0.`)
     * Field set = lake-owned subset only (no element / river fields),
     * per OQ1 doc docs/p1d/p1d_first_touch_design.md §"LAKE-owned 子集".
     * Excluded: yLakeStg[i] (written = Y[iLAKE] non-zero at L210) and
     * y2LakeArea[i] (written = lake[i].u_toparea non-zero at L213) —
     * both are persistent state writers, not pure-zero candidates.
     * Gated by g_numa_first_touch_enabled to match allocation-time +
     * PR-C + PR-D gate (preserves deterministic-vs-serial behavior when
     * OMP_PROC_BIND is unset; see shud.cpp L41-L77). */
    if (g_numa_first_touch_enabled) {
        int i;
#pragma omp parallel for schedule(static) default(none) shared(QLakeSub, QLakeSurf, qLakeEvap, qLakePrcp, QLakeRivIn, QLakeRivOut) private(i)
        for (i = 0; i < NumLake; i++) {
            QLakeSub[i]    = 0.0;
            QLakeSurf[i]   = 0.0;
            qLakeEvap[i]   = 0.0;
            qLakePrcp[i]   = 0.0;
            QLakeRivIn[i]  = 0.0;
            QLakeRivOut[i] = 0.0;
        }
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
 * includes both #ifdef SHUD_DUMP_RHS probe blocks + after-PassValue_legacy
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
 *   6. PassValue_legacy()     (zero-reset + segment re-accumulate)
 * Followed by the after-PassValue_legacy no-op SHUD_DUMP_RHS hook (legacy).
 *
 * Dump tag strings `"f_loop_before_passvalue"` (MD_f.cpp:67) and
 * `"f_loop"` (MD_f.cpp:72) are preserved verbatim so PR #54's 12
 * `_before_passvalue.bin` + 12 unsuffixed snapshot goldens (4 case
 * × 3 t_values) stay addressable. Renaming = breaks goldens.
 *
 * No sub-function split (rhs_element_vertical / rhs_segment_compute
 * etc. belong to S3, not S1b — per spec Scenario "No sub-function
 * split inside rhs_flux"). */

/* P1c PR-B (#245): Fixed-shape pairwise tree reduction over an index
 * list, B0 serial traversal order preserved. Tree shape depends only on
 * the list length, NOT on NUM_OPENMP — used at per-owner accumulation
 * sites to obtain a thread-count-independent canonical sum. Spec
 * p1c-deterministic-reduction "fixed-shape pairwise canonical reduction".
 *
 * Empty list returns 0.0 (preserves prior explicit zero-init semantics);
 * single element returns src[idx[0]]. Range-pointer variant avoids
 * per-call vector copies (O(n log n) work, O(log n) stack). Stack depth
 * = ceil(log2(n)); n is bounded by NumEle so recursion is safe. */
static inline double fixed_pairwise_sum_range(
        const int* idx, std::size_t n, const double* src) {
    if (n == 0) return 0.0;
    if (n == 1) return src[idx[0]];
    if (n == 2) return src[idx[0]] + src[idx[1]];
    const std::size_t mid = n / 2;
    const double lo = fixed_pairwise_sum_range(idx, mid, src);
    const double hi = fixed_pairwise_sum_range(idx + mid, n - mid, src);
    return lo + hi;
}

static inline double fixed_pairwise_sum_indexed(
        const std::vector<int>& idx, const double* src) {
    return fixed_pairwise_sum_range(idx.data(), idx.size(), src);
}

/* P1c PR-C (#246): Fixed-shape leftfold canonical reduction over an
 * index list (B0 traversal order preserved). For each src element
 * src[idx[k]] in list order, accumulator runs
 *   `acc = 0; acc += src[idx[0]]; acc += src[idx[1]]; ...`.
 * Result is bitwise-identical to legacy serial `for x in list: dst +=
 * src[x]` accumulation (assuming dst starts at 0.0) since the operation
 * order is identical to a left-fold pattern. Use this for sites 3-7
 * (segment -> river / element gathers, upstream river gather) where
 * keliya N=1 bitwise vs B0 must be preserved. Sites 1-2 (lake
 * aggregation) use the tree variant `fixed_pairwise_sum_indexed`
 * instead (different shape trade-off documented in PR-B). */
static inline double fixed_leftfold_sum_indexed(
        const std::vector<int>& idx, const double* src) {
    double acc = 0.0;
    for (int i : idx) acc += src[i];
    return acc;
}

/* Fixed-shape leftfold canonical reduction over an index-pair list
 * (B0 traversal order preserved). For each (ie, j) pair in B0
 * traversal order, accumulator runs `acc += src[ie * stride + j]`.
 * Bitwise-identical to legacy `for (auto& ej : list) acc += src[ej.first
 * * stride + ej.second]`. Used for site 8 (lake_bank_edge_by_lake,
 * S4.6) where adjacency list contains (ie, j) pairs; lookup is the
 * standard stride=3 element-edge convention. Constant-stride keeps the
 * helper trivially inlinable.
 */
static inline double fixed_leftfold_sum_pair_indexed(
        const std::vector<std::pair<int,int>>& pairs,
        const double* src,
        int stride) {
    double acc = 0.0;
    for (const auto& ej : pairs) acc += src[ej.first * stride + ej.second];
    return acc;
}

void Model_Data:: rhs_flux(double t){
    int i;
    /* P1d.2.2 (#278) — steady-state river first-touch warm-up.
     * Mirrors allocation-time first-touch in Model_Data.cpp::malloc_EleRiv
     * + PR-C rhs_update element first-touch (MD_rhs_core.cpp L65-98);
     * complements (does NOT replace) the allocation-time zero-write by
     * re-touching the river-owned hot flux arrays on each CVODE RHS
     * evaluation so each owner thread page-faults its own NUMA-local
     * pages every step. Pure zero-write, no reduction / accumulation /
     * non-zero write (design D2 + R1). Every field below is re-zeroed
     * and then overwritten LATER in this same rhs_flux body (via the
     * rhs_deterministic_gather() call at L502) before any read:
     *   QrivSurf[i] -> re-zeroed L545 (`QrivSurf[i] = 0.` in gather pre-zero)
     *                  -> overwritten L561 (`QrivSurf[ir] = fixed_leftfold_sum_indexed(...)`)
     *   QrivSub[i]  -> re-zeroed L546 (`QrivSub[i] = 0.`)
     *                  -> overwritten L562 (`QrivSub[ir]  = fixed_leftfold_sum_indexed(...)`)
     *   QrivUp[i]   -> re-zeroed L547 (`QrivUp[i] = 0.`)
     *                  -> overwritten L582 (`QrivUp[ir]  = -fixed_leftfold_sum_indexed(...)`)
     * All 3 fields are first READ in rhs_apply L719
     * (`DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + ...)`),
     * which runs AFTER rhs_flux returns. Field set = river-owned subset
     * only (no element / lake fields) per OQ1 doc
     * docs/p1d/p1d_first_touch_design.md §"RIVER-owned 子集". Gated by
     * g_numa_first_touch_enabled to match allocation-time + PR-C gate
     * (preserves deterministic-vs-serial behavior when OMP_PROC_BIND is
     * unset; see shud.cpp L41-L77). */
    if (g_numa_first_touch_enabled) {
#pragma omp parallel for schedule(static) default(none) shared(QrivSurf, QrivSub, QrivUp) private(i)
        for (i = 0; i < NumRiv; i++) {
            QrivSurf[i] = 0.0;
            QrivSub[i]  = 0.0;
            QrivUp[i]   = 0.0;
        }
    }

    /* S5c-B (#174): 5 inner buckets (ET / lateral / segment / river /
     * gather) are scoped via `shud_diag::ScopeTimer`. The block braces
     * are pre-existing in some loops (none here) so we add explicit
     * `{ }` to scope each ScopeTimer to exactly one phase. Under
     * `SHUD_ENABLE_DIAGNOSTICS` undefined, the timer macros expand to
     * nothing — the loop bodies and brace nesting are unchanged. */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_et(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_ET]);
#endif
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            Ele[i].updateLakeElement();
            /* S6b.4 (#205): SoA/AoS sync drift fix. updateLakeElement()
             * mutates AoS hot fields (Ele[i].u_effKH = KsatH, etc.) but
             * does NOT refresh hot.<field>[i]; sync_hot_dynamic(i) is
             * required so subsequent consumers reading hot.u_effKH[i]
             * see the post-update value rather than the stale forcing
             * blend left by updateforcing(). Pattern mirrors the
             * dead-code legacy MD_f.cpp::f_loop L27-28. */
            sync_hot_dynamic(i);
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
            /* S6b.4 (#205): same SoA refresh pattern after updateElement
             * mutation. Mirrors MD_f.cpp::f_loop L40-41 dead-code
             * pattern; Flux_Infiltration / Flux_Recharge consumers read
             * hot.u_effkInfi[i] etc. and must see post-update value. */
            sync_hot_dynamic(i);
            fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
            fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
        }
    }
    } /* end ET bucket */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_lat(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_LATERAL]);
#endif
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
    } /* end lateral bucket */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_seg(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_SEGMENT]);
#endif
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    } /* end segment bucket */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_riv(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_RIVER]);
#endif
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    /* S3b.4 (PR-9): transitional gather per-element -> per-lake.
     * Must run BEFORE the lake clamp below (clamp reads qLakeEvap /
     * qLakePrcp). Cannot live in PassValue_legacy because PassValue_legacy() is
     * called AFTER the clamp. Will be replaced by
     * rhs_deterministic_gather() in S3c (PR-11).
     *
     * P1c PR-B (#245): L278/L279 inline serial += replaced with per-lake
     * fixed-shape pairwise tree reduction over `ele_by_lake[ilake]`
     * (S4.5; populated in MD_adjacency.cpp by the same
     * `for (i = 0; i < NumEle; i++) if Ele[i].iLake > 0` traversal so the
     * list is already in B0 canonical order — do NOT re-sort). Tree
     * shape determined by list length only, NOT NUM_OPENMP. Prior
     * explicit zero-init is removed because fixed_pairwise_sum_indexed
     * returns 0.0 on empty lists, preserving the empty-lake semantics.
     * Fork-join structure unchanged; no schedule / atomic / reduction
     * pragmas added (P9 owns parallel attribution). */
    if(lakeon){
        for (i = 0; i < NumLake; i++) {
            qLakeEvap[i] = fixed_pairwise_sum_indexed(
                    ele_by_lake[i], qEleEvapo_lake);
            qLakePrcp[i] = fixed_pairwise_sum_indexed(
                    ele_by_lake[i], qElePrep_lake);
        }
    }
    for (i = 0; i < NumLake; i++) {
        qLakeEvap[i] = min(qLakeEvap[i], qLakePrcp[i] + yLakeStg[i]);
        qLakeEvap[i] = max(0, qLakeEvap[i]);
    }
    } /* end river bucket (incl. lake transitional gather + clamp) */
    /* #43 (S1-pre-B): before-PassValue_legacy probe. Dumps Qe2r_Surf
     * (length NumEle), a PassValue_legacy() write-set member that carries
     * the previous iteration's element-to-river surface flux state.
     * PassValue_legacy (see body at L182-205) zero-resets Qe2r_Surf[0..NumEle-1]
     * and then accumulates QsegSurf over NumSegmt segments; capturing
     * Qe2r_Surf HERE gives a deterministic snapshot of the value that
     * is about to be cleared + re-derived by PassValue_legacy. PR #54 round-1
     * fix F4 replaced the prior QeleSurfTot probe payload (which
     * f_update zero-resets so the snapshot was always all zeros and
     * thus useless as a before-vs-after PassValue_legacy diff) with this
     * write-set member. Site tag "f_loop_before_passvalue" is distinct
     * from the no-op "f_loop" hook below + the "f_update" hook in
     * MD_update.cpp:151; the writer SHUD_DUMP_FNAME_SUFFIX env
     * disambiguates output files (`snapshot_t<v>_before_passvalue.bin`
     * vs `snapshot_t<v>.bin`). SHUD_DUMP_RHS=0 builds emit zero code
     * (compile-switch neutrality contract). */
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop_before_passvalue", t, Qe2r_Surf, NumEle);
#endif
    /* Shared for both OpenMP and Serial, to update.
     * S3c.3 (PR-11 #155): retired the legacy in-line gather function
     * (formerly PassValue_legacy, MD_f.cpp); the replacement
     * rhs_deterministic_gather() consumes the 7 S4 adjacency lists
     * (PR-10) to do all segment->river/element + downstream river +
     * lake river-in/surf/sub gathering, with bitwise-preserved
     * iteration order. */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_gather(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_GATHER]);
#endif
        rhs_deterministic_gather();
    } /* end gather bucket */
#ifdef SHUD_DUMP_RHS
    shud_rhs_dump_point("f_loop", t, NULL, 0);
#endif
}

/* S3c.3 (PR-11 #155) -- `Model_Data::rhs_deterministic_gather` is the
 * unified deterministic gather called from `rhs_flux` at the prior
 * legacy in-line gather site (formerly `PassValue_legacy`). Per design.md
 * D12 the function body lives WITH the RHS core in MD_rhs_core.cpp
 * (NOT in a separate MD_gather.cpp file).
 *
 * Body content:
 *   1. Pre-zero all river / element / lake accumulators (River:
 *      QrivSurf, QrivSub, QrivUp; Element: Qe2r_Surf, Qe2r_Sub;
 *      Lake: QLakeRivIn, QLakeSurf, QLakeSub).
 *   2. S3c.1 -- segment -> river gather via S4.1 seg_by_riv.
 *   3. S3c.1 -- segment -> element gather via S4.2 seg_by_ele.
 *   4. S3c.2 -- downstream river -> upstream gather via S4.3
 *      upstream_by_down (predicate `toLake<=0` already baked into the
 *      list at build time, MD_adjacency.cpp L108-115).
 *   5. S3b.1 -- per-river-down -> per-lake gather via S4.4
 *      riv_in_by_lake (lake-only branch).
 *   6. S3b.2 -- per-element-edge surf -> per-lake gather via S4.6
 *      lake_bank_edge_by_lake.
 *   7. S3b.3 -- per-element-edge sub -> per-lake gather via S4.6
 *      lake_bank_edge_by_lake.
 *
 * NOT included here (and intentionally so):
 *   - S3b.4 qLakeEvap / qLakePrcp per-element->per-lake gather; that
 *     stays in `rhs_flux` BEFORE the lake clamp pass (the clamp reads
 *     the gathered values; constraint documented in PR-9 commit log).
 *
 * Bitwise reproducibility: each S4 list iterates in B0 ascending
 * array-index order (MD_adjacency.cpp), so the per-accumulator `+=`
 * sequence is identical to the legacy serial loops. No `#pragma omp
 * parallel` directive is added; OpenMP parallelization of the gather
 * is deferred to P1+ per master plan.
 */
void Model_Data::rhs_deterministic_gather(){
    /* -------- pre-zeros -------- */
    for (int i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    for (int i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }

    /* -------- S3c.1: segment -> river gather (S4.1, fixed-shape leftfold
     * canonical reduction over seg_by_riv, B0 traversal order preserved).
     * P1c PR-C (#246): sites 3-4 wrapped with fixed_leftfold_sum_indexed
     * (byte-equivalent to original `acc += src[iseg]` since the operation
     * order is identical). Pre-zero loop at function head kept for
     * explicitness; helper assignment overwrites (zero-cost behavior). */
    for (int ir = 0; ir < NumRiv; ir++) {
        QrivSurf[ir] = fixed_leftfold_sum_indexed(seg_by_riv[ir], QsegSurf);
        QrivSub[ir]  = fixed_leftfold_sum_indexed(seg_by_riv[ir], QsegSub);
    }

    /* -------- S3c.1: segment -> element gather (S4.2, leftfold +
     * post-negate; `-leftfold(...)` is IEEE-754 bitwise-equivalent to
     * original left-to-right `acc += -src[iseg]` since negation is exact
     * sign-bit flip). P1c PR-C (#246): sites 5-6. Pre-zero kept for
     * explicitness; helper assignment overwrites. */
    for (int ie = 0; ie < NumEle; ie++) {
        Qe2r_Surf[ie] = -fixed_leftfold_sum_indexed(seg_by_ele[ie], QsegSurf);
        Qe2r_Sub[ie]  = -fixed_leftfold_sum_indexed(seg_by_ele[ie], QsegSub);
    }

    /* -------- S3c.2: downstream river -> upstream gather (S4.3,
     * leftfold + post-negate; `-leftfold(...)` is IEEE-754
     * bitwise-equivalent to left-to-right `acc += -src[up]` since
     * negation is exact sign-bit flip). upstream_by_down[ir] was
     * built with both the `iDownStrm>=0` and `Riv[i].toLake<=0`
     * predicates already baked in. -------- */
    for (int ir = 0; ir < NumRiv; ir++) {
        QrivUp[ir] = -fixed_leftfold_sum_indexed(upstream_by_down[ir], QrivDown);
    }

    /* -------- lake-side gathers (lakeon-gated) -------- */
    if (lakeon) {
        /* S3b.1: per-river-down -> per-lake (S4.4 riv_in_by_lake,
         * fixed-shape leftfold canonical reduction over B0 ascending
         * iriv order). */
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeRivIn[ilake] = 0.;
        }
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeRivIn[ilake] = fixed_leftfold_sum_indexed(
                    riv_in_by_lake[ilake], QrivDown);
        }

        /* S3b.2: per-element-edge surface -> per-lake (S4.6
         * lake_bank_edge_by_lake, fixed-shape leftfold over B0 ascending
         * (iele, j) pairs; stride=3 element-edge convention). */
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSurf[ilake] = 0.;
        }
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSurf[ilake] = fixed_leftfold_sum_pair_indexed(
                    lake_bank_edge_by_lake[ilake], QeleSurf_lake, 3);
        }

        /* S3b.3: per-element-edge subsurface -> per-lake (same S4.6
         * list, separate accumulator). */
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSub[ilake] = 0.;
        }
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSub[ilake] = fixed_leftfold_sum_pair_indexed(
                    lake_bank_edge_by_lake[ilake], QeleSub_lake, 3);
        }
    }
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
        /* S5d.2-5a (#179) — flat read via accessor. */
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurfAt(i, j);
            QeleSubTot[i] += QeleSubAt(i, j);
            CheckNANij(QeleSurfAt(i, j), i, "QeleSurfAt(i, j)");
            CheckNANij(QeleSubAt(i, j), i, "QeleSubAt(i, j)");
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
            /* S5c-B (#174): bucket 0 (update) + bucket 6 (applyDY) are
             * wrapped at the dispatch seam; buckets 1-5 are wrapped
             * inside rhs_flux(). rhs_flux as a whole is NOT timed as a
             * single bucket — its 5 inner sub-phases collectively cover
             * it, so the sum of the 7 buckets equals the wall time of
             * rhs_core's Serial branch (modulo std::chrono overhead). */
            {
#ifdef SHUD_ENABLE_DIAGNOSTICS
                shud_diag::ScopeTimer _t_upd(
                    &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_UPDATE]);
#endif
                rhs_update(Y, DY, t);
            }
            rhs_flux(t);
            {
#ifdef SHUD_ENABLE_DIAGNOSTICS
                shud_diag::ScopeTimer _t_app(
                    &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_APPLYDY]);
#endif
                rhs_apply(DY, t);
            }
            break;
#ifdef SHUD_ENABLE_OPENMP_RHS
        case ExecPolicy::StrictOMP: {
            /* P1e PR-F (#314) — replaces S1d.1 std::abort() stub with the
             * design D2 single-region OpenMP impl. Three sequential phases
             * (rhs_update -> rhs_flux -> rhs_apply) share a single outer
             * `#pragma omp parallel`, with implicit barriers between phases
             * supplied by `omp single` (each `single` ends with a barrier
             * because no `nowait` is requested). Phase ordering and data
             * dependencies are identical to the Serial branch above.
             *
             * Scope contract per p1e-strict-omp-rhs tasks.md §3 + design D2:
             *   - Single outer parallel region (master plan C7 fork-join
             *     minimization)
             *   - 3 phases with implicit barriers between (Phase 1 update
             *     -> Phase 2 flux -> Phase 3 apply)
             *   - `default(none) shared(Y, DY, t) private()` — explicit
             *     data-sharing (Y/DY/t are the only crossing scalars; all
             *     entity state is reached through `this->` member access
             *     which is implicit-shared via the enclosing method)
             *   - `schedule(static)` is reserved for the `omp for`
             *     decomposition introduced in PR-G/I (per task 3.4 the
             *     PR-F-stage Makefile does not yet add -fopenmp, so
             *     pragmas compile to serial code; cross-N parallel
             *     verification is deferred to PR-G after the -fopenmp
             *     wiring in tasks 3.5/3.6)
             *   - PR-F structural minimum: phase boundary expressed at
             *     function-call granularity via `omp single` (only one
             *     thread runs the call; other threads stall at the
             *     implicit barrier). PR-G/I will progressively convert
             *     these `single` blocks into per-entity `omp for` loops
             *     once the steady-state first-touch loops (design D4)
             *     have been removed in PR-G — avoiding nested `omp
             *     parallel` regions that the existing inner first-touch
             *     directives would create from inside an outer parallel.
             *
             * Bitwise contract: with -fopenmp ON or OFF, the work executes
             * sequentially in the calling thread (without -fopenmp the
             * directives are no-ops; with -fopenmp the `single` clause
             * serializes the call). N=1 mode C output therefore equals
             * mode A bit-for-bit.
             *
             * P1e PR-G (#315 / design D3) — `num_threads(omp_get_max_threads())`
             * pins the team size to the value `omp_set_num_threads()` left
             * in the OpenMP ICVs at startup (driven by `SHUD_RHS_THREADS`,
             * defaulting to `omp_get_max_threads()` itself). The explicit
             * clause is defense against accidental nesting (e.g. a future
             * caller wrapping rhs_core in its own outer parallel region):
             * without the clause the inner team could inherit the outer
             * team size and over-subscribe. The PR-G design forbids
             * `omp_set_num_threads` inside this hot path (it has been
             * called exactly once at shud.cpp startup); grep guard
             * `grep -nE 'omp_set_num_threads' MD_rhs_core.cpp == 0`.
             */
            #pragma omp parallel default(none) shared(Y, DY, t) \
                num_threads(omp_get_max_threads())
            {
                /* Phase 1: rhs_update — element / river / lake owner-local
                 * update + DY zero-init. Implicit barrier at end of
                 * `single` (no nowait) ensures Phase 2 sees fully updated
                 * uYsf/uYus/uYgw/uYriv/yLakeStg + zero-initialized DY. */
                #pragma omp single
                {
                    rhs_update(Y, DY, t);
                }

                /* Phase 2: rhs_flux — element/river/lake flux computation
                 * + deterministic gather. Implicit barrier ensures Phase 3
                 * sees finalized Qe2r_Surf/Sub, QrivSurf/Sub/Up/Down, and
                 * lake gather buffers. */
                #pragma omp single
                {
                    rhs_flux(t);
                }

                /* Phase 3: rhs_apply — DY accumulation over NumY. */
                #pragma omp single
                {
                    rhs_apply(DY, t);
                }
            } /* end parallel region */
            break;
        }
        case ExecPolicy::ProductionOMP:
            /* ProductionOMP backend remains a P2+ scope stub. */
            std::abort();
#endif
        default:
            /* Catches enumerator additions not yet wired + ABI drift. */
            std::abort();
    }
}
