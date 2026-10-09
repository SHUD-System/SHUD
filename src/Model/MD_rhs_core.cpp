/* MD_rhs_core.cpp — the coupled right-hand side (RHS) of the ODE system.
 *
 * `Model_Data::rhs_core()` runs three phases in order:
 *   1. `rhs_update()` — copy the state vector Y into the working
 *      arrays, apply boundary conditions, zero the flux accumulators
 *      and DY.
 *   2. `rhs_flux()`   — compute element / segment / river / lake fluxes
 *      and gather them per owner.
 *   3. `rhs_apply()`  — assemble DY from the fluxes.
 *
 * The floating-point operation order inside each phase is fixed; the
 * results must not depend on the build variant or the thread count.
 *
 * The `SHUD_DUMP_RHS` hook tags (`"f_update"`,
 * `"f_loop_before_passvalue"`, `"f_loop"`, `"f_applyDY"`) are keys used
 * by external snapshot-comparison tooling; do not rename them to match
 * the function names.
 *
 * `timeNow` is NOT written here — f() in f.cpp assigns it once before
 * calling `rhs_core()`.
 */
#include "MD_rhs_core.hpp"
#include "MD_adjacency.hpp"  /* Adjacency lists consumed by the gathers
                              * in rhs_flux() and
                              * rhs_deterministic_gather(). They are
                              * built once from
                              * `Model_Data::initialize()`
                              * (MD_initialize.cpp). */
#include <cstdlib>   /* std::abort -- unimplemented ExecPolicy cases */
#ifdef SHUD_DUMP_RHS
#include "MD_rhs_dump.h"
#endif
/* RHS 7-bucket diagnostic timer. The header is empty unless
 * SHUD_ENABLE_DIAGNOSTICS is defined, so the default build has no
 * extra code in the hot path and no floating-point change. */
#include "MD_diagnostics.hpp"

#ifdef SHUD_ENABLE_DIAGNOSTICS
namespace shud_diag {
/* Definition of the 7-bucket nanosecond accumulators declared in
 * MD_diagnostics.hpp. Plain globals: race-free when the RHS runs
 * serially, racy under the OpenMP RHS (see rhs_flux()). */
long long g_rhs_timer_ns[RHS_BUCKET_COUNT] = {0, 0, 0, 0, 0, 0, 0};
}  // namespace shud_diag
#endif

/* NUMA first-touch gate, defined and set once at startup in shud.cpp.
 * It is declared here but not consulted in this file: first-touch
 * happens once at allocation time (Model_Data.cpp::malloc_EleRiv) and
 * at load time (MD_initialize.cpp::LoadIC), outside the RHS hot path.
 * Do NOT add first-touch `#pragma omp parallel for` loops inside the
 * RHS: under `ExecPolicy::StrictOMP` the whole RHS body already runs
 * inside one parallel region, so they would create a nested region. */
extern int g_numa_first_touch_enabled;

void Model_Data::rhs_update(double *Y, double *DY, double t){
    /* Each top-level loop below carries `#pragma omp for
     * schedule(static)`. When invoked from `ExecPolicy::StrictOMP` the
     * directives work-share the iteration space across the enclosing
     * team; when invoked from `ExecPolicy::Serial` (no enclosing
     * parallel region) each `omp for` is an orphaned construct that
     * runs in the encountering thread, so the serial result is
     * unchanged bit-for-bit. The five upstream loops use `nowait`
     * because their owner-local writes target disjoint slots; only
     * the final `for (i < NumY) DY[i] = 0.` omits `nowait` so its
     * implicit barrier closes Phase 1 before rhs_flux reads DY. */
    #pragma omp for schedule(static) nowait
    for (int i = 0; i < NumEle; i++) {
//        uYsf[i] = (Y[iSF] >= 0.) ? Y[iSF] : 0.;
//        uYus[i] = (Y[iUS] >= 0.) ? Y[iUS] : 0.;
        /* Zero the per-edge fluxes through the flat-array accessors. */
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

    #pragma omp for schedule(static) nowait
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

    #pragma omp for schedule(static) nowait
    for (int i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    #pragma omp for schedule(static) nowait
    for (int i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }

    #pragma omp for schedule(static) nowait
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
    /* Final owner-local omp for in Phase 1. Carries no `nowait` so its
     * implicit barrier synchronises all upstream Phase-1 owner-local
     * writes before Phase 2 (rhs_flux) reads them. The post-loop
     * SHUD_DUMP_RHS hook below is wrapped in `omp single` (its own
     * implicit barrier) when the hook is enabled. */
    #pragma omp for schedule(static)
    for (int i = 0; i < NumY; i++){
        DY[i] = 0.;
    }
#ifdef SHUD_DUMP_RHS
    /* The SHUD_DUMP_RHS hook performs file I/O; `omp single` makes
     * exactly one thread run it under StrictOMP. */
    #pragma omp single
    {
        shud_rhs_dump_point("f_update", t, DY, NumY);
    }
#endif
}

/* `Model_Data::rhs_flux` computes all fluxes for one RHS evaluation.
 *
 * The process order is a strict sequence; each step reads what the
 * previous ones wrote:
 *   1. Element pass 1  (lake updateLakeElement / fun_Ele_lakeVertical
 *                       + per-element qLakeEvap/qLakePrcp shares;
 *                       non-lake f_etFlux + updateElement +
 *                       fun_Ele_Infiltraion + fun_Ele_Recharge)
 *   2. Element pass 2  (lake fun_Ele_lakeHorizon; non-lake
 *                       fun_Ele_surface + fun_Ele_sub)
 *   3. Segment pass    (fun_Seg_surface + fun_Seg_sub)
 *   4. River pass      (Flux_RiverDown)
 *   5. Lake pass       (per-lake gather of qLakeEvap/qLakePrcp, then
 *                       min/max clamp on qLakeEvap)
 *   6. rhs_deterministic_gather() (zero-reset + per-owner gathers)
 * The two SHUD_DUMP_RHS hooks sit immediately before and after step 6.
 *
 * Dump tag strings `"f_loop_before_passvalue"` and `"f_loop"` are keys
 * used by external snapshot-comparison tooling; do not rename them. */

/* Fixed-shape pairwise tree reduction over an index list, visited in
 * list order. Tree shape depends only on the list length, NOT on the
 * thread count — used at per-owner accumulation sites to obtain a
 * thread-count-independent canonical sum.
 *
 * Empty list returns 0.0 (so callers need no explicit zero-init);
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

/* Fixed-shape leftfold canonical reduction over an index list. For
 * each src element src[idx[k]] in list order, accumulator runs
 *   `acc = 0; acc += src[idx[0]]; acc += src[idx[1]]; ...`.
 * Result is bitwise-identical to a plain serial `for x in list: dst +=
 * src[x]` accumulation (with dst starting at 0.0) since the operation
 * order is the same. Used for the segment -> river / element gathers,
 * the upstream river gather and the river -> lake gather, which must
 * reproduce the serial loop bit-for-bit. The per-lake qLakeEvap /
 * qLakePrcp aggregation uses the tree variant
 * `fixed_pairwise_sum_indexed` instead. */
static inline double fixed_leftfold_sum_indexed(
        const std::vector<int>& idx, const double* src) {
    double acc = 0.0;
    for (int i : idx) acc += src[i];
    return acc;
}

/* Fixed-shape leftfold canonical reduction over an index-pair list.
 * For each (ie, j) pair in list order, accumulator runs
 * `acc += src[ie * stride + j]`. Bitwise-identical to a plain serial
 * `for (auto& ej : list) acc += src[ej.first * stride + ej.second]`.
 * Used for the lake-bank gathers (lake_bank_edge_by_lake), whose
 * adjacency list contains (ie, j) pairs; lookup is the standard
 * stride=3 element-edge convention. Constant-stride keeps the helper
 * trivially inlinable.
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
    /* The 5 inner diagnostic buckets (ET / lateral / segment / river /
     * gather) are scoped via `shud_diag::ScopeTimer`; the explicit
     * `{ }` blocks scope each ScopeTimer to exactly one phase. With
     * `SHUD_ENABLE_DIAGNOSTICS` undefined the timers are compiled out
     * and the loop bodies are unchanged.
     *
     * Under StrictOMP the ScopeTimer constructions race on
     * `g_rhs_timer_ns[bucket]` (each team thread runs the ctor/dtor
     * and += the global). This is accepted: diagnostics builds are
     * for offline profiling only and their timings are not exact with
     * the OpenMP RHS. */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_et(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_ET]);
#endif
    /* Each iteration writes only to Ele[i] (AoS + hot SoA slot i) and
     * to per-element scratch (qEleEvapo_lake[i], qElePrep_lake[i],
     * qEleEvapo[i], qElePrep[i] via f_etFlux / fun_Ele_* — all
     * owner-local to i). The trailing implicit barrier (no `nowait`)
     * is REQUIRED: the lateral bucket below reads `hot.u_effKH[inabr]`
     * (neighbour element's hot SoA slot, populated here by
     * sync_hot_dynamic(i) after updateElement). Without the barrier
     * some neighbour SoA slots could still hold stale values when the
     * lateral pass reads them — a cross-thread data race. */
    #pragma omp for schedule(static)
    for (int i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            Ele[i].updateLakeElement();
            /* updateLakeElement() mutates AoS hot fields
             * (Ele[i].u_effKH = KsatH, etc.) but does NOT refresh
             * hot.<field>[i]; sync_hot_dynamic(i) is required so
             * subsequent consumers reading hot.u_effKH[i] see the
             * post-update value rather than the stale forcing blend
             * left by updateforcing(). */
            sync_hot_dynamic(i);
            fun_Ele_lakeVertical(i, t);
            /* The per-lake sums
             *   qLakeEvap[Ele[i].iLake-1] += qEleEvapo[i] / NumEleLake
             *   qLakePrcp[Ele[i].iLake-1] += qElePrep[i]  / NumEleLake
             * would be shared writes, so each element stores its share
             * (already divided by lake.NumEleLake) in its own slot. The
             * gather after the river loop, BEFORE the lake clamp below,
             * sums them to per-lake qLakeEvap / qLakePrcp. */
            qEleEvapo_lake[i] = qEleEvapo[i] / lake[Ele[i].iLake - 1].NumEleLake;
            qElePrep_lake[i]  = qElePrep[i]  / lake[Ele[i].iLake - 1].NumEleLake;
        }else{
            f_etFlux(i, t);
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            /* Same SoA refresh after updateElement(): the infiltration
             * / recharge code reads hot.u_effkInfi[i] etc. and must see
             * the post-update value. */
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
    /* Owner-local writes to QeleSurfAt(i,j) / QeleSubAt(i,j).
     * Trailing implicit barrier (no `nowait`) — segment / river
     * buckets below cannot start until lateral writes are visible
     * (rhs_deterministic_gather() further down reads QsegSurf via
     * seg_by_ele indirection that fans out to neighbour element
     * slots, and we want a clean phase boundary). */
    #pragma omp for schedule(static)
    for (int i = 0; i < NumEle; i++) {
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
    /* Owner-local writes to QsegSurf[i] / QsegSub[i]. Trailing
     * implicit barrier (no `nowait`) — rhs_deterministic_gather()
     * below reads QsegSurf / QsegSub via seg_by_riv / seg_by_ele
     * adjacency lookups; the gather entry must observe all segment
     * writes. */
    #pragma omp for schedule(static)
    for (int i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    } /* end segment bucket */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_riv(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_RIVER]);
#endif
    /* Flux_RiverDown(t, i) writes QrivDown[i] (owner-local) and
     * reads neighbour `uYriv[iDown]` / `Riv[iDown].depth` (Phase 1
     * -> Phase 2 barrier already synchronised those). Trailing
     * implicit barrier (no `nowait`) — rhs_deterministic_gather()
     * below reads QrivDown via upstream_by_down + riv_in_by_lake. */
    #pragma omp for schedule(static)
    for (int i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    /* Per-element -> per-lake gather of qLakeEvap / qLakePrcp. It must
     * run BEFORE the lake clamp below (the clamp reads both), so it
     * cannot live in rhs_deterministic_gather(), which is called
     * AFTER the clamp.
     *
     * Each lake sum is a fixed-shape pairwise tree reduction over
     * `ele_by_lake[ilake]`. MD_adjacency.cpp fills that list by an
     * ascending `for (i = 0; i < NumEle; i++) if Ele[i].iLake > 0`
     * traversal, which is the canonical order — do NOT re-sort. Tree
     * shape is determined by list length only, NOT by the thread
     * count. No explicit zero-init is needed because
     * fixed_pairwise_sum_indexed returns 0.0 on an empty list.
     *
     * The gather + clamp is a sequential block that depends on the
     * NumEle ET pass above (qEleEvapo_lake / qElePrep_lake). It is
     * bracketed in `#pragma omp single` so exactly one team thread
     * runs it; the trailing implicit barrier ensures all team
     * threads observe the gathered + clamped per-lake values before
     * the rhs_deterministic_gather() block below. */
    #pragma omp single
    {
        if(lakeon){
            for (int i = 0; i < NumLake; i++) {
                qLakeEvap[i] = fixed_pairwise_sum_indexed(
                        ele_by_lake[i], qEleEvapo_lake);
                qLakePrcp[i] = fixed_pairwise_sum_indexed(
                        ele_by_lake[i], qElePrep_lake);
            }
        }
        for (int i = 0; i < NumLake; i++) {
            qLakeEvap[i] = min(qLakeEvap[i], qLakePrcp[i] + yLakeStg[i]);
            qLakeEvap[i] = max(0, qLakeEvap[i]);
        }
    } /* end omp single — lake transitional gather + clamp */
    } /* end river bucket (incl. lake transitional gather + clamp) */
    /* Before-gather probe. Dumps Qe2r_Surf (length NumEle), which at
     * this point still holds the element-to-river surface flux of the
     * previous RHS evaluation: rhs_deterministic_gather() below
     * zero-resets it and re-accumulates QsegSurf into it, so this is a
     * snapshot of the value about to be cleared and re-derived.
     * (QeleSurfTot would be useless here — rhs_update zero-resets it,
     * so it is always all zeros at this point.) Site tag
     * "f_loop_before_passvalue" is distinct from the no-op "f_loop"
     * hook below and the "f_update" hook in rhs_update(); the
     * SHUD_DUMP_FNAME_SUFFIX environment variable read by the writer
     * disambiguates output files (`snapshot_t<v>_before_passvalue.bin`
     * vs `snapshot_t<v>.bin`).
     *
     * The hook performs file I/O, so `#pragma omp single` makes
     * exactly one thread run the dump under StrictOMP. Without
     * SHUD_DUMP_RHS both the directive and the hook are compiled
     * out. */
#ifdef SHUD_DUMP_RHS
    #pragma omp single
    {
        shud_rhs_dump_point("f_loop_before_passvalue", t, Qe2r_Surf, NumEle);
    }
#endif
    /* Shared for both OpenMP and Serial, to update.
     * rhs_deterministic_gather() consumes the adjacency lists to do
     * all segment->river/element + downstream river + lake
     * river-in/surf/sub gathering in a fixed iteration order. */
    {
#ifdef SHUD_ENABLE_DIAGNOSTICS
        shud_diag::ScopeTimer _t_gather(
            &shud_diag::g_rhs_timer_ns[shud_diag::RHS_BUCKET_GATHER]);
#endif
        /* rhs_deterministic_gather() contains its own
         * `#pragma omp for schedule(static)` directives; do NOT wrap
         * the call site in `omp single` (would defeat parallelisation). */
        rhs_deterministic_gather();
    } /* end gather bucket */
#ifdef SHUD_DUMP_RHS
    #pragma omp single
    {
        shud_rhs_dump_point("f_loop", t, NULL, 0);
    }
#endif
}

/* `Model_Data::rhs_deterministic_gather` is the unified per-owner
 * gather called at the end of `rhs_flux`.
 *
 * Body content:
 *   1. Pre-zero all river / element / lake accumulators (River:
 *      QrivSurf, QrivSub, QrivUp; Element: Qe2r_Surf, Qe2r_Sub;
 *      Lake: QLakeRivIn, QLakeSurf, QLakeSub).
 *   2. Segment -> river gather via seg_by_riv.
 *   3. Segment -> element gather via seg_by_ele.
 *   4. Downstream river -> upstream gather via upstream_by_down
 *      (predicate `toLake<=0` already baked into the list at build
 *      time, MD_adjacency.cpp).
 *   5. Per-river-down -> per-lake gather via riv_in_by_lake
 *      (lake-only branch).
 *   6. Per-element-edge surf -> per-lake gather via
 *      lake_bank_edge_by_lake.
 *   7. Per-element-edge sub -> per-lake gather via
 *      lake_bank_edge_by_lake.
 *
 * NOT included here (and intentionally so):
 *   - the qLakeEvap / qLakePrcp per-element->per-lake gather; that
 *     stays in `rhs_flux` BEFORE the lake clamp pass, because the
 *     clamp reads the gathered values.
 *
 * Bitwise reproducibility: each adjacency list is stored in ascending
 * array-index order (MD_adjacency.cpp), so the per-accumulator `+=`
 * sequence is identical to a plain serial loop over the source array.
 */
void Model_Data::rhs_deterministic_gather(){
    /* Each top-level loop carries `#pragma omp for schedule(static)`.
     * Only the loop over owners is work-shared; each owner's sum is
     * computed by one thread through the helper reductions
     * (fixed_leftfold_sum_indexed / fixed_leftfold_sum_pair_indexed),
     * which walk that owner's adjacency list in its stored order —
     * independent of which team thread evaluates it — so results are
     * bit-identical for any thread count.
     *
     * The lake-side branch (`if (lakeon)`) wraps three pairs of
     * `omp for` directives. The conditional itself is uniform across
     * all team threads (lakeon is a Model_Data member set during
     * initialisation), so all threads either enter or skip the
     * branch together — `omp for` semantics inside a conditional
     * that all threads encounter identically are well-defined.
     *
     * The trailing `nowait` annotations let upstream gathers retire
     * as soon as their work-share finishes; the final loop in the
     * function (sub-gather over lake_bank_edge_by_lake) carries the
     * implicit barrier that synchronises team state before
     * rhs_flux's caller (the StrictOMP case body) proceeds to
     * Phase 3 rhs_apply. */

    /* -------- pre-zeros -------- */
    #pragma omp for schedule(static) nowait
    for (int i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    #pragma omp for schedule(static)
    for (int i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }
    /* Implicit barrier above — downstream gathers READ pre-zeroed
     * slots, so the team must observe the zero writes before
     * proceeding. */

    /* -------- segment -> river gather (fixed-shape leftfold canonical
     * reduction over seg_by_riv, in list order; bitwise-equivalent to
     * a serial `acc += src[iseg]` loop). The pre-zero loop at function
     * head is kept for explicitness; the helper assignment overwrites
     * it. */
    #pragma omp for schedule(static) nowait
    for (int ir = 0; ir < NumRiv; ir++) {
        QrivSurf[ir] = fixed_leftfold_sum_indexed(seg_by_riv[ir], QsegSurf);
        QrivSub[ir]  = fixed_leftfold_sum_indexed(seg_by_riv[ir], QsegSub);
    }

    /* -------- segment -> element gather (seg_by_ele, leftfold +
     * post-negate; `-leftfold(...)` is IEEE-754 bitwise-equivalent to
     * a left-to-right `acc += -src[iseg]` since negation is exact
     * sign-bit flip). Pre-zero kept for explicitness; helper
     * assignment overwrites. */
    #pragma omp for schedule(static)
    for (int ie = 0; ie < NumEle; ie++) {
        Qe2r_Surf[ie] = -fixed_leftfold_sum_indexed(seg_by_ele[ie], QsegSurf);
        Qe2r_Sub[ie]  = -fixed_leftfold_sum_indexed(seg_by_ele[ie], QsegSub);
    }
    /* Implicit barrier above — the downstream gather below reads
     * QrivDown (populated by Flux_RiverDown upstream + finalised by
     * the prior river bucket); we must also synchronise across the
     * QrivUp writes before any lake-side branch reads them. */

    /* -------- downstream river -> upstream gather (upstream_by_down,
     * leftfold + post-negate; `-leftfold(...)` is IEEE-754
     * bitwise-equivalent to left-to-right `acc += -src[up]` since
     * negation is exact sign-bit flip). upstream_by_down[ir] was
     * built with both the `iDownStrm>=0` and `Riv[i].toLake<=0`
     * predicates already baked in. -------- */
    #pragma omp for schedule(static)
    for (int ir = 0; ir < NumRiv; ir++) {
        QrivUp[ir] = -fixed_leftfold_sum_indexed(upstream_by_down[ir], QrivDown);
    }

    /* -------- lake-side gathers (lakeon-gated) -------- */
    if (lakeon) {
        /* Per-river-down -> per-lake (riv_in_by_lake, fixed-shape
         * leftfold canonical reduction over ascending iriv order). */
        #pragma omp for schedule(static) nowait
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeRivIn[ilake] = 0.;
        }
        #pragma omp for schedule(static)
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeRivIn[ilake] = fixed_leftfold_sum_indexed(
                    riv_in_by_lake[ilake], QrivDown);
        }

        /* Per-element-edge surface -> per-lake
         * (lake_bank_edge_by_lake, fixed-shape leftfold over ascending
         * (iele, j) pairs; stride=3 element-edge convention). */
        #pragma omp for schedule(static) nowait
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSurf[ilake] = 0.;
        }
        #pragma omp for schedule(static)
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSurf[ilake] = fixed_leftfold_sum_pair_indexed(
                    lake_bank_edge_by_lake[ilake], QeleSurf_lake, 3);
        }

        /* Per-element-edge subsurface -> per-lake (same list,
         * separate accumulator). Final `omp for` in the function
         * — must NOT carry `nowait` so its implicit barrier
         * synchronises team state before returning to rhs_flux's
         * caller. */
        #pragma omp for schedule(static) nowait
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSub[ilake] = 0.;
        }
        #pragma omp for schedule(static)
        for (int ilake = 0; ilake < NumLake; ilake++) {
            QLakeSub[ilake] = fixed_leftfold_sum_pair_indexed(
                    lake_bank_edge_by_lake[ilake], QeleSub_lake, 3);
        }
    }
    /* When lakeon == false the final `omp for` in the function
     * is the downstream-gather loop above (which omits `nowait`
     * already), so the implicit-barrier requirement still holds.
     * When lakeon == true the final loop is the lake-sub gather
     * directly above (also no `nowait`). */
}

/* `Model_Data::rhs_apply` assembles DY from the fluxes computed by
 * `rhs_flux`.
 *
 * River DY uses the 3-step formula (divide by Riv[i].Length, clamp
 * to the available cross-section area, then `fun_dAtodY()`); a
 * simpler 1-step "/ u_TopArea" formula MUST NOT be used. Lake DY
 * divides by y2LakeArea[i] with no zero guard / epsilon; a lake with
 * zero top area is therefore not handled here.
 *
 * `SHUD_DUMP_RHS` hook tag string is `"f_applyDY"` (NOT `"rhs_apply"`)
 * — it is a key used by external snapshot-comparison tooling and must
 * not be renamed. */
void Model_Data::rhs_apply(double *DY, double t){
    /* `area / isf / ius / igw` are declared inside the NumEle loop
     * body, not at method scope, so each StrictOMP team thread holds
     * its own per-iteration locals; method-scope declarations would be
     * shared across the team and race under the `omp for` below. */
    #pragma omp for schedule(static) nowait
    for (int i = 0; i < NumEle; i++) {
        int isf = iSF;
        int ius = iUS;
        int igw = iGW;
        double area = Ele[i].area;
        QeleSurfTot[i] = Qe2r_Surf[i];
        QeleSubTot[i] = Qe2r_Sub[i];
        /* Read the per-edge fluxes through the flat-array accessors. */
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
    /* NumRiv loop writes only DY[iRIV] (owner-local). `nowait`
     * since the next NumLake loop touches a disjoint DY range
     * (iLAKE = i + 3 * NumEle + NumRiv). */
    #pragma omp for schedule(static) nowait
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
    /* NumLake loop writes only DY[iLAKE]. Final `omp for` in
     * rhs_apply — no `nowait` so its implicit barrier closes Phase 3
     * before the StrictOMP parallel region exits. */
    #pragma omp for schedule(static)
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
    /* The SHUD_DUMP_RHS hook performs file I/O; `omp single` makes
     * exactly one thread run it under StrictOMP. */
    #pragma omp single
    {
        shud_rhs_dump_point("f_applyDY", t, DY, 3 * NumEle + NumRiv + NumLake);
    }
#endif
}

/* `rhs_core` ExecPolicy dispatch.
 *
 *   - `ExecPolicy::Serial` runs `rhs_update` -> `rhs_flux` ->
 *     `rhs_apply` sequentially.
 *   - `ExecPolicy::StrictOMP` runs the same three phases inside one
 *     `#pragma omp parallel` region (see the case body).
 *   - `ExecPolicy::ProductionOMP` is not implemented and calls
 *     `std::abort()`. `assert(false)` must not be used instead:
 *     `-DNDEBUG` strips assert to a no-op and would let execution
 *     fall through silently to the next statement.
 *   - Without SHUD_ENABLE_OPENMP_RHS the OMP cases are `#ifdef`ed
 *     out of the translation unit, so the binary contains no
 *     OMP-path code.
 *   - `default:` branch also aborts to catch ABI drift / future
 *     enumerator additions that haven't been wired up here.
 *
 * No template specialization, no virtual dispatch — plain switch.
 * The policy is a compile-time constant at the only caller (f() in
 * f.cpp), so the compiler can eliminate the branch in optimized
 * builds.
 */
void Model_Data::rhs_core(double *Y, double *DY, double t, ExecPolicy policy){
    switch (policy) {
        case ExecPolicy::Serial:
            /* Diagnostic bucket 0 (update) + bucket 6 (applyDY) are
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
            /* Single-region OpenMP implementation. The three phases
             * (rhs_update -> rhs_flux -> rhs_apply) share ONE
             * `#pragma omp parallel` region, to avoid a fork-join per
             * loop. Phase ordering and data dependencies are identical
             * to the Serial branch above.
             *
             *   - `default(none) shared(Y, DY, t)` — explicit
             *     data-sharing (Y/DY/t are the only crossing scalars; all
             *     entity state is reached through `this->` member access
             *     which is implicit-shared via the enclosing method).
             *   - All threads enter the three RHS methods together. The
             *     per-entity top-level loops inside them carry
             *     `#pragma omp for schedule(static)`; `dynamic` /
             *     `guided` schedules are not allowed. The last `omp for`
             *     in each method omits `nowait`, so its implicit barrier
             *     synchronises team state at method return.
             *   - No `#pragma omp parallel` may appear inside the three
             *     methods: it would create a nested parallel region.
             *
             * Reproducibility: with -fopenmp OFF the `omp parallel`
             * and `omp for` directives are no-ops and the methods run
             * fully serial. With -fopenmp ON and `SHUD_RHS_THREADS=1`
             * each `omp for` runs on the single team member in the
             * canonical serial order, so the output equals the serial
             * build bit-for-bit. For N>1 the per-iteration owner-local
             * writes (each thread updates a disjoint Ele[i] / Riv[i]
             * / Lake[i] / DY[i] slot) and the canonical leftfold /
             * pairwise reductions in rhs_flux() and
             * rhs_deterministic_gather() (which iterate over the
             * adjacency lists, NOT over OMP-decomposed ranges) keep
             * the result bit-identical across thread counts.
             *
             * `num_threads(omp_get_max_threads())` pins the team size
             * to the value `omp_set_num_threads()` set at startup in
             * shud.cpp (driven by `SHUD_RHS_THREADS`). The explicit
             * clause is defense against accidental nesting (e.g. a
             * future caller wrapping rhs_core in its own parallel
             * region): without the clause the inner team could
             * inherit the enclosing team size and over-subscribe. Do
             * not call `omp_set_num_threads` in this hot path; it is
             * called once at startup.
             */
            #pragma omp parallel default(none) shared(Y, DY, t) \
                num_threads(omp_get_max_threads())
            {
                /* Phase 1: rhs_update — element / river / lake owner-local
                 * update + DY zero-init. All threads enter the call; the
                 * top-level loops inside rhs_update carry `#pragma omp for
                 * schedule(static)`. The trailing implicit barrier at the
                 * end of the last `omp for` (the DY zero-init NumY loop)
                 * synchronises team state before Phase 2 reads it. */
                rhs_update(Y, DY, t);

                /* Phase 2: rhs_flux — element/river/lake flux computation
                 * + deterministic gather. Inner `omp for` loops + trailing
                 * implicit barrier in rhs_deterministic_gather()
                 * synchronise team state before Phase 3 reads the finalized
                 * Qe2r_Surf/Sub, QrivSurf/Sub/Up/Down, and lake gather
                 * buffers. */
                rhs_flux(t);

                /* Phase 3: rhs_apply — DY accumulation over NumY. Inner
                 * `omp for` loops decompose the per-element / per-river /
                 * per-lake DY writes; trailing implicit barrier on the
                 * final lake loop closes the parallel region. */
                rhs_apply(DY, t);
            } /* end parallel region */
            break;
        }
        case ExecPolicy::ProductionOMP:
            /* ProductionOMP is not implemented. */
            std::abort();
#endif
        default:
            /* Catches enumerator additions not yet wired + ABI drift. */
            std::abort();
    }
}
