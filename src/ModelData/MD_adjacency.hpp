/* MD_adjacency.hpp
 *
 * Declares the 7 adjacency lists used by the deterministic gather.
 * Each list stores, per target (river / element / lake), the source
 * indices in the order of the plain serial loop (ascending array
 * index). `Model_Data::rhs_deterministic_gather()` sums fluxes via
 * `for k in list[i]` traversal, so the order of every sum is fixed
 * and the results are reproducible for any thread count.
 *
 * The lists are built once at init time by `build_adjacency_lists()`
 * (called from `Model_Data::initialize()`) and are read-only
 * afterwards.
 */

#ifndef MD_ADJACENCY_HPP
#define MD_ADJACENCY_HPP

#include <vector>
#include <utility>

class Model_Data;  // forward decl — full include would re-enter classes/

/* ---- seg_by_riv[ir] ----------------------------------------------------
 * Per-river segment indices; used to sum segment fluxes into each river.
 * Iteration: `for i in [0, NumSegmt): if RivSeg[i].iRiv - 1 == ir yield i`.
 * Sort rule: segment array index ascending.
 */
extern std::vector<std::vector<int>> seg_by_riv;

/* ---- seg_by_ele[ie] ----------------------------------------------------
 * Per-element segment indices; used to sum segment fluxes into each
 * element.
 * Iteration: `for i in [0, NumSegmt): if RivSeg[i].iEle - 1 == ie yield i`.
 * Sort rule: segment array index ascending.
 */
extern std::vector<std::vector<int>> seg_by_ele;

/* ---- upstream_by_down[ir] ----------------------------------------------
 * Per-downstream-river upstream river indices. Uses Riv[i].down
 * (= iDownStrm macro); rivers that drain into a lake are excluded
 * (they are listed in riv_in_by_lake instead).
 * Iteration: `for i in [0, NumRiv): if Riv[i].down - 1 == ir and
 *             Riv[i].toLake <= 0 yield i`.
 * Sort rule: river array index ascending.
 */
extern std::vector<std::vector<int>> upstream_by_down;

/* ---- riv_in_by_lake[ilake] ---------------------------------------------
 * Per-lake river-inflow indices; used to sum QrivDown into QLakeRivIn.
 * Pitfall: Riv[i].toLake is used as a 0-based lake index — no `-1` is
 * applied, and `Riv[i].toLake >= 0` is the "flows into a lake" test —
 * unlike Ele[i].iLake and Ele[i].lakenabr[j], which are 1-based.
 * Iteration: `for i in [0, NumRiv): if Riv[i].toLake >= 0 yield i into bucket toLake`.
 * Sort rule: river array index ascending.
 */
extern std::vector<std::vector<int>> riv_in_by_lake;

/* ---- ele_by_lake[ilake] ------------------------------------------------
 * Per-lake element indices (the lake-cell elements); used in rhs_flux
 * to sum per-element evaporation / precipitation into each lake.
 * Iteration: `for i in [0, NumEle): if Ele[i].iLake > 0 yield i into bucket iLake-1`.
 * Sort rule: element array index ascending.
 */
extern std::vector<std::vector<int>> ele_by_lake;

/* ---- lake_bank_edge_by_lake[ilake] -------------------------------------
 * Per-lake (element-index, edge-index) pairs: the element edges that
 * border the lake. The lake branches of fun_Ele_surface / fun_Ele_sub
 * write the element->lake flux of each such edge into a per-edge slot;
 * this list is used to sum those slots into each lake.
 * Iteration: `for i in [0, NumEle): for j in {0,1,2}: if Ele[i].lakenabr[j]-1 >= 0
 *             yield (i, j) into bucket lakenabr[j]-1`.
 * Sort rule: element array index ascending, j inner = 0,1,2.
 */
extern std::vector<std::vector<std::pair<int,int>>> lake_bank_edge_by_lake;

/* ---- edge_by_ele[ie] ---------------------------------------------------
 * Per-element 3-neighbor list. Used by per-element lateral flux gather.
 * Each element stores its 3 edge neighbors at fixed j = 0,1,2 order:
 * value = Ele[ie].nabr[j] - 1, or -1 if nabr[j] == 0 (boundary).
 * Sort rule: fixed j = 0,1,2.
 */
extern std::vector<std::vector<int>> edge_by_ele;

/* ---- assert / fallback bookkeeping -------------------------------------
 * 3 asserts (RivSeg / Riv / Ele all satisfy `index == array_index + 1`).
 * If any assert fails, the lists are still ordered by array index, not
 * sorted by id. The booleans below are populated by
 * `build_adjacency_lists()` so the host code + unit test can inspect
 * post-build status.
 *
 * The identifier checked is the class member `.index` — the SHUD entity
 * classes (RiverSegement, _River, _Element) all expose `index` as the
 * canonical 1-based identifier.
 */
extern bool adjacency_assert_rivseg_pass;
extern bool adjacency_assert_riv_pass;
extern bool adjacency_assert_ele_pass;
extern bool adjacency_fallback_triggered;  // true if any of the above failed

/* Build all 7 adjacency lists from MD's current entity arrays.
 *
 * Always uses array-index ordering for list construction (the equivalent
 * of id-sort if `index == array_index + 1` holds, which is the
 * documented invariant; otherwise array-index ordering is the required
 * fallback, because it is the order of the plain serial loop).
 * The 3 asserts are evaluated and exposed via the booleans above so
 * callers (and unit tests) can verify the invariants without
 * re-walking the entity arrays.
 *
 * Idempotent — calling multiple times produces the same result;
 * subsequent calls clear-and-rebuild. It is called exactly once, from
 * `Model_Data::initialize()` (after the entity arrays + counts,
 * including NumLake, are populated).
 *
 * Returns true if all 3 asserts pass; false if any assert failed (in
 * which case `adjacency_fallback_triggered = true`).
 */
bool build_adjacency_lists(Model_Data* MD);

#endif  /* MD_ADJACENCY_HPP */
