/* MD_adjacency.hpp — S4 PR-10 (openMP issue #154).
 *
 * Declares the 7 B1a deterministic-gather adjacency lists (S4.1–S4.7
 * per `openspec/changes/b1a-finalization/specs/s4-adjacency-topology/spec.md`
 * + master plan §5 L1283–L1291). Each list captures the B0 serial loop
 * iteration order so PR-11 (S3c) can drive `rhs_deterministic_gather()`
 * via `for k in list[i]` style traversal and preserve bitwise
 * reproducibility vs B0.
 *
 * PR-10 BUILDS the lists at init time but does NOT YET USE them at
 * runtime (PR-11 will replace PassValue_legacy's gather). PR-10's bitwise
 * vs B0 guarantee is therefore trivial: adding init-time data
 * structures cannot change runtime behavior.
 *
 * design.md D12 — file landing.
 * design.md D9  — manifest YAML schema (docs/topology_manifest.yaml).
 * master plan §5 L1304–L1313 — assert + fallback contract.
 */

#ifndef MD_ADJACENCY_HPP
#define MD_ADJACENCY_HPP

#include <vector>
#include <utility>

class Model_Data;  // forward decl — full include would re-enter classes/

/* ---- S4.1 — seg_by_riv[ir] ---------------------------------------------
 * Per-river segment indices. B0 source: MD_f.cpp:170-171 (legacy)
 * == MD_f.cpp PassValue_legacy() body L214-221 (post-PR-9 = post-S3a).
 * Iteration: `for i in [0, NumSegmt): if RivSeg[i].iRiv - 1 == ir yield i`.
 * Sort rule: B0 iseg array index ascending.
 */
extern std::vector<std::vector<int>> seg_by_riv;

/* ---- S4.2 — seg_by_ele[ie] ---------------------------------------------
 * Per-element segment indices. B0 source: MD_f.cpp:172-173 (legacy)
 * == same loop body L214-221.
 * Iteration: `for i in [0, NumSegmt): if RivSeg[i].iEle - 1 == ie yield i`.
 * Sort rule: B0 iseg array index ascending.
 */
extern std::vector<std::vector<int>> seg_by_ele;

/* ---- S4.3 — upstream_by_down[ir] ---------------------------------------
 * Per-downstream-river upstream river indices. B0 source: MD_f.cpp:177
 * == PassValue_legacy() L222-226 (post-PR-9). Uses Riv[i].down (= iDownStrm macro).
 * Iteration: `for i in [0, NumRiv): if Riv[i].down - 1 == ir yield i`.
 * Sort rule: B0 iriv array index ascending.
 */
extern std::vector<std::vector<int>> upstream_by_down;

/* ---- S4.4 — riv_in_by_lake[ilake] --------------------------------------
 * Per-lake river-inflow indices. B0 source: MD_RiverFlux.cpp Flux_RiverDown
 * old shared-write (`QLakeRivIn[Riv[i].toLake] += QrivDown[i]`); post-PR-9
 * this lives in PassValue_legacy (MD_f.cpp L233-241). Riv[i].toLake is currently
 * 0-indexed in the post-PR-9 code (PassValue_legacy L238-239 uses it as direct
 * index into QLakeRivIn — no `-1` applied, with the guard
 * `Riv[i].toLake >= 0` for "no lake"). The spec L65 + master plan §5 L1288
 * describe the conceptual mapping as `riv_in_by_lake[ilake]` where
 * `ilake = toLake`. PR-10 list build follows the SAME predicate as the
 * post-PR-9 runtime — `Riv[i].toLake >= 0` yields i into bucket
 * `riv_in_by_lake[Riv[i].toLake]` (no `-1` applied).
 * Iteration: `for i in [0, NumRiv): if Riv[i].toLake >= 0 yield i into bucket toLake`.
 * Sort rule: B0 iriv array index ascending.
 */
extern std::vector<std::vector<int>> riv_in_by_lake;

/* ---- S4.5 — ele_by_lake[ilake] -----------------------------------------
 * Per-lake element indices. B0 source: MD_f.cpp:15-16 (legacy) ==
 * MD_rhs_core.cpp::rhs_flux L160-161 + L210-213 (post-PR-9).
 * Iteration: `for i in [0, NumEle): if Ele[i].iLake > 0 yield i into bucket iLake-1`.
 * Sort rule: B0 iele array index ascending.
 */
extern std::vector<std::vector<int>> ele_by_lake;

/* ---- S4.6 — lake_bank_edge_by_lake[ilake] ------------------------------
 * Per-lake (element-index, edge-index) pairs. B0 source: fun_Ele_surface +
 * fun_Ele_sub lake branches; the gather equivalent post-PR-9 lives in
 * PassValue_legacy (MD_f.cpp L248-271).
 * Iteration: `for i in [0, NumEle): for j in {0,1,2}: if Ele[i].lakenabr[j]-1 >= 0
 *             yield (i, j) into bucket lakenabr[j]-1`.
 * Sort rule: B0 iele ascending, j inner = 0,1,2.
 */
extern std::vector<std::vector<std::pair<int,int>>> lake_bank_edge_by_lake;

/* ---- S4.7 — edge_by_ele[ie] --------------------------------------------
 * Per-element 3-neighbor list. Used by per-element lateral flux gather.
 * Each element stores its 3 edge neighbors at fixed j = 0,1,2 order:
 * value = Ele[ie].nabr[j] - 1, or -1 if nabr[j] == 0 (boundary).
 * Sort rule: fixed j = 0,1,2.
 */
extern std::vector<std::vector<int>> edge_by_ele;

/* ---- assert / fallback bookkeeping -------------------------------------
 * 3 asserts (RivSeg / Riv / Ele all satisfy `index == array_index + 1`).
 * If any assert fails, build falls back to array-index ordering rather
 * than id-sort (spec L106-118 + master plan §5 L1313). The booleans
 * below are populated by `build_adjacency_lists()` so the host code +
 * unit test can inspect post-build status.
 *
 * S4 spec uses class member `.index` (not `.id`) — the SHUD entity
 * classes (RiverSegement, _River, _Element) all expose `index` as the
 * canonical 1-based identifier (River.hpp L48 / L97, Element.hpp L70).
 */
extern bool adjacency_assert_rivseg_pass;
extern bool adjacency_assert_riv_pass;
extern bool adjacency_assert_ele_pass;
extern bool adjacency_fallback_triggered;  // true if any of the above failed

/* Build all 7 adjacency lists from MD's current entity arrays.
 *
 * Always uses array-index ordering for list construction (the equivalent
 * of id-sort if `index == array_index + 1` holds, which is the
 * documented invariant; otherwise array-index ordering is the spec's
 * mandated fallback path). The 3 asserts are evaluated and exposed via
 * the booleans above so callers (and unit tests) can verify per-case
 * invariants without re-walking the entity arrays.
 *
 * Idempotent — calling multiple times produces the same result;
 * subsequent calls clear-and-rebuild. PR-10 calls this exactly once
 * from `Model_Data::initialize()` (after `malloc_EleRiv()` so the entity
 * arrays + counts are populated).
 *
 * Returns true if all 3 asserts pass; false if any assert failed (in
 * which case `adjacency_fallback_triggered = true`).
 */
bool build_adjacency_lists(Model_Data* MD);

#endif  /* MD_ADJACENCY_HPP */
