/* MD_adjacency.cpp — S4 PR-10 (openMP issue #154).
 *
 * Builds 7 B1a deterministic-gather adjacency lists at init time. See
 * MD_adjacency.hpp top-of-file comment for per-list semantics +
 * b0_source references.
 *
 * PR-10 BUILDS the lists but does NOT YET USE them at runtime. PR-11
 * (S3c) will replace `Model_Data::PassValue()`'s in-loop gather with
 * `rhs_deterministic_gather()` that iterates by list index. PR-10's
 * bitwise vs B0 guarantee is trivial: init-time-only changes cannot
 * affect runtime output.
 */

#include "MD_adjacency.hpp"
#include "Model_Data.hpp"
#include <cstdio>
#include <cstdlib>  /* getenv */

/* --- definitions of the 7 extern lists + assert booleans --------------- */
std::vector<std::vector<int>> seg_by_riv;
std::vector<std::vector<int>> seg_by_ele;
std::vector<std::vector<int>> upstream_by_down;
std::vector<std::vector<int>> riv_in_by_lake;
std::vector<std::vector<int>> ele_by_lake;
std::vector<std::vector<std::pair<int,int>>> lake_bank_edge_by_lake;
std::vector<std::vector<int>> edge_by_ele;

bool adjacency_assert_rivseg_pass     = true;
bool adjacency_assert_riv_pass        = true;
bool adjacency_assert_ele_pass        = true;
bool adjacency_fallback_triggered     = false;

bool build_adjacency_lists(Model_Data* MD){
    /* Reset all 7 lists + bookkeeping (idempotent). */
    seg_by_riv.assign((MD->NumRiv > 0 ? MD->NumRiv : 0), std::vector<int>());
    seg_by_ele.assign((MD->NumEle > 0 ? MD->NumEle : 0), std::vector<int>());
    upstream_by_down.assign((MD->NumRiv > 0 ? MD->NumRiv : 0), std::vector<int>());
    riv_in_by_lake.assign((MD->NumLake > 0 ? MD->NumLake : 0), std::vector<int>());
    ele_by_lake.assign((MD->NumLake > 0 ? MD->NumLake : 0), std::vector<int>());
    lake_bank_edge_by_lake.assign((MD->NumLake > 0 ? MD->NumLake : 0),
                                  std::vector<std::pair<int,int>>());
    edge_by_ele.assign((MD->NumEle > 0 ? MD->NumEle : 0), std::vector<int>());

    adjacency_assert_rivseg_pass = true;
    adjacency_assert_riv_pass    = true;
    adjacency_assert_ele_pass    = true;
    adjacency_fallback_triggered = false;

    /* ---- 3 asserts (id == array_index + 1) ---------------------------- */
    /* SHUD entity classes use `.index` as the canonical 1-based identifier
     * (the spec's `id` field). The check verifies the invariant referenced
     * by master plan §5 L1304-1313. Note: lists are ALWAYS built in
     * array-index order regardless of assert pass/fail — when the invariant
     * holds, array-index order == id-sort order; when it fails (theoretical
     * fallback path), array-index order is the mandated correct ordering
     * (spec L106-118). */
    for (int i = 0; i < MD->NumSegmt; ++i){
        if (MD->RivSeg[i].index != i + 1){
            adjacency_assert_rivseg_pass = false;
            adjacency_fallback_triggered = true;
            break;  /* one failure marks the assert; loop body iteration
                     * cost is already paid up to here. */
        }
    }
    for (int i = 0; i < MD->NumRiv; ++i){
        if (MD->Riv[i].index != i + 1){
            adjacency_assert_riv_pass = false;
            adjacency_fallback_triggered = true;
            break;
        }
    }
    for (int i = 0; i < MD->NumEle; ++i){
        if (MD->Ele[i].index != i + 1){
            adjacency_assert_ele_pass = false;
            adjacency_fallback_triggered = true;
            break;
        }
    }

    /* ---- S4.1 — seg_by_riv -------------------------------------------- */
    /* B0 source: legacy MD_f.cpp:170-171 (deleted in PR-9 S3a but the
     * conceptual gather lives on in PassValue L214-221). Iteration: outer
     * `for i in [0, NumSegmt)` produces ascending array-index order. */
    for (int i = 0; i < MD->NumSegmt; ++i){
        int ir = MD->RivSeg[i].iRiv - 1;
        if (ir >= 0 && ir < MD->NumRiv){
            seg_by_riv[ir].push_back(i);
        }
    }

    /* ---- S4.2 — seg_by_ele -------------------------------------------- */
    /* B0 source: legacy MD_f.cpp:172-173. Same outer iteration as S4.1. */
    for (int i = 0; i < MD->NumSegmt; ++i){
        int ie = MD->RivSeg[i].iEle - 1;
        if (ie >= 0 && ie < MD->NumEle){
            seg_by_ele[ie].push_back(i);
        }
    }

    /* ---- S4.3 — upstream_by_down -------------------------------------- */
    /* B0 source: legacy MD_f.cpp:177 == PassValue() L222-226 (post-PR-9).
     * Macros.hpp:49 `#define iDownStrm Riv[i].down - 1`; predicate
     * `iDownStrm >= 0` ⇒ skip rivers whose down sentinel is -1 / 0
     * boundary. Iteration: outer `for i in [0, NumRiv)`. The legacy code
     * also gated on `Riv[i].toLake <= 0` for the upstream sum (lakes are
     * accumulated separately by S4.4); we replicate that here so the list
     * mirrors B0 iteration scope. */
    for (int i = 0; i < MD->NumRiv; ++i){
        int idown = MD->Riv[i].down - 1;
        if (idown >= 0 && MD->Riv[i].toLake <= 0){
            if (idown < MD->NumRiv){
                upstream_by_down[idown].push_back(i);
            }
        }
    }

    /* ---- S4.4 — riv_in_by_lake ---------------------------------------- */
    /* B0 source: PassValue() L237-241 (post-PR-9). The post-PR-9 code uses
     * `Riv[i].toLake` as a 0-based lake index (no `-1` applied; guard
     * `Riv[i].toLake >= 0` filters non-lake-bound rivers). Replicate
     * exactly so the list iteration order matches the runtime gather. */
    for (int i = 0; i < MD->NumRiv; ++i){
        int ilake = MD->Riv[i].toLake;
        if (ilake >= 0 && ilake < MD->NumLake){
            riv_in_by_lake[ilake].push_back(i);
        }
    }

    /* ---- S4.5 — ele_by_lake ------------------------------------------- */
    /* B0 source: rhs_flux L160-161 + L210-213 (post-PR-9). Predicate
     * `Ele[i].iLake > 0` ⇒ this element is a lake interior cell; bucket
     * key = iLake - 1 (1-based to 0-based). Outer iteration: array-index
     * ascending. */
    for (int i = 0; i < MD->NumEle; ++i){
        if (MD->Ele[i].iLake > 0){
            int ilake = MD->Ele[i].iLake - 1;
            if (ilake < MD->NumLake){
                ele_by_lake[ilake].push_back(i);
            }
        }
    }

    /* ---- S4.6 — lake_bank_edge_by_lake -------------------------------- */
    /* B0 source: PassValue() L251-258 + L264-271 (post-PR-9 lake-bank
     * gather pattern). 2-D iteration: outer i ∈ [0, NumEle), inner j ∈
     * {0,1,2}; predicate `Ele[i].lakenabr[j] - 1 >= 0` ⇒ append (i, j)
     * to bucket key `lakenabr[j] - 1`. */
    for (int i = 0; i < MD->NumEle; ++i){
        for (int j = 0; j < 3; ++j){
            int ilake = MD->Ele[i].lakenabr[j] - 1;
            if (ilake >= 0 && ilake < MD->NumLake){
                lake_bank_edge_by_lake[ilake].push_back(std::make_pair(i, j));
            }
        }
    }

    /* ---- S4.7 — edge_by_ele ------------------------------------------- */
    /* B0 source: the canonical 3-neighbor loop used inside fun_Ele_* /
     * f_applyDY. Each element gets a fixed-length 3 list; value is the
     * neighbor element index (1-based `nabr[j] - 1`); j = 0 is on the
     * boundary if `nabr[j] == 0`, encoded as -1. */
    for (int i = 0; i < MD->NumEle; ++i){
        edge_by_ele[i].reserve(3);
        for (int j = 0; j < 3; ++j){
            int inabr = MD->Ele[i].nabr[j] - 1;
            edge_by_ele[i].push_back(inabr);  /* -1 when nabr[j] == 0 (boundary) */
        }
    }

    /* Optional one-shot manifest probe: when SHUD_ADJACENCY_LOG is set
     * in the env, emit a single line to stderr that tools can grep for
     * to populate `docs/topology_manifest.yaml` `asserts` rows. Init-time
     * only — does NOT affect runtime hot-loops; bitwise neutrality
     * preserved (no stdout, no .dat side-effect). Used in PR-10
     * verification flow to confirm the 4 Mac cases + qhh all have
     * `index == i + 1` (asserts pass, no fallback). */
    if (std::getenv("SHUD_ADJACENCY_LOG") != nullptr){
        std::fprintf(stderr,
            "[SHUD_ADJACENCY] NumSegmt=%d NumRiv=%d NumEle=%d NumLake=%d "
            "rivseg_pass=%d riv_pass=%d ele_pass=%d fallback=%d\n",
            MD->NumSegmt, MD->NumRiv, MD->NumEle, MD->NumLake,
            adjacency_assert_rivseg_pass ? 1 : 0,
            adjacency_assert_riv_pass    ? 1 : 0,
            adjacency_assert_ele_pass    ? 1 : 0,
            adjacency_fallback_triggered ? 1 : 0);
    }

    return !adjacency_fallback_triggered;
}
