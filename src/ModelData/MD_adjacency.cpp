/* MD_adjacency.cpp
 *
 * Builds the 7 deterministic-gather adjacency lists at init time. See
 * MD_adjacency.hpp for per-list semantics.
 *
 * The lists are consumed at runtime by
 * `Model_Data::rhs_deterministic_gather()` and `rhs_flux()`
 * (MD_rhs_core.cpp), which sum fluxes by iterating each list in its
 * stored order.
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
    /* SHUD entity classes use `.index` as the canonical 1-based identifier.
     * The check verifies the invariant `index == array_index + 1`.
     * Note: lists are ALWAYS built in array-index order regardless of
     * assert pass/fail — when the invariant holds, array-index order ==
     * id-sort order; when it fails (theoretical fallback path),
     * array-index order is still the correct ordering. */
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

    /* ---- seg_by_riv --------------------------------------------------- */
    /* Iteration: outer `for i in [0, NumSegmt)` produces ascending
     * array-index order. */
    for (int i = 0; i < MD->NumSegmt; ++i){
        int ir = MD->RivSeg[i].iRiv - 1;
        if (ir >= 0 && ir < MD->NumRiv){
            seg_by_riv[ir].push_back(i);
        }
    }

    /* ---- seg_by_ele --------------------------------------------------- */
    /* Same outer iteration as seg_by_riv. */
    for (int i = 0; i < MD->NumSegmt; ++i){
        int ie = MD->RivSeg[i].iEle - 1;
        if (ie >= 0 && ie < MD->NumEle){
            seg_by_ele[ie].push_back(i);
        }
    }

    /* ---- upstream_by_down --------------------------------------------- */
    /* Macros.hpp `#define iDownStrm Riv[i].down - 1`; predicate
     * `iDownStrm >= 0` ⇒ skip rivers whose down sentinel is -1 / 0
     * boundary. Iteration: outer `for i in [0, NumRiv)`. A river is
     * counted as upstream inflow only if `Riv[i].toLake <= 0`; rivers
     * draining into a lake are accumulated separately via
     * riv_in_by_lake. */
    for (int i = 0; i < MD->NumRiv; ++i){
        int idown = MD->Riv[i].down - 1;
        if (idown >= 0 && MD->Riv[i].toLake <= 0){
            if (idown < MD->NumRiv){
                upstream_by_down[idown].push_back(i);
            }
        }
    }

    /* ---- riv_in_by_lake ----------------------------------------------- */
    /* `Riv[i].toLake` is a 0-based lake index (no `-1` applied; guard
     * `Riv[i].toLake >= 0` filters non-lake-bound rivers). */
    for (int i = 0; i < MD->NumRiv; ++i){
        int ilake = MD->Riv[i].toLake;
        if (ilake >= 0 && ilake < MD->NumLake){
            riv_in_by_lake[ilake].push_back(i);
        }
    }

    /* ---- ele_by_lake -------------------------------------------------- */
    /* Predicate
     * `Ele[i].iLake > 0` ⇒ this element is a lake cell; bucket
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

    /* ---- lake_bank_edge_by_lake --------------------------------------- */
    /* 2-D iteration: outer i ∈ [0, NumEle), inner j ∈
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

    /* ---- edge_by_ele -------------------------------------------------- */
    /* Same order as the 3-neighbor loop used inside fun_Ele_* /
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

    /* Optional one-shot diagnostic: when SHUD_ADJACENCY_LOG is set
     * in the env, emit a single line to stderr with the entity counts
     * and the assert / fallback status, so a project can be checked for
     * `index == i + 1` (asserts pass, no fallback). Init-time only —
     * does NOT affect runtime hot-loops or model output (no stdout, no
     * .dat side-effect). */
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
