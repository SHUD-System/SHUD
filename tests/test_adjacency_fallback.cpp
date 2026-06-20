/* tests/test_adjacency_fallback.cpp — S4 PR-10 (openMP issue #154).
 *
 * Standalone unit test for the `id != index + 1` fallback path inside
 * `build_adjacency_lists()` (MD_adjacency.cpp).
 *
 * Purpose: enforce the spec contract (s4-adjacency-topology Scenario
 * "三条 assert 在所有 6 case 都 pass + fallback 单测 PASS") that when
 * the `index == array_index + 1` invariant is violated on any entity
 * type, the build sets `adjacency_fallback_triggered = true` AND falls
 * back to the array-index ordering (rather than id-sort) so the
 * adjacency list output remains equivalent to the array-index-ordered
 * expected output.
 *
 * Mock strategy:
 *   - construct a `Model_Data` on the heap (`new`); default ctor body is
 *     empty (Model_Data.cpp:4-5), so all members hold their in-class
 *     default values (no UB)
 *   - manually set NumSegmt / NumRiv / NumEle / NumLake on the mock
 *   - manually allocate RivSeg / Riv / Ele arrays (_River / _Element /
 *     RiverSegement default-construct cleanly because of their in-class
 *     `= NA_VALUE` initializers)
 *   - set `.index` fields with a deliberate offset (`= i + 100`) so the
 *     spec's `index == i + 1` assert fires; for each case also set
 *     iRiv / iEle / iLake / down / toLake / nabr / lakenabr so the
 *     adjacency lists have a known expected shape
 *   - call `build_adjacency_lists(mock)`; assert
 *     `adjacency_fallback_triggered == true` for each of the 3
 *     synthetic mocks (one per assertable entity)
 *   - verify the lists match array-index ordering (NOT id-sort) using
 *     the canonical reference algorithm (= `for (i = 0; i < N; ++i)`
 *     with predicate match yields i in ascending array order)
 *   - destructor invocation: we DO NOT delete the mock — Model_Data's
 *     destructor calls FreeData() which deletes many unrelated arrays
 *     that we have NOT allocated (would UB on uninitialized pointers).
 *     `exit(0)` on success is the canonical pattern (matches
 *     s1d_strictomp_assert_smoke.cpp child-process path).
 *
 * Build / run: wired via the `test_adjacency_fallback` Makefile target.
 * CI: `.github/workflows/serial-baseline.yml` calls
 *     `make test_adjacency_fallback && ./tests/test_adjacency_fallback`.
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <utility>

#include "Model_Data.hpp"
#include "MD_adjacency.hpp"

namespace {

/* Test scaffold: build a synthetic Model_Data with a single non-canonical
 * `index` (= array_index + 100) and pre-populated entity arrays that
 * encode a known adjacency pattern. Returns the heap-allocated mock so
 * the caller can call `build_adjacency_lists()` and inspect results.
 *
 * Mesh shape (3 elements / 2 rivers / 4 segments / 1 lake):
 *
 *   ele 0 (iLake=0, lakenabr={1,0,0}, nabr={2,3,0})  ← bank-edge to lake 0 on j=0
 *   ele 1 (iLake=1, lakenabr={0,0,0}, nabr={3,0,1})  ← interior of lake 0
 *   ele 2 (iLake=0, lakenabr={0,1,0}, nabr={1,0,0})  ← bank-edge to lake 0 on j=1
 *
 *   riv 0 (down=2, toLake=-1)                         ← upstream of riv 1
 *   riv 1 (down=-1, toLake=0)                         ← flows into lake 0
 *
 *   seg 0 → riv 0, ele 1
 *   seg 1 → riv 0, ele 0
 *   seg 2 → riv 1, ele 2
 *   seg 3 → riv 1, ele 0
 *
 * Expected lists (array-index ordering):
 *   seg_by_riv[0] = [0, 1]
 *   seg_by_riv[1] = [2, 3]
 *   seg_by_ele[0] = [1, 3]
 *   seg_by_ele[1] = [0]
 *   seg_by_ele[2] = [2]
 *   upstream_by_down[1] = [0]      (riv 0's down = 2 ⇒ idown=1; toLake<=0)
 *   upstream_by_down[0] = []       (no one points down to riv 0)
 *   riv_in_by_lake[0] = [1]        (riv 1's toLake=0)
 *   ele_by_lake[0] = [1]           (ele 1's iLake=1 ⇒ bucket 0)
 *   lake_bank_edge_by_lake[0] = [(0,0), (2,1)]
 *   edge_by_ele[0] = [1, 2, -1]    (nabr 2,3,0 ⇒ 1, 2, -1)
 *   edge_by_ele[1] = [2, -1, 0]    (nabr 3,0,1 ⇒ 2, -1, 0)
 *   edge_by_ele[2] = [0, -1, -1]   (nabr 1,0,0 ⇒ 0, -1, -1)
 */
struct ExpectedSnapshot {
    bool rivseg_assert_should_be;
    bool riv_assert_should_be;
    bool ele_assert_should_be;
    std::vector<std::vector<int>> seg_by_riv_exp;
    std::vector<std::vector<int>> seg_by_ele_exp;
    std::vector<std::vector<int>> upstream_by_down_exp;
    std::vector<std::vector<int>> riv_in_by_lake_exp;
    std::vector<std::vector<int>> ele_by_lake_exp;
    std::vector<std::vector<std::pair<int,int>>> lake_bank_edge_by_lake_exp;
    std::vector<std::vector<int>> edge_by_ele_exp;
};

ExpectedSnapshot canonical_expected(){
    ExpectedSnapshot e;
    e.rivseg_assert_should_be = false;  /* test sets RivSeg[i].index = i+100 */
    e.riv_assert_should_be    = false;
    e.ele_assert_should_be    = false;
    e.seg_by_riv_exp = {{0, 1}, {2, 3}};
    e.seg_by_ele_exp = {{1, 3}, {0}, {2}};
    e.upstream_by_down_exp = {{}, {0}};
    e.riv_in_by_lake_exp = {{1}};
    e.ele_by_lake_exp = {{1}};
    e.lake_bank_edge_by_lake_exp = {{{0, 0}, {2, 1}}};
    e.edge_by_ele_exp = {{1, 2, -1}, {2, -1, 0}, {0, -1, -1}};
    return e;
}

Model_Data* make_mock_with_bad_indices(){
    Model_Data* mock = new Model_Data();
    mock->NumSegmt = 4;
    mock->NumRiv   = 2;
    mock->NumEle   = 3;
    mock->NumLake  = 1;

    mock->RivSeg = new RiverSegement[mock->NumSegmt];
    mock->Riv    = new _River[mock->NumRiv];
    mock->Ele    = new _Element[mock->NumEle];

    /* seg 0: iRiv=1, iEle=2; index OFFSET by +100 ⇒ fallback triggered */
    mock->RivSeg[0].index = 100; mock->RivSeg[0].iRiv = 1; mock->RivSeg[0].iEle = 2;
    mock->RivSeg[1].index = 101; mock->RivSeg[1].iRiv = 1; mock->RivSeg[1].iEle = 1;
    mock->RivSeg[2].index = 102; mock->RivSeg[2].iRiv = 2; mock->RivSeg[2].iEle = 3;
    mock->RivSeg[3].index = 103; mock->RivSeg[3].iRiv = 2; mock->RivSeg[3].iEle = 1;

    /* riv 0: down=2 (⇒ idown=1), toLake=-1 (no lake) */
    mock->Riv[0].index = 200; mock->Riv[0].down = 2;  mock->Riv[0].toLake = -1;
    /* riv 1: down=-1 (boundary), toLake=0 (lake 0) */
    mock->Riv[1].index = 201; mock->Riv[1].down = -1; mock->Riv[1].toLake = 0;

    /* ele 0: iLake=0, lakenabr={1,0,0}, nabr={2,3,0} */
    mock->Ele[0].index = 300; mock->Ele[0].iLake = 0;
    mock->Ele[0].lakenabr[0] = 1; mock->Ele[0].lakenabr[1] = 0; mock->Ele[0].lakenabr[2] = 0;
    mock->Ele[0].nabr[0] = 2; mock->Ele[0].nabr[1] = 3; mock->Ele[0].nabr[2] = 0;
    /* ele 1: iLake=1 (⇒ ele_by_lake bucket 0), lakenabr={0,0,0}, nabr={3,0,1} */
    mock->Ele[1].index = 301; mock->Ele[1].iLake = 1;
    mock->Ele[1].lakenabr[0] = 0; mock->Ele[1].lakenabr[1] = 0; mock->Ele[1].lakenabr[2] = 0;
    mock->Ele[1].nabr[0] = 3; mock->Ele[1].nabr[1] = 0; mock->Ele[1].nabr[2] = 1;
    /* ele 2: iLake=0, lakenabr={0,1,0}, nabr={1,0,0} */
    mock->Ele[2].index = 302; mock->Ele[2].iLake = 0;
    mock->Ele[2].lakenabr[0] = 0; mock->Ele[2].lakenabr[1] = 1; mock->Ele[2].lakenabr[2] = 0;
    mock->Ele[2].nabr[0] = 1; mock->Ele[2].nabr[1] = 0; mock->Ele[2].nabr[2] = 0;

    return mock;
}

#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s (file %s line %d)\n", msg, __FILE__, __LINE__); \
        return false; \
    } \
} while (0)

bool vectors_equal(const std::vector<int>& a, const std::vector<int>& b){
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return false;
    return true;
}

bool vectors_equal(const std::vector<std::pair<int,int>>& a,
                   const std::vector<std::pair<int,int>>& b){
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i){
        if (a[i].first != b[i].first || a[i].second != b[i].second) return false;
    }
    return true;
}

template <typename T>
bool lists_equal(const std::vector<std::vector<T>>& got,
                 const std::vector<std::vector<T>>& exp,
                 const char* name){
    if (got.size() != exp.size()){
        std::fprintf(stderr, "FAIL: %s outer size got=%zu exp=%zu\n",
                     name, got.size(), exp.size());
        return false;
    }
    for (size_t i = 0; i < exp.size(); ++i){
        if (!vectors_equal(got[i], exp[i])){
            std::fprintf(stderr, "FAIL: %s[%zu] differs\n", name, i);
            return false;
        }
    }
    return true;
}

bool run_fallback_test(){
    Model_Data* mock = make_mock_with_bad_indices();
    bool ok = build_adjacency_lists(mock);
    ExpectedSnapshot e = canonical_expected();

    /* Spec contract: when any entity violates `index == i + 1`, the
     * corresponding assert boolean flips to false AND the global
     * fallback flag is set. */
    CHECK(adjacency_assert_rivseg_pass == e.rivseg_assert_should_be,
          "adjacency_assert_rivseg_pass mismatch");
    CHECK(adjacency_assert_riv_pass == e.riv_assert_should_be,
          "adjacency_assert_riv_pass mismatch");
    CHECK(adjacency_assert_ele_pass == e.ele_assert_should_be,
          "adjacency_assert_ele_pass mismatch");
    CHECK(adjacency_fallback_triggered == true,
          "adjacency_fallback_triggered should be true on bad-index mock");
    CHECK(ok == false, "build_adjacency_lists return should be false on fallback");

    /* Lists must match the expected array-index-ordered output. The
     * fallback contract: when assert fails, build uses array index
     * order (NOT id-sort) — verified by checking the lists match
     * `for i in [0, N): if predicate(i) yield i` reference. */
    CHECK(lists_equal(seg_by_riv,        e.seg_by_riv_exp,        "seg_by_riv"),
          "seg_by_riv list does not match array-index ordering");
    CHECK(lists_equal(seg_by_ele,        e.seg_by_ele_exp,        "seg_by_ele"),
          "seg_by_ele list does not match array-index ordering");
    CHECK(lists_equal(upstream_by_down,  e.upstream_by_down_exp,  "upstream_by_down"),
          "upstream_by_down list does not match array-index ordering");
    CHECK(lists_equal(riv_in_by_lake,    e.riv_in_by_lake_exp,    "riv_in_by_lake"),
          "riv_in_by_lake list does not match array-index ordering");
    CHECK(lists_equal(ele_by_lake,       e.ele_by_lake_exp,       "ele_by_lake"),
          "ele_by_lake list does not match array-index ordering");
    CHECK(lists_equal(lake_bank_edge_by_lake, e.lake_bank_edge_by_lake_exp,
                     "lake_bank_edge_by_lake"),
          "lake_bank_edge_by_lake list does not match array-index ordering");
    CHECK(lists_equal(edge_by_ele,       e.edge_by_ele_exp,       "edge_by_ele"),
          "edge_by_ele list does not match array-index ordering");

    /* Intentionally do NOT delete mock: Model_Data's dtor calls
     * FreeData() which deletes 30+ unrelated heap arrays that the test
     * never allocated. Process exits cleanly on the assertion path. */
    return true;
}

}  /* anonymous namespace */

int main(){
    std::printf("=== test_adjacency_fallback (S4 PR-10 #154) ===\n");

    if (!run_fallback_test()){
        std::fprintf(stderr, "FAIL: fallback unit test\n");
        std::_Exit(1);
    }

    std::printf("PASS: adjacency_fallback_triggered + lists in array-index order\n");
    std::_Exit(0);
}
