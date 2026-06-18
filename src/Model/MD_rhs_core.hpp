/* MD_rhs_core.hpp — S1a RHS core scaffolding (openMP issue #44).
 *
 * Stage: S1a (Group 1, tasks 1.1–1.10) — extract `Model_Data::f_update`
 * → `Model_Data::rhs_update` and introduce `Model_Data::rhs_core`
 * dispatch scaffold. Per design.md D7 + spec rhs-core-scaffolding
 * ADDED Requirement "S1a 阶段签名约定", `rhs_core()` is the
 * THREE-PARAM data-flow form `(Y, DY, t)` — no `ExecPolicy` parameter.
 * The four-arg `ExecPolicy` overload is introduced in S1d.1 by
 * capability exec-policy-enum.
 *
 * `rhs_update` and `rhs_core` are `Model_Data::` members rather than
 * free functions because:
 *   - the legacy `Model_Data::f_update` body uses the index macros
 *     `iSF` / `iUS` / `iGW` / `iRIV` / `iLAKE` (Macros.hpp:21-25)
 *     which expand to expressions containing `NumEle` / `NumRiv` —
 *     those are `Model_Data` members and only resolve correctly with
 *     `this` in scope;
 *   - spec Scenario "Source carry-over diff is structural-only"
 *     explicitly allows `Model_Data::` qualifier in the diff.
 *
 * Definitions live in MD_rhs_core.cpp. Declarations also appear on
 * `Model_Data` (Model_Data.hpp) per C++ member-function rules.
 *
 * Header is included ONLY by `SHUD/src/Model/f.cpp` and
 * `SHUD/src/Model/MD_rhs_core.cpp` per spec Scenario
 * "Header is included only by `f.cpp` and `MD_rhs_core.cpp`" — keeps
 * compile-surface from spreading.
 */
#ifndef MD_RHS_CORE_HPP
#define MD_RHS_CORE_HPP

#include "Model_Data.hpp"

#endif /* MD_RHS_CORE_HPP */
