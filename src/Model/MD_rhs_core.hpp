/* MD_rhs_core.hpp — RHS core scaffolding + ExecPolicy dispatch.
 *
 * Stage history:
 *   - S1a (openMP #44): introduced `Model_Data::rhs_update` (pure
 *     carry-over of `f_update`) and the three-arg `rhs_core(Y, DY, t)`
 *     dispatch scaffold, gated by an S1a scaffold Makefile flag.
 *   - S1b (openMP #45): added `Model_Data::rhs_flux` (carry-over of
 *     `f_loop`); rhs_core became `rhs_update + rhs_flux + legacy f_applyDY`.
 *   - S1c (openMP #46): added `Model_Data::rhs_apply` (carry-over of
 *     `f_applyDY`); rhs_core became zero-fallback `rhs_update + rhs_flux
 *     + rhs_apply`. The scaffold flag still gated f.cpp dispatch.
 *   - S1d.1 (openMP #47): replaced the three-arg signature with the
 *     four-arg ExecPolicy form `rhs_core(Y, DY, t, ExecPolicy)`. The
 *     S1a scaffold Makefile / source flag is RETIRED in the same
 *     atomic commit.
 *   - S2 capstone (PR-8): the legacy-vs-rhs_core gating macro that
 *     selected between the original B0 path (`f_update/f_loop/
 *     f_applyDY`) and the B1a path (`rhs_core(..., ExecPolicy::Serial)`)
 *     has been retired; f.cpp now unconditionally calls
 *     `rhs_core(..., ExecPolicy::Serial)`.
 *
 * ExecPolicy design (per openspec/changes/s1-rhs-core-extraction/
 * design.md D7 + spec exec-policy-enum):
 *   - Plain `enum class` + `switch` inside `rhs_core` — NO template
 *     specialization, NO virtual dispatch. Compile-time-known policy
 *     keeps the branch predictable; the switch is optimized away on
 *     fixed callers and matches SHUD's existing C-style enum style
 *     (Model_Control etc.).
 *   - StrictOMP / ProductionOMP cases are S1-phase stubs: each calls
 *     `std::abort()` immediately. `assert(false)` is forbidden because
 *     `-DNDEBUG` (release builds, `EXTRA_CXXFLAGS=-DNDEBUG` smoke
 *     compile) strips assert to a no-op and would let execution
 *     silently fall through to the next statement — destroying the
 *     contract that OMP paths cannot impersonate Serial.
 *   - SHUD_ENABLE_OPENMP_RHS=0 (default) `#ifdef`s the OMP cases out
 *     of the translation unit entirely; the binary contains no OMP
 *     path symbols. SHUD_ENABLE_OPENMP_RHS=1 includes the cases for
 *     smoke compile + runtime SIGABRT verification.
 *
 * `rhs_update` / `rhs_flux` / `rhs_apply` / `rhs_core` are
 * `Model_Data::` members rather than free functions because the legacy
 * function bodies use the index macros `iSF` / `iUS` / `iGW` / `iRIV`
 * / `iLAKE` (Macros.hpp:21-25) which expand to expressions containing
 * `NumEle` / `NumRiv` — those resolve only with `this` in scope.
 *
 * Definitions live in MD_rhs_core.cpp. Member declarations also
 * appear on `Model_Data` (Model_Data.hpp) per C++ member-function
 * rules.
 *
 * Header is included ONLY by `SHUD/src/Model/f.cpp` and
 * `SHUD/src/Model/MD_rhs_core.cpp` per spec Scenario
 * "Header is included only by `f.cpp` and `MD_rhs_core.cpp`" — keeps
 * compile-surface from spreading.
 */
#ifndef MD_RHS_CORE_HPP
#define MD_RHS_CORE_HPP

/* `enum class ExecPolicy` is defined in Model_Data.hpp (see comment
 * there) to break the circular-dependency that would arise if it
 * lived here — Model_Data::rhs_core takes it by value, and
 * MD_rhs_core.hpp #include's Model_Data.hpp. Including Model_Data.hpp
 * here re-exposes the enum for f.cpp / MD_rhs_core.cpp consumers. */
#include "Model_Data.hpp"

#endif /* MD_RHS_CORE_HPP */
