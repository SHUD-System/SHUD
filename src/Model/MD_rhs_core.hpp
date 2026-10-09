/* MD_rhs_core.hpp — RHS core + ExecPolicy dispatch.
 *
 * `Model_Data::rhs_core(Y, DY, t, ExecPolicy)` evaluates the coupled
 * right-hand side as `rhs_update` -> `rhs_flux` -> `rhs_apply`. f() in
 * f.cpp is its only caller.
 *
 * ExecPolicy design:
 *   - Plain `enum class` + `switch` inside `rhs_core` — NO template
 *     specialization, NO virtual dispatch. Compile-time-known policy
 *     keeps the branch predictable; the switch is optimized away on
 *     fixed callers and matches SHUD's existing C-style enum style
 *     (Model_Control etc.).
 *   - Serial runs the three phases sequentially. StrictOMP runs them
 *     inside one OpenMP parallel region and gives results that are
 *     bit-identical to Serial for any thread count.
 *   - ProductionOMP is not implemented: it calls `std::abort()`.
 *     `assert(false)` must not be used instead, because `-DNDEBUG`
 *     strips assert to a no-op and would let execution silently fall
 *     through to the next statement.
 *   - Without SHUD_ENABLE_OPENMP_RHS the OMP cases are `#ifdef`ed out
 *     of the translation unit entirely; SHUD_ENABLE_OPENMP_RHS=1
 *     (the default for `make shud_omp`) compiles them in and f()
 *     selects StrictOMP.
 *
 * `rhs_update` / `rhs_flux` / `rhs_apply` / `rhs_core` are
 * `Model_Data::` members rather than free functions because their
 * bodies use the index macros `iSF` / `iUS` / `iGW` / `iRIV` / `iLAKE`
 * (Macros.hpp) which expand to expressions containing `NumEle` /
 * `NumRiv` — those resolve only with `this` in scope.
 *
 * Definitions live in MD_rhs_core.cpp. Member declarations also
 * appear on `Model_Data` (Model_Data.hpp) per C++ member-function
 * rules.
 *
 * Include this header only from `f.cpp` and `MD_rhs_core.cpp`, to
 * keep the compile surface from spreading.
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
