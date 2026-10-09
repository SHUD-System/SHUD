#ifndef MD_NVEC_HYBRID_HPP
#define MD_NVEC_HYBRID_HPP
/* =====================================================================
 * Hybrid NVector reduction overrides ("Config E", the default
 * `make shud_omp` build). Background: OpenMP_NVector_Determinism.md.
 *
 * Config E = parallel RHS + OpenMP NVector element-wise ops + SHUD-owned
 * SERIAL reduction overrides. The stock OpenMP NVector backend
 * (`SHUD_USE_OPENMP_NVECTOR=1`) parallelizes EVERY op, including the
 * reductions, whose `reduction(+:sum)` / `omp critical` partial-combine
 * order varies with thread count, so without these overrides ("Config D")
 * the results depend on the thread count.
 *
 * MECHANISM: after `N_VNew_OpenMP` creates the coupled
 * `udata`/`du`, overwrite the reduction entries of `v->ops` with
 * SHUD-owned functions written ONLY against the generic NVector API
 * (`N_VGetArrayPointer` + `N_VGetLength` + a plain serial loop). The
 * element-wise entries keep their stock OpenMP implementations (each
 * output element is a fixed per-element expression of same-index inputs
 * → identical FP order for any thread count AND vs the Serial backend).
 * `N_VClone`/`N_VCloneEmpty` copy the ops table (`N_VCopyOps`), so the
 * overrides propagate to every CVODE-internal temporary automatically.
 *
 * HARD RULES:
 *   - NEVER call an `N_V*_Serial` backend function on the OpenMP-content
 *     vector, and NEVER apply `NV_CONTENT_S`/`NV_DATA_S`/`NV_CONTENT_OMP`
 *     content-struct macros here: the OpenMP content struct carries a
 *     `num_threads` tail field, so prefix-layout aliasing is UB-adjacent.
 *     The generic `N_VGetArrayPointer`/`N_VGetLength` dispatch through the
 *     ops table and are correct for any backend.
 *   - Each override reproduces the corresponding `nvector_serial.c` body
 *     bit-for-bit (left-to-right accumulation, `min=x[0]` seed,
 *     `SUNSQR(x·w)` per element, `SUNRabs`), so Config E results equal
 *     those of the Serial NVector backend (Config C).
 *   - Install BEFORE `CVodeInit` AND BEFORE the `SHUD_NVEC_PROF` profiler
 *     wrap (overrides first, shims outermost — the profiler then delegates
 *     to the override pointer and reports `backend=hybrid`).
 *   - Standard reduction slots and their `*local` siblings are ALIASED to
 *     the same stock pointer; `install()` writes the serial override into
 *     BOTH slots for every reduction so no aliased sibling keeps the stock
 *     parallel body. The write is idempotent (a SHUD address is stored,
 *     the stock pointer is never read).
 *
 * SCOPE: only the coupled `udata`/`du` family is overridden (install() is
 * called once per vector at the coupled creation site). The decoupled
 * 5-solver loop vectors (SHUD_uncouple in shud.cpp) stay Serial and are
 * untouched. No code here lives on the RHS f() path. Compiled in ONLY
 * under `SHUD_NVEC_HYBRID` (the whole TU is `#ifdef`-guarded), so builds
 * without that flag are unaffected.
 * ===================================================================== */

#include <sundials/sundials_nvector.h>

/* Install SHUD-owned serial reduction overrides onto the ops table of
 * `v`, IN PLACE. Overrides every populated reduction/accumulating slot
 * (dotprod, maxnorm, wrmsnorm, wrmsnormmask, min, wl2norm, l1norm,
 * invtest, constrmask, minquotient, wsqrsumlocal, wsqrsummasklocal,
 * dotprodmultilocal + each aliased `*local` sibling); element-wise slots
 * are left stock. NULL slots (fused/vector-array, disabled by default)
 * stay NULL. Safe to call on udata AND du (idempotent — each call writes
 * the same SHUD addresses). Prints a one-line install summary to stdout.
 *
 * Config E2 (`SHUD_NVEC_DETRED=1`): the SUMMATION reduction slots
 * (dotprod / wsqrsum[mask] / wl2norm / l1norm + dotprodmultilocal + their
 * aliased `*local` siblings) are installed with FIXED-TREE deterministic
 * bodies instead of the plain serial fold; the non-summation reductions
 * (min / maxnorm / invtest / constrmask / minquotient) keep the plain
 * serial bodies (already cross-thread deterministic — no combine order to
 * fix). Selected entirely at compile time; with `SHUD_NVEC_DETRED` unset
 * this function is exactly the Config E install. */
void nvec_hybrid_install(N_Vector v);

/* Config E2 identity + parameters for the startup banner.
 * nvec_hybrid_detred_active() returns 1 iff built with `SHUD_NVEC_DETRED=1`
 * (fixed-tree deterministic reductions active), else 0 (plain Config E);
 * ..._block_size() returns the compile-time block size B; ..._neumaier()
 * returns 1 iff Neumaier compensation is compiled in. All three are
 * meaningful only under a hybrid build; the non-hybrid fallback returns
 * 0 / 0 / 0. */
int nvec_hybrid_detred_active(void);
int nvec_hybrid_detred_block_size(void);
int nvec_hybrid_detred_neumaier(void);

/* Smoke assert: clone `v` via BOTH N_VClone and N_VCloneEmpty and verify
 * each clone's 10 standard reduction op pointers equal the override
 * addresses installed on `v` (i.e. the ops-table copy carried the
 * overrides to the clone).
 * Prints a PASS/FAIL line to stdout and returns true on PASS. */
bool nvec_hybrid_clone_carries_overrides(N_Vector v);

/* Predicate: is `fp` one of the SHUD serial reduction override addresses
 * installed by nvec_hybrid_install()? Used by the PROF×HYBRID composition
 * assert: after the profiler wraps the (already-overridden) table, each
 * reduction shim's
 * captured delegate must equal the override address — the profiler calls
 * this predicate on each captured delegate. Returns false under a non-hybrid
 * build (the no-op fallback), so the composition assert is hybrid-only. */
bool nvec_hybrid_addr_is_override(void *fp);

#endif /* MD_NVEC_HYBRID_HPP */
