#ifndef MD_NVEC_HYBRID_HPP
#define MD_NVEC_HYBRID_HPP
/* =====================================================================
 * P12-nvec PR-N1 (#443) — Config E hybrid NVector reduction overrides.
 *
 * Config E = Config C StrictOMP RHS + OpenMP NVector element-wise ops +
 * SHUD-owned SERIAL reduction overrides. The OpenMP NVector backend
 * (`SHUD_USE_OPENMP_NVECTOR=1`) parallelizes EVERY op, including the
 * reductions, whose `reduction(+:sum)` / `omp critical` partial-combine
 * order varies with thread count → the Config D trajectory drift.
 *
 * MECHANISM (design D2): after `N_VNew_OpenMP` creates the coupled
 * `udata`/`du`, overwrite the reduction entries of `v->ops` with
 * SHUD-owned functions written ONLY against the generic NVector API
 * (`N_VGetArrayPointer` + `N_VGetLength` + a plain serial loop). The
 * element-wise entries keep their stock OpenMP implementations (each
 * output element is a fixed per-element expression of same-index inputs
 * → identical FP order for any thread count AND vs the Serial backend).
 * `N_VClone`/`N_VCloneEmpty` copy the ops table (`N_VCopyOps`), so the
 * overrides propagate to every CVODE-internal temporary automatically.
 *
 * HARD RULES (design D2 §Hard rules, spec "generic-API serial override"):
 *   - NEVER call an `N_V*_Serial` backend function on the OpenMP-content
 *     vector, and NEVER apply `NV_CONTENT_S`/`NV_DATA_S`/`NV_CONTENT_OMP`
 *     content-struct macros here: the OpenMP content struct carries a
 *     `num_threads` tail field, so prefix-layout aliasing is UB-adjacent.
 *     The generic `N_VGetArrayPointer`/`N_VGetLength` dispatch through the
 *     ops table and are correct for any backend.
 *   - Each override reproduces the corresponding `nvector_serial.c` body
 *     bit-for-bit (left-to-right accumulation, `min=x[0]` seed,
 *     `SUNSQR(x·w)` per element, `SUNRabs`), so Config E == Serial == C.
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
 * 5-solver loop vectors (shud.cpp:495+) stay Serial and are untouched.
 * No code here lives on the RHS f() path; f.cpp / MD_rhs_core.cpp are
 * byte-unchanged. Compiled in ONLY under `SHUD_NVEC_HYBRID` (the whole TU
 * is `#ifdef`-guarded), so default builds are preprocessor-identical.
 * ===================================================================== */

#include <sundials/sundials_nvector.h>

/* Install SHUD-owned serial reduction overrides onto the ops table of
 * `v`, IN PLACE. Overrides every populated reduction/accumulating slot
 * (dotprod, maxnorm, wrmsnorm, wrmsnormmask, min, wl2norm, l1norm,
 * invtest, constrmask, minquotient, wsqrsumlocal, wsqrsummasklocal,
 * dotprodmultilocal + each aliased `*local` sibling); element-wise slots
 * are left stock. NULL slots (fused/vector-array, disabled by default)
 * stay NULL. Safe to call on udata AND du (idempotent — each call writes
 * the same SHUD addresses). Prints a one-line install summary to stdout
 * for the evidence log. */
void nvec_hybrid_install(N_Vector v);

/* Smoke assert (spec scenario "propagation assert"): clone `v` via BOTH
 * N_VClone and N_VCloneEmpty and verify each clone's reduction op pointers
 * equal the override addresses installed on `v` (i.e. the ops-table copy
 * carried the overrides to the clone), and DIFFER from the stock OpenMP
 * addresses captured from a fresh stock vector of the same length/threads.
 * Prints a PASS/FAIL line to stdout and returns true on PASS. */
bool nvec_hybrid_clone_carries_overrides(N_Vector v);

/* Predicate: is `fp` one of the SHUD serial reduction override addresses
 * installed by nvec_hybrid_install()? Used by the PROF×HYBRID composition
 * assert (spec nvec-op-profile "composition with hybrid overrides"): after
 * the profiler wraps the (already-overridden) table, each reduction shim's
 * captured delegate must equal the override address — the profiler calls
 * this predicate on each captured delegate. Returns false under a non-hybrid
 * build (the no-op fallback), so the composition assert is hybrid-only. */
bool nvec_hybrid_addr_is_override(void *fp);

#endif /* MD_NVEC_HYBRID_HPP */
