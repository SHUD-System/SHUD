#ifndef MD_NVEC_PROF_HPP
#define MD_NVEC_PROF_HPP
/* =====================================================================
 * P12-nvec PR-N0 (#442) — env-gated NVector op-share profiler.
 *
 * SHUD_NVEC_PROF=1 (STRICT `=1` string compare, P11-osc predicate
 * discipline) → wrap every POPULATED function-pointer of the coupled-
 * solver vector's `ops` table at creation time (before CVodeInit) with a
 * shim that (a) increments a per-op call counter, (b) accumulates
 * monotonic-clock elapsed nanoseconds, then (c) delegates to the ORIGINAL
 * implementation — pure delegation, zero reordering of any computation.
 *
 * WHY an ops-table wrapper (design D1): existing profile buckets stop at
 * t_CVODE_raw (the whole CVode() call) and cannot see inside. The wrapper
 * measures EXACTLY the NVector boundary we would parallelize, and (unlike
 * a sampling profiler) needs no symbol attribution through the SUNDIALS
 * static libs. The dump decomposes t_CVODE_raw into element-wise NVector
 * work vs reduction NVector work vs the non-NVector remainder — the
 * measured input to gates G-E2 context and G-E3(ii)/(iii).
 *
 * COMPOSITION ORDER (spec "ops-table wrapper instrumentation" + design D1
 * / D2): the profiler wraps LAST / OUTERMOST. When the PR-N1 Config E
 * hybrid reduction overrides (SHUD_NVEC_HYBRID=1) are also active they
 * install FIRST, so each shim delegates to the EFFECTIVE (overridden)
 * pointer and no reduction measurement is clobbered. PR-N0 ships no
 * hybrid path yet; install() is written so the override step (PR-N1) can
 * run against the same vector BEFORE this call with no change here.
 *
 * CLONE PROPAGATION: wrapping happens ONCE at creation, before CVodeInit.
 * The generic N_VClone / N_VCloneEmpty copy the ops table (N_VCopyOps),
 * so every CVODE-internal temporary inherits the shim table automatically
 * → all internal vector work is measured. nvec_prof_clone_carries_shims()
 * is the smoke assert (spec scenario "clone propagation").
 *
 * SCOPE: only the coupled `udata` / `du` family is wrapped (the single
 * coupled-solver vector family). The decoupled 5-solver loop vectors
 * (shud.cpp:495+, N_VNew_Serial u1..u5 / du1..du5) are NOT wrapped —
 * install() is called only on the coupled creation site.
 *
 * BITWISE NEUTRALITY: gate OFF (unset / "" / "0" / any non-"1") → shims
 * are NEVER installed and the creation path is byte-identical to the
 * unpatched build. Gate ON → shims reorder nothing (pure delegation), so
 * trajectories stay identical; a keliya SHA leg proves it regardless
 * (P11-osc style). No code here lives on the RHS f() call path; f.cpp /
 * MD_rhs_core.cpp are untouched.
 * ===================================================================== */

#include <sundials/sundials_nvector.h>

/* Strict `=1` env predicate — the single gate semantics. getenv present
 * AND exact string "1". NULL, "", "0", "true", "1 " etc. all → false.
 * Mirrors osc_diag_env_is_1 (MD_osc_diag.hpp). */
bool nvec_prof_env_is_1(const char *name);

/* True iff SHUD_NVEC_PROF=1 at the moment of the call (read once and
 * cached inside install()). Callers gate the dump on the SAME state. */
bool nvec_prof_is_on();

/* Install counting/timing shims onto the ops table of `v`, IN PLACE.
 * No-op (and `v` untouched) unless SHUD_NVEC_PROF=1. Only POPULATED
 * (non-NULL) computational op slots are wrapped; NULL slots (e.g. the
 * fused/array ops CVODE leaves disabled by default) stay NULL, so there
 * is no unwrapped-op leakage and no NULL dispatch. Pure-utility slots
 * that are not computational work (nvgetarraypointer, nvsetarraypointer,
 * nvgetlength, nvspace, nvgetvectorid, nvgetcommunicator,
 * nvgetdevicearraypointer, nvdestroy, XBraid buf ops, debug print) are
 * intentionally left un-wrapped — they carry no FP work to attribute and
 * (nvgetarraypointer) sit on unrelated hot paths.
 *
 * `backend` is one of "serial" / "openmp" / "hybrid" and is recorded in
 * the CSV header. Safe to call on multiple vectors of the same family
 * (udata AND du): the shim table is a per-process singleton keyed by op
 * slot, so wrapping du after udata is idempotent (the original pointers
 * are identical and already captured). */
void nvec_prof_install(N_Vector v, const char *backend);

/* Smoke assert (spec scenario "clone propagation"): clone `v` via BOTH
 * N_VClone and N_VCloneEmpty and verify each clone's wrapped op pointers
 * equal the shims installed on `v` (i.e. the ops table copy carried the
 * shims to the clone). No-op unless the profiler is on. Prints a PASS/
 * FAIL line to stdout and returns true on PASS. */
bool nvec_prof_clone_carries_shims(N_Vector v);

/* Dump nvec_prof.csv to `outpath` (the project output dir). No-op unless
 * the profiler is on. Header lines carry project_name / NY / nthreads /
 * backend; one data row per wrapped op: op_name,op_class,calls,total_ns
 * with op_class ∈ {elementwise, reduction, other} from the fixed source-
 * committed mapping table below. Also emits, to stdout, the wrapped-entry
 * debug count for the spec's "no unwrapped-op leakage" cross-check. */
void nvec_prof_dump(const char *project_name, int NY, int nthreads,
                    const char *backend, const char *outpath);

#endif /* MD_NVEC_PROF_HPP */
