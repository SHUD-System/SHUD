/* sunlinsol_hypre.h — P8-tune.G0 PR-0 (openspec change
 *   p8tune-g0-instrumented-amg-smoke). Custom SUNDIALS 6.0
 *   `SUNLinearSolver` wrapper around HYPRE's BoomerAMG.
 *
 * Purpose:
 *   Provide a matrix-free iterative LS plug-in for CVODE 6.0 that
 *   uses BoomerAMG as the inner solver. The wrapper exposes the full
 *   15-callback `SUNLinearSolver_Ops` ABI required by SUNDIALS 6.0
 *   (`SHUD/InstallSundials/include/sundials/sundials_linearsolver.h`
 *   L108-127) with NO NULL slots.
 *
 * Activation: opt-in via `SHUD_LINSOL=amg` environment variable. When
 *   the env var is unset OR set to "spgmr" (default), the
 *   `cvode_config.cpp` factory dispatches `SUNLinSol_SPGMR(...)`
 *   instead and this wrapper is never instantiated. The default-path
 *   bit-identical SPGMR baseline (G0-1 anchor) is therefore
 *   preserved.
 *
 * G0 scope:
 *   - Hardcoded `(interp_type=6, coarsen_type=8)` per spec REQ
 *     "Reject-non-(6,8) constructor guard". G1 may relax.
 *   - Hypre threads pinned to 1 under `shud_omp` (`HYPRE_SetGlobalOptions
 *     "default_thread_count=1"`) — defense-in-depth alongside the
 *     sbatch script `export OMP_NUM_THREADS=1` (design.md D10).
 *   - Setup-call cadence: PR-0 spike (task 1.4) measured CVODE
 *     `nsetups=0` on keliya 90-day SHORT under SPGMR `PREC_NONE`.
 *     The wrapper therefore lazy-builds its AMG hierarchy at the
 *     FIRST `Solve` invocation rather than relying on CVODE-issued
 *     Setup events.
 *
 * Telemetry:
 *   - Ring buffer of per-Solve entries
 *     (`hypre_iters`, `hypre_op_count`, `setup_wall_sec`,
 *      `solve_wall_sec`, plus per-step context populated via
 *      `SUNLinSol_Hypre_SetStepContext`).
 *   - `MARKER:AMG_TELEMETRY_REAL` emitted to stdout on first
 *     successful Solve that captures both `Get*` calls
 *     (idempotent — single emission per process lifetime).
 *   - `SUNLinSol_Hypre_DrainTelemetry(LS, fp)` writes TSV rows.
 *
 * Failure semantics (driver-visible):
 *   - `(interp,coarsen) != (6,8)` constructor → returns NULL with
 *     stderr error.
 *   - NULL `y` constructor → returns NULL with stderr error.
 *   - HYPRE missing at runtime → factory in `cvode_config.cpp`
 *     fatal-exits BEFORE wrapper allocation.
 *   - BoomerAMG divergence in Solve → returns `SUNLS_CONV_FAIL`
 *     (recoverable for CVODE retry) and writes
 *     `MARKER:AMG_SOLVE_DIVERGE_DETECTED` to stderr.
 *   - Structural Hypre failure → returns `SUNLS_PACKAGE_FAIL_UNREC`.
 */

#ifndef SUNLINSOL_HYPRE_H
#define SUNLINSOL_HYPRE_H

#include <stdio.h>
#include "sundials/sundials_types.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_linearsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Forward decl of SHUD Model_Data so the header is includable from
 * pure-C consumers without dragging the C++ class definition in.
 * Setup obtains topology via the void* (cast to Model_Data*)
 * stashed at construction time. */
struct Model_Data_fwd;

/* Constructor — must be called with the full N_Vector `y` (matching
 * the cvode_config.cpp SetCVODE site, where `udata` is itself an
 * N_Vector). `MD` is the SHUD topology / state context (passed
 * through via void* so this header doesn't need <Model_Data.hpp>).
 *
 * Hardcoded `(interp_type=6, coarsen_type=8)` per G0 spec; any other
 * pair returns NULL with stderr error (constructor guard).
 *
 * Returns NULL on:
 *   - `y == NULL`
 *   - `(interp_type, coarsen_type) != (6, 8)`
 *   - `SUNLinSolNewEmpty(sunctx)` failure (rare)
 *   - content struct allocation failure
 *
 * On success, the returned `SUNLinearSolver` has all 15
 * `SUNLinearSolver_Ops` slots populated (build-time assertion). */
SUNLinearSolver SUNLinSol_Hypre(N_Vector y,
                                void *MD,
                                int interp_type,
                                int coarsen_type,
                                SUNContext sunctx);

/* Per-step context plumbing (driver-side telemetry path per PR-0
 * spike task 1.5). SHUD driver loop calls this BEFORE each
 * `CVode(...)` invocation to populate the next ring-buffer entry's
 * per-step context fields. Computes deltas internally against the
 * previous call's cumulative counters. No-op + returns SUNLS_SUCCESS
 * when `LS` is NULL or not a Hypre wrapper.
 *
 * Parameters:
 *   `step_idx`        — monotone step index (caller-defined).
 *   `t_sim`           — simulation time at this step.
 *   `cvode_nli_cum`   — `CVodeGetNumLinIters` cumulative at call time.
 *   `cvode_nfeLS_cum` — `CVodeGetNumLinRhsEvals` cumulative at call. */
int SUNLinSol_Hypre_SetStepContext(SUNLinearSolver LS,
                                   long step_idx,
                                   realtype t_sim,
                                   long cvode_nli_cum,
                                   long cvode_nfeLS_cum);

/* Drain the per-Solve telemetry ring buffer to `out` as TSV rows.
 * Schema (one row per Solve):
 *   step_idx  t_sim  setup_called  hypre_iters  hypre_op_count
 *   setup_wall_sec  solve_wall_sec  cvode_nli_step  cvode_nfeLS_step
 * Resets head/tail/overflow-counter after draining. Returns the
 * number of rows written (≥0) or a negative SUNLS_* return code on
 * error. No-op + returns 0 when `LS` is NULL or out is NULL. */
int SUNLinSol_Hypre_DrainTelemetry(SUNLinearSolver LS, FILE *out);

#ifdef __cplusplus
}
#endif

#endif  /* SUNLINSOL_HYPRE_H */
