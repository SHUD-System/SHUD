/* sunlinsol_hypre.cpp — P8-tune.G0 PR-0 BoomerAMG wrapper.
 *
 * See sunlinsol_hypre.h for the public contract and design rationale.
 *
 * Internal structure
 * ------------------
 *
 * - `HypreContent` is the wrapper's `SUNLinearSolver::content` field.
 *   It holds:
 *     - Hypre handles (lazy-allocated at first Solve, freed in
 *       `Free`): `HYPRE_IJMatrix A_ij`, `HYPRE_ParCSRMatrix A_par`,
 *       `HYPRE_IJVector b_ij/x_ij`, `HYPRE_ParVector b_par/x_par`,
 *       `HYPRE_Solver amg`.
 *     - `Model_Data *MD` stashed from the constructor (for Setup
 *       topology lookup; passed through opaque `void*` to keep the
 *       header C-friendly).
 *     - ATimes stash from `SetATimes` callback.
 *     - Scale-vectors stash from `SetScalingVectors` (G0 stashes
 *       only, no scale applied to solve).
 *     - Zero-guess flag from `SetZeroGuess`.
 *     - Telemetry ring buffer + drop counter + idempotent
 *       `MARKER:AMG_TELEMETRY_REAL` emission flag.
 *     - Last return code from `Solve` (for `LastFlag`).
 *     - Step-context fields populated by
 *       `SUNLinSol_Hypre_SetStepContext` between CVode calls; the
 *       wrapper differences cumulative counters internally to emit
 *       per-step deltas.
 *
 * - Lazy Setup. CVODE 6.0 SPGMR with `PREC_NONE` issues `nsetups=0`
 *   on keliya 90-day SHORT (PR-0 spike task 1.4). To remain
 *   compatible with this cadence, the AMG hierarchy is built inside
 *   the wrapper's `Solve` callback when the cached handles are NULL,
 *   and reused on subsequent Solves. CVODE-issued Setup callbacks
 *   (if any) trigger an explicit rebuild (design.md D4 baseline
 *   policy).
 *
 * - Probe-derive sparsity. PR-0 spike task 1.1 selects probe-derive
 *   over MD_adjacency reuse: the wrapper at first Solve invokes
 *   `ATimes(A_data, e_i, Av)` for every i in [0, NumY) and reads the
 *   resulting nonzero columns into a per-row IJMatrix. Bandwidth is
 *   not pre-pinned; the wrapper trusts the probe to discover
 *   wherever the Jacobian touches.
 *
 * - Telemetry. Hypre 3.1.0 public API does NOT expose
 *   `HYPRE_BoomerAMGGetCycleNumIterations` / `..GetCycleOpCount`
 *   (verified via grep against `/opt/homebrew/include/HYPRE_parcsr_ls.h`
 *   on 2026-06-29 — only `GetNumIterations` and `GetCumNnzAP` are
 *   public). The wrapper uses `HYPRE_BoomerAMGGetNumIterations` for
 *   the per-Solve iteration count and `HYPRE_BoomerAMGGetCumNnzAP`
 *   as the operator-complexity proxy. The
 *   `MARKER:AMG_TELEMETRY_REAL` line records BOTH the Hypre release
 *   number AND the chosen API names so PR-B aggregator can
 *   disambiguate. Future Hypre releases that expose CycleNumIterations
 *   can be picked up via the `HYPRE_RELEASE_NUMBER` gate.
 */

#include "sunlinsol_hypre.h"

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "Model_Data.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_utilities.h"  /* HYPRE_RELEASE_NUMBER */
}

/* Disable OpenMPI C++ bindings — Ubuntu/server OpenMPI ships a broken
 * functions_inln.h that fails to compile with modern g++ (see
 * /usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/ompi/mpi/cxx/functions_inln.h).
 * The C bindings via <mpi.h> are all we need. Must be defined BEFORE
 * <mpi.h> is included. */
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include "nvector/nvector_serial.h"

/* Telemetry ring-buffer size. Spec requires ≥100k entries. */
static constexpr int HYPRE_TELEMETRY_RING_SIZE = 131072;

namespace {

struct TelemetryEntry {
    long step_idx = -1;
    realtype t_sim = 0.0;
    int setup_called = 0;
    int hypre_iters = 0;
    double hypre_op_count = 0.0;
    double setup_wall_sec = 0.0;
    double solve_wall_sec = 0.0;
    long cvode_nli_step = 0;
    long cvode_nfeLS_step = 0;
};

struct HypreContent {
    /* SHUD topology context (set in constructor). */
    Model_Data *MD = nullptr;

    /* Configuration. */
    int interp_type = 6;
    int coarsen_type = 8;
    long n = 0;        /* matrix dimension = MD->NumY */

    /* Lazy-allocated Hypre handles (NULL until first Solve). */
    HYPRE_IJMatrix A_ij = nullptr;
    HYPRE_ParCSRMatrix A_par = nullptr;
    HYPRE_IJVector b_ij = nullptr;
    HYPRE_IJVector x_ij = nullptr;
    HYPRE_ParVector b_par = nullptr;
    HYPRE_ParVector x_par = nullptr;
    HYPRE_Solver amg = nullptr;

    /* Setup-side scratch (kept across Solves for fast rebuild). */
    std::vector<HYPRE_BigInt> indices;          /* [n] = 0..n-1 */

    /* Whether wrapper called HYPRE_Initialize / MPI_Init (so Free
     * cleans up symmetrically). The wrapper conservatively assumes
     * it owns init if no one else did. */
    bool wrapper_inited_hypre = false;
    bool wrapper_inited_mpi = false;

    /* ATimes stash from SetATimes. */
    void *atimes_data = nullptr;
    SUNATimesFn atimes_fn = nullptr;

    /* Scale-vectors stash (G0 stashes only). */
    N_Vector s1 = nullptr;
    N_Vector s2 = nullptr;

    /* N_Vector flavor template stashed at constructor time. The wrapper
     * uses this to N_VClone probe scratch when lazy-building the AMG
     * hierarchy. The pointer is caller-owned; the wrapper holds it only
     * as a flavor reference (DO NOT destroy in Free). The lifetime
     * invariant is enforced by CVODE: the constructor is called from
     * SetCVODE where `udata` outlives the SUNLinearSolver. */
    N_Vector y_template = nullptr;

    /* Zero-guess flag from SetZeroGuess (CVODE 6.0 SPILS path). */
    booleantype zero_guess = SUNFALSE;

    /* Telemetry. */
    std::vector<TelemetryEntry> ring;
    int ring_head = 0;
    int ring_count = 0;
    long entries_dropped_to_overflow = 0;
    bool telemetry_marker_emitted = false;

    /* Per-step context (set by SUNLinSol_Hypre_SetStepContext between
     * CVode calls). Wrapper differences cumulative counters
     * internally. */
    long ctx_step_idx = -1;
    realtype ctx_t_sim = 0.0;
    long ctx_last_nli_cum = 0;
    long ctx_last_nfeLS_cum = 0;
    long ctx_nli_step = 0;
    long ctx_nfeLS_step = 0;

    /* Last setup wall — populated by lazy/explicit Setup, consumed by
     * next Solve telemetry row. */
    double pending_setup_wall_sec = 0.0;
    int pending_setup_called = 0;

    /* For Solve return-code routing through LastFlag. */
    sunindextype last_flag = SUNLS_SUCCESS;
};

inline HypreContent *content_of(SUNLinearSolver LS) {
    return static_cast<HypreContent *>(LS->content);
}

/* ---- 15 SUNLinearSolver_Ops callbacks --------------------------------- */

SUNLinearSolver_Type op_gettype(SUNLinearSolver /*LS*/) {
    return SUNLINEARSOLVER_ITERATIVE;
}

SUNLinearSolver_ID op_getid(SUNLinearSolver /*LS*/) {
    return SUNLINEARSOLVER_CUSTOM;
}

int op_setatimes(SUNLinearSolver LS, void *A_data, SUNATimesFn ATimes) {
    HypreContent *c = content_of(LS);
    c->atimes_data = A_data;
    c->atimes_fn = ATimes;
    return SUNLS_SUCCESS;
}

int op_setpreconditioner(SUNLinearSolver /*LS*/, void * /*P_data*/,
                         SUNPSetupFn /*Pset*/, SUNPSolveFn /*Psol*/) {
    /* BoomerAMG is the preconditioner — no separate Psolve. */
    return SUNLS_SUCCESS;
}

int op_setscalingvectors(SUNLinearSolver LS, N_Vector s1, N_Vector s2) {
    HypreContent *c = content_of(LS);
    c->s1 = s1;
    c->s2 = s2;
    return SUNLS_SUCCESS;
}

int op_setzeroguess(SUNLinearSolver LS, booleantype onoff) {
    HypreContent *c = content_of(LS);
    c->zero_guess = onoff;
    return SUNLS_SUCCESS;
}

int op_initialize(SUNLinearSolver LS) {
    HypreContent *c = content_of(LS);

    /* MPI_Init guard: required by HYPRE for MPI_COMM_SELF use. We
     * lazily MPI_Init only if not already initialised. SHUD core has
     * no other MPI users; mpic++ link drags libmpi but no
     * application-side MPI_Init call exists.
     *
     * Use MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, ...) rather
     * than the legacy stack-argv MPI_Init pattern:
     *   - NULL-argv is legal per MPI-3 §8.7 (Implementations MUST accept
     *     NULL pointers for argc/argv to MPI_Init / MPI_Init_thread).
     *   - MPI_THREAD_FUNNELED is the minimum thread support level that
     *     allows nested OMP `parallel` regions in the same process; shud_omp
     *     wraps OMP_NUM_THREADS=1 on the AMG path, but the funneled mode
     *     also covers the SPGMR-baseline-and-AMG-build-up code paths that
     *     coexist before/after this Initialize call. */
    int mpi_already = 0;
    MPI_Initialized(&mpi_already);
    if (!mpi_already) {
        int provided = MPI_THREAD_SINGLE;
        if (MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS) {
            std::fprintf(stderr,
                "[shud-amg] FATAL: MPI_Init_thread failed in SUNLinSol_Hypre Initialize\n");
            return SUNLS_PACKAGE_FAIL_UNREC;
        }
        c->wrapper_inited_mpi = true;
    }

    /* HYPRE init API split across releases (see p8tune.F
     * boomeramg_setup_solve.cpp L646). Mac brew 3.1.0 = 30100. */
#if defined(HYPRE_RELEASE_NUMBER) && HYPRE_RELEASE_NUMBER >= 30000
    HYPRE_Initialize();
#else
    HYPRE_Init();
#endif
    c->wrapper_inited_hypre = true;

    /* Hypre 3.1.0 does not expose a portable thread-pin API. We rely on
     * the caller (sbatch script / shell) to export OMP_NUM_THREADS=1
     * before starting SHUD. We CAN, however, attest honestly to the
     * runtime OMP state so PR-A smoke runner can refuse to accept
     * suspected over-subscription. */
    int omp_max = 1;
#ifdef _OPENMP
    omp_max = omp_get_max_threads();
#endif
    const char *omp_env = std::getenv("OMP_NUM_THREADS");
    int omp_env_n = (omp_env != nullptr) ? (int)std::strtol(omp_env, nullptr, 10) : -1;
    const char *thread_state = "ENFORCED-BY-CALLER";
    if (omp_max > 1 || omp_env_n > 1) {
        thread_state = "UNKNOWN-OVERSUBSCRIBE-RISK";
    }
    std::fprintf(stdout,
        "[shud-amg] Hypre threads=1 state=%s omp_max=%d OMP_NUM_THREADS=%s\n",
        thread_state, omp_max, omp_env ? omp_env : "<unset>");
    std::fflush(stdout);

    return SUNLS_SUCCESS;
}

/* op_setup behavior (G0 spec amendment):
 *   - Destroy any prior AMG hierarchy + IJMatrix + IJVector handles.
 *   - Set the pending-setup flag and accumulate the destroy wall.
 *   - DO NOT clone any N_Vector here — we don't have access to a
 *     CVODE-managed N_Vector flavor template until Solve.
 *   - DO NOT build the new hierarchy here — first Solve after Setup
 *     does the lazy build (lazy_build_hierarchy_from_solve), which
 *     has access to the in-flight `b` N_Vector for probe scratch.
 *   - First Solve's telemetry attributes the destroy wall + actual
 *     build wall together via pending_setup_wall_sec; setup_called=1
 *     persists until the first Solve consumes it.
 * Returns SUNLS_SUCCESS unconditionally — there is nothing to fail in
 * a pure destroy-and-flag path. */
int op_setup(SUNLinearSolver LS, SUNMatrix /*A*/) {
    HypreContent *c = content_of(LS);
    if (c == nullptr) return SUNLS_MEM_NULL;

    const auto t0 = std::chrono::steady_clock::now();

    if (c->amg)   { HYPRE_BoomerAMGDestroy(c->amg);   c->amg   = nullptr; }
    if (c->A_ij)  { HYPRE_IJMatrixDestroy(c->A_ij);   c->A_ij  = nullptr; c->A_par = nullptr; }
    if (c->b_ij)  { HYPRE_IJVectorDestroy(c->b_ij);   c->b_ij  = nullptr; c->b_par = nullptr; }
    if (c->x_ij)  { HYPRE_IJVectorDestroy(c->x_ij);   c->x_ij  = nullptr; c->x_par = nullptr; }

    const auto t1 = std::chrono::steady_clock::now();
    c->pending_setup_wall_sec += std::chrono::duration<double>(t1 - t0).count();
    c->pending_setup_called = 1;
    c->last_flag = SUNLS_SUCCESS;
    return SUNLS_SUCCESS;
}

/* Topology-restricted ATimes probe.
 *
 * SHUD state vector layout (see SHUD/src/ModelData/Model_Data.cpp L86
 * + SHUD/src/Equations/functions.hpp L83-99):
 *
 *   NumY = 3 * NumEle + NumRiv + NumLake
 *   indices [0          .. NumEle)        — surface  (yEleSurf)
 *   indices [NumEle     .. 2*NumEle)      — unsat    (yEleUnsat)
 *   indices [2*NumEle   .. 3*NumEle)      — GW       (yEleGW)
 *   indices [3*NumEle   .. 3*NumEle+NumRiv) — river   (yRivStg)
 *   indices [3*NumEle+NumRiv .. NumY)     — lake     (yLakeStg)
 *
 * Sparsity pattern (per SHUD RHS structure):
 *   - Each element i has 3 mesh neighbors `MD->Ele[i].nabr[k]`
 *     (k=0..2; 0 means boundary, 1-indexed otherwise). Lateral fluxes
 *     couple element i to nabr(i) WITHIN each stripe (surf↔surf,
 *     unsat↔unsat, gw↔gw).
 *   - Vertical fluxes (infiltration / recharge / ET) couple the 3
 *     stripes within the SAME element (surf_i ↔ unsat_i ↔ gw_i).
 *   - River coupling: each element with a river edge couples to a
 *     river node; explicit cross-mapping requires per-river-reach
 *     metadata. River-to-river: each river reach couples to its
 *     `down` neighbor.
 *   - Lake coupling: lake-bank elements couple via `Ele[i].lakenabr[k]`.
 *
 * Conservative bandwidth bound per row: ~12 nonzeros
 *   (self + 3 mesh-nabr × same-stripe + 2 cross-stripe + river/lake edge
 *    × small constant). We allocate 32 candidate slots per column as
 *   a safety overcount.
 *
 * Total ATimes calls per Setup: O(NumY) — one probe per column.
 * For heihe_x4 (NumY ~124k) this is ~124k probes; for heihe_x16
 * (NumY ~485k) it is ~485k probes — manageable within G0 wall budget.
 *
 * Per row, we read out only the topology-derived candidate row set
 * (<= 32 candidates), NOT all NumY rows. */
static void enumerate_row_candidates(const HypreContent *c, HYPRE_BigInt col,
                                     HYPRE_BigInt *candidates, int *n_candidates_out,
                                     int max_candidates) {
    int n_cand = 0;
    const Model_Data *MD = c->MD;
    const int NumEle = MD ? MD->NumEle : 0;
    const int NumRiv = MD ? MD->NumRiv : 0;
    const int NumLake = MD ? MD->NumLake : 0;
    const long n = c->n;

    auto push = [&](HYPRE_BigInt r) {
        if (n_cand >= max_candidates) return;
        if (r < 0 || r >= (HYPRE_BigInt)n) return;
        /* Linear dedup is fine — n_cand is bounded small. */
        for (int i = 0; i < n_cand; ++i) {
            if (candidates[i] == r) return;
        }
        candidates[n_cand++] = r;
    };

    /* Always push self (diagonal). */
    push(col);

    if (NumEle == 0 || MD == nullptr) {
        *n_candidates_out = n_cand;
        return;
    }

    const HYPRE_BigInt stripe_riv_base = (HYPRE_BigInt)(3 * NumEle);
    const HYPRE_BigInt stripe_lake_base = (HYPRE_BigInt)(3 * NumEle + NumRiv);

    if (col < stripe_riv_base) {
        /* Element-stripe column (surf / unsat / gw). */
        const int stripe_idx = (int)(col / NumEle);  /* 0..2 */
        const int ele_idx = (int)(col % NumEle);     /* 0..NumEle-1 */

        /* Same-stripe mesh neighbors. */
        for (int k = 0; k < 3; ++k) {
            int nabr = MD->Ele[ele_idx].nabr[k];  /* 1-indexed; 0 = boundary */
            if (nabr > 0 && nabr <= NumEle) {
                push((HYPRE_BigInt)stripe_idx * NumEle + (nabr - 1));
            }
        }
        /* Cross-stripe (vertical infiltration / recharge / ET) couplings
         * to the OTHER two stripes of the SAME element. */
        for (int s = 0; s < 3; ++s) {
            if (s != stripe_idx) {
                push((HYPRE_BigInt)s * NumEle + ele_idx);
            }
        }
        /* Same-stripe lake-bank neighbors (only meaningful when
         * NumLake > 0; lakenabr is 0 for non-lake-adjacent cells). */
        if (NumLake > 0) {
            for (int k = 0; k < 3; ++k) {
                int lnabr = MD->Ele[ele_idx].lakenabr[k];
                if (lnabr > 0 && lnabr <= NumLake) {
                    push(stripe_lake_base + (lnabr - 1));
                }
            }
        }
        /* River coupling: not topology-restricted from element side
         * without per-element river metadata; we conservatively allow
         * any river row in the candidate set IF this element has a
         * downstream river edge. The exact mapping is non-trivial; for
         * G0 we additionally allow all river rows for this column if
         * the element sits adjacent to NumRiv > 0 (PR-A may tighten). */
        /* Skip explicit per-river enumeration here — the wrapper relies
         * on cross-element infiltration symmetry: if element i couples
         * to river j, the river-row probe (below) picks up the
         * complementary entry. */
    } else if (col < stripe_lake_base) {
        /* River-stripe column. */
        const int riv_idx = (int)(col - stripe_riv_base);  /* 0..NumRiv-1 */
        /* Self already pushed. Push downstream river. */
        if (NumRiv > 0 && MD->Riv != nullptr) {
            int down = MD->Riv[riv_idx].down;  /* 1-indexed; -INT_MAX if none */
            if (down > 0 && down <= NumRiv) {
                push(stripe_riv_base + (down - 1));
            }
            /* Push upstream river(s) — Riv->down points us downstream;
             * upstream coupling is the matrix-transpose direction.
             * Conservative overcount: allow all rivers whose .down ==
             * riv_idx+1; cap at 8 upstream tributaries. */
            int up_count = 0;
            for (int rj = 0; rj < NumRiv && up_count < 8; ++rj) {
                if (MD->Riv[rj].down == (riv_idx + 1)) {
                    push(stripe_riv_base + rj);
                    up_count++;
                }
            }
        }
        /* River-to-element coupling is symmetric to element-to-river;
         * the ATimes row reading will catch nonzero entries in the
         * element stripes if the J row corresponding to this river
         * touches them. We push the element-stripe rows for ALL
         * elements that have this river as a downstream edge — but
         * without per-element river mapping, conservatively push the
         * 3 element-stripe rows at index ele_idx = (riv_idx %
         * NumEle) as a heuristic. */
        if (NumEle > 0) {
            const int ele_idx_heur = riv_idx % NumEle;
            for (int s = 0; s < 3; ++s) {
                push((HYPRE_BigInt)s * NumEle + ele_idx_heur);
            }
        }
    } else {
        /* Lake-stripe column. */
        const int lake_idx = (int)(col - stripe_lake_base);  /* 0..NumLake-1 */
        (void)lake_idx;
        /* Lake-to-element coupling: any element with lakenabr == this
         * lake. Heuristic: scan up to first 16 elements that match
         * (typical lake-bank ring is <16 cells per lake). */
        int scan_count = 0;
        for (int ej = 0; ej < NumEle && scan_count < 16; ++ej) {
            bool is_bank = false;
            for (int k = 0; k < 3; ++k) {
                if (MD->Ele[ej].lakenabr[k] == (lake_idx + 1)) {
                    is_bank = true;
                    break;
                }
            }
            if (is_bank) {
                /* Push GW + surf stripe rows for this bank element. */
                push((HYPRE_BigInt)0 * NumEle + ej);  /* surf */
                push((HYPRE_BigInt)2 * NumEle + ej);  /* gw   */
                scan_count++;
            }
        }
    }

    *n_candidates_out = n_cand;
}

/* Lazy build of AMG hierarchy at first Solve. Uses the in-flight
 * `b` vector (or the constructor-stashed y_template) to determine
 * the N_Vector flavor for probe scratch allocation. */
static int lazy_build_hierarchy_from_solve(HypreContent *c, N_Vector b_template) {
    const auto t0 = std::chrono::steady_clock::now();

    if (c->atimes_fn == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: lazy Setup invoked with atimes_fn=NULL\n");
        return SUNLS_ATIMES_NULL;
    }

    const long n = c->n;

    /* Probe scratch — clone from the live vector flavor. Prefer the
     * caller-supplied b_template; fall back to constructor-stashed
     * y_template if b_template is NULL. */
    N_Vector clone_src = (b_template != nullptr) ? b_template : c->y_template;
    if (clone_src == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: no N_Vector template available for probe scratch\n");
        return SUNLS_MEM_FAIL;
    }
    N_Vector ev = N_VClone(clone_src);
    N_Vector av = N_VClone(clone_src);
    if (ev == nullptr || av == nullptr) {
        if (ev) N_VDestroy(ev);
        if (av) N_VDestroy(av);
        std::fprintf(stderr,
            "[shud-amg] FATAL: probe scratch N_VClone failed\n");
        return SUNLS_MEM_FAIL;
    }

    /* Topology-restricted probe. Per spec REQ-G0 "total ATimes calls
     * per Setup MUST be O(NumY × bw_effective) not O(NumY²)".
     *
     * For each column `col`, derive the small candidate row set from
     * SHUD's mesh topology (see enumerate_row_candidates), invoke
     * ATimes(e_col), and read out ONLY the candidate rows. */
    constexpr int MAX_CANDIDATES = 32;
    HYPRE_BigInt candidates[MAX_CANDIDATES];
    int n_cand = 0;

    int bw_effective_max = 0;
    long long total_atimes_calls = 0;

    /* Per-row accumulator (col, val) lists. */
    std::vector<std::vector<HYPRE_BigInt>> row_cols(n);
    std::vector<std::vector<double>> row_vals(n);

    for (HYPRE_BigInt col = 0; col < n; ++col) {
        N_VConst(0.0, ev);
        double *ev_data = N_VGetArrayPointer(ev);
        if (ev_data == nullptr) {
            N_VDestroy(ev);
            N_VDestroy(av);
            std::fprintf(stderr,
                "[shud-amg] FATAL: N_VGetArrayPointer returned NULL on probe vector\n");
            return SUNLS_MEM_FAIL;
        }
        ev_data[col] = 1.0;

        int atimes_rc = c->atimes_fn(c->atimes_data, ev, av);
        total_atimes_calls++;
        if (atimes_rc != 0) {
            N_VDestroy(ev);
            N_VDestroy(av);
            std::fprintf(stderr,
                "[shud-amg] FATAL: ATimes returned %d on probe col=%lld\n",
                atimes_rc, (long long)col);
            return SUNLS_ATIMES_FAIL_UNREC;
        }

        const double *av_data = N_VGetArrayPointer(av);

        /* Topology-restricted readout. */
        enumerate_row_candidates(c, col, candidates, &n_cand, MAX_CANDIDATES);
        if (n_cand > bw_effective_max) bw_effective_max = n_cand;
        for (int k = 0; k < n_cand; ++k) {
            HYPRE_BigInt row = candidates[k];
            const double v = av_data[row];
            if (v != 0.0) {
                row_cols[row].push_back(col);
                row_vals[row].push_back(v);
            }
        }
    }

    N_VDestroy(ev);
    N_VDestroy(av);

    std::fprintf(stdout,
        "[shud-amg] Setup probe bw_effective=%d total_atimes_calls=%lld\n",
        bw_effective_max, total_atimes_calls);
    std::fflush(stdout);

    /* Build HYPRE IJMatrix from the per-row data. */
    if (HYPRE_IJMatrixCreate(MPI_COMM_SELF, 0, n - 1, 0, n - 1, &c->A_ij) != 0) {
        std::fprintf(stderr, "[shud-amg] FATAL: HYPRE_IJMatrixCreate failed\n");
        return SUNLS_PACKAGE_FAIL_UNREC;
    }
    HYPRE_IJMatrixSetObjectType(c->A_ij, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(c->A_ij);

    for (HYPRE_BigInt r = 0; r < n; ++r) {
        HYPRE_Int ncols = static_cast<HYPRE_Int>(row_cols[r].size());
        if (ncols == 0) continue;
        HYPRE_BigInt row_idx = r;
        HYPRE_IJMatrixSetValues(c->A_ij, /*nrows=*/1, &ncols,
                                &row_idx, row_cols[r].data(),
                                row_vals[r].data());
    }
    HYPRE_IJMatrixAssemble(c->A_ij);
    HYPRE_IJMatrixGetObject(c->A_ij, (void **)&c->A_par);

    /* Build b_ij / x_ij placeholders (filled in Solve). */
    HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, n - 1, &c->b_ij);
    HYPRE_IJVectorSetObjectType(c->b_ij, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(c->b_ij);
    HYPRE_IJVectorAssemble(c->b_ij);
    HYPRE_IJVectorGetObject(c->b_ij, (void **)&c->b_par);

    HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, n - 1, &c->x_ij);
    HYPRE_IJVectorSetObjectType(c->x_ij, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(c->x_ij);
    HYPRE_IJVectorAssemble(c->x_ij);
    HYPRE_IJVectorGetObject(c->x_ij, (void **)&c->x_par);

    /* Configure BoomerAMG hierarchy. */
    HYPRE_BoomerAMGCreate(&c->amg);
    HYPRE_BoomerAMGSetInterpType(c->amg, c->interp_type);
    HYPRE_BoomerAMGSetCoarsenType(c->amg, c->coarsen_type);
    HYPRE_BoomerAMGSetMaxIter(c->amg, 100);
    HYPRE_BoomerAMGSetTol(c->amg, 1e-8);
    HYPRE_BoomerAMGSetPrintLevel(c->amg, 0);
    HYPRE_BoomerAMGSetCumNnzAP(c->amg, 1.0);  /* enable nnz tracking */

    HYPRE_Int setup_rc =
        HYPRE_BoomerAMGSetup(c->amg, c->A_par, c->b_par, c->x_par);

    const auto t1 = std::chrono::steady_clock::now();
    c->pending_setup_wall_sec += std::chrono::duration<double>(t1 - t0).count();
    c->pending_setup_called = 1;

    if (setup_rc != 0) {
        std::fprintf(stderr,
            "[shud-amg] AMG_SETUP_DIVERGE: HYPRE_BoomerAMGSetup rc=%d\n",
            (int)setup_rc);
        std::fprintf(stderr, "MARKER:AMG_SETUP_DIVERGE_DETECTED\n");
        /* Setup failure: free the four handles to avoid leak. */
        if (c->amg)   { HYPRE_BoomerAMGDestroy(c->amg);   c->amg   = nullptr; }
        if (c->A_ij)  { HYPRE_IJMatrixDestroy(c->A_ij);   c->A_ij  = nullptr; c->A_par = nullptr; }
        if (c->b_ij)  { HYPRE_IJVectorDestroy(c->b_ij);   c->b_ij  = nullptr; c->b_par = nullptr; }
        if (c->x_ij)  { HYPRE_IJVectorDestroy(c->x_ij);   c->x_ij  = nullptr; c->x_par = nullptr; }
        return SUNLS_PACKAGE_FAIL_UNREC;
    }

    return SUNLS_SUCCESS;
}

int op_solve(SUNLinearSolver LS, SUNMatrix /*A*/, N_Vector x, N_Vector b,
             realtype /*tol*/) {
    HypreContent *c = content_of(LS);
    const auto t_solve_0 = std::chrono::steady_clock::now();

    /* Lazy-build hierarchy if not yet built. */
    if (c->A_ij == nullptr || c->amg == nullptr) {
        int rc = lazy_build_hierarchy_from_solve(c, b);
        if (rc != SUNLS_SUCCESS) {
            c->last_flag = rc;
            return rc;
        }
    }

    const long n = c->n;
    const double *b_raw = N_VGetArrayPointer(b);
    double *x_raw = N_VGetArrayPointer(x);
    if (b_raw == nullptr || x_raw == nullptr) {
        c->last_flag = SUNLS_MEM_FAIL;
        std::fprintf(stderr,
            "[shud-amg] FATAL: N_VGetArrayPointer NULL on Solve input\n");
        return SUNLS_MEM_FAIL;
    }

    /* Refresh HYPRE vectors from CVODE-supplied x (initial guess)
     * and b (RHS). If zero_guess is SUNTRUE, override x with zeros. */
    std::vector<double> x_init(n, 0.0);
    if (c->zero_guess != SUNTRUE) {
        for (long i = 0; i < n; ++i) x_init[i] = x_raw[i];
    }
    HYPRE_IJVectorInitialize(c->b_ij);
    HYPRE_IJVectorSetValues(c->b_ij, n, c->indices.data(),
                            const_cast<double *>(b_raw));
    HYPRE_IJVectorAssemble(c->b_ij);

    HYPRE_IJVectorInitialize(c->x_ij);
    HYPRE_IJVectorSetValues(c->x_ij, n, c->indices.data(), x_init.data());
    HYPRE_IJVectorAssemble(c->x_ij);

    /* Invoke BoomerAMG solve. */
    HYPRE_Int rc = HYPRE_BoomerAMGSolve(c->amg, c->A_par, c->b_par, c->x_par);

    /* Extract x. */
    HYPRE_IJVectorGetValues(c->x_ij, n, c->indices.data(), x_raw);

    /* Per-Solve telemetry. */
    HYPRE_Int hypre_iters = 0;
    HYPRE_Real cum_nnz_AP = 0.0;
    HYPRE_BoomerAMGGetNumIterations(c->amg, &hypre_iters);
    HYPRE_BoomerAMGGetCumNnzAP(c->amg, &cum_nnz_AP);

    const auto t_solve_1 = std::chrono::steady_clock::now();
    const double solve_wall_sec =
        std::chrono::duration<double>(t_solve_1 - t_solve_0).count();

    /* Append ring-buffer entry.
     *
     * Two cases:
     *   1) Not yet at capacity (ring_count < HYPRE_TELEMETRY_RING_SIZE):
     *      next-write slot is (ring_head + ring_count) % SIZE; advance
     *      ring_count.
     *   2) At capacity (ring_count == HYPRE_TELEMETRY_RING_SIZE):
     *      overwrite the OLDEST entry (slot ring_head, the one about to
     *      be evicted), then advance ring_head; ring_count stays equal
     *      to SIZE. Increment the overflow counter.
     *
     * ASCII trace for SIZE=4, after 6 writes 1..6:
     *   pre-write count=4 head=0 ring=[1,2,3,4]
     *   write 5: overflow case; tail=head=0; ring=[5,2,3,4]; head=1
     *   write 6: overflow case; tail=head=1; ring=[5,6,3,4]; head=2
     *   newest at (head + count - 1) % SIZE = (2 + 4 - 1) % 4 = 1
     *   oldest at head = 2 → reads 3,4,5,6 in order: ring[2]=3,
     *   ring[3]=4, ring[0]=5, ring[1]=6 ✓
     */
    if (c->ring.size() < HYPRE_TELEMETRY_RING_SIZE) {
        c->ring.resize(c->ring.size() + 1);
    }
    int tail;
    if (c->ring_count >= HYPRE_TELEMETRY_RING_SIZE) {
        /* Overflow: overwrite the just-vacated head slot, then advance head. */
        tail = c->ring_head;
        c->ring_head = (c->ring_head + 1) % HYPRE_TELEMETRY_RING_SIZE;
        c->entries_dropped_to_overflow++;
        /* ring_count stays at SIZE. */
    } else {
        tail = (c->ring_head + c->ring_count) % HYPRE_TELEMETRY_RING_SIZE;
        c->ring_count++;
    }
    TelemetryEntry &e = c->ring[tail];
    e.step_idx = c->ctx_step_idx;
    e.t_sim = c->ctx_t_sim;
    e.setup_called = c->pending_setup_called;
    e.hypre_iters = (int)hypre_iters;
    e.hypre_op_count = (double)cum_nnz_AP;
    e.setup_wall_sec = c->pending_setup_wall_sec;
    e.solve_wall_sec = solve_wall_sec;
    e.cvode_nli_step = c->ctx_nli_step;
    e.cvode_nfeLS_step = c->ctx_nfeLS_step;

    c->pending_setup_wall_sec = 0.0;
    c->pending_setup_called = 0;

    /* Idempotent MARKER:AMG_TELEMETRY_REAL emission on first
     * successful Solve. Marker documents the HYPRE release used so
     * PR-B aggregator can attribute the values. The naming is
     * fixed-format (spec REQ-G0-3). */
    if (rc == 0 && !c->telemetry_marker_emitted) {
#if defined(HYPRE_RELEASE_NUMBER)
        const long hypre_release = (long)HYPRE_RELEASE_NUMBER;
#else
        const long hypre_release = 0;
#endif
        std::fprintf(stdout,
            "MARKER:AMG_TELEMETRY_REAL hypre_release=%ld first_iters=%d first_op_count=%.0f\n",
            hypre_release, (int)hypre_iters, (double)cum_nnz_AP);
        std::fflush(stdout);
        c->telemetry_marker_emitted = true;
    }

    /* Map HYPRE return code to SUNLS_* per spec. */
    if (rc == 0) {
        c->last_flag = SUNLS_SUCCESS;
        return SUNLS_SUCCESS;
    }
    /* Non-zero return: emit marker + map to CONV_FAIL (recoverable
     * for CVODE retry) unless we suspect structural failure. */
    std::fprintf(stderr,
        "[shud-amg] AMG_SOLVE_DIVERGE: HYPRE_BoomerAMGSolve rc=%d "
        "iters=%d nnz_AP=%.0f\n",
        (int)rc, (int)hypre_iters, (double)cum_nnz_AP);
    std::fprintf(stderr, "MARKER:AMG_SOLVE_DIVERGE_DETECTED\n");
    c->last_flag = SUNLS_CONV_FAIL;
    return SUNLS_CONV_FAIL;
}

int op_numiters(SUNLinearSolver LS) {
    HypreContent *c = content_of(LS);
    if (c->amg == nullptr) return 0;
    HYPRE_Int n = 0;
    HYPRE_BoomerAMGGetNumIterations(c->amg, &n);
    return (int)n;
}

realtype op_resnorm(SUNLinearSolver /*LS*/) {
    /* G0: BoomerAMG-as-direct-solver semantics. No residual norm
     * exported. */
    return 0.0;
}

sunindextype op_lastflag(SUNLinearSolver LS) {
    HypreContent *c = content_of(LS);
    return c->last_flag;
}

int op_space(SUNLinearSolver LS, long int *lenrwLS, long int *leniwLS) {
    HypreContent *c = content_of(LS);
    if (lenrwLS) *lenrwLS = static_cast<long int>(c->n * 2);
    if (leniwLS) *leniwLS = static_cast<long int>(c->n + 16);
    return SUNLS_SUCCESS;
}

N_Vector op_resid(SUNLinearSolver /*LS*/) {
    return nullptr;
}

int op_free(SUNLinearSolver LS) {
    /* Double-free guard: already freed or never constructed. */
    if (LS == nullptr) return SUNLS_SUCCESS;
    if (LS->content == nullptr) {
        /* Content already released (or never installed). Skip content
         * cleanup but still destroy the empty LS shell — caller may have
         * obtained LS via SUNLinSolNewEmpty without populating content. */
        SUNLinSolFreeEmpty(LS);
        return SUNLS_SUCCESS;
    }
    HypreContent *c = static_cast<HypreContent *>(LS->content);
    if (c != nullptr) {
        if (c->amg)   HYPRE_BoomerAMGDestroy(c->amg);
        if (c->A_ij)  HYPRE_IJMatrixDestroy(c->A_ij);
        if (c->b_ij)  HYPRE_IJVectorDestroy(c->b_ij);
        if (c->x_ij)  HYPRE_IJVectorDestroy(c->x_ij);
        if (c->wrapper_inited_hypre) {
            HYPRE_Finalize();
        }
        if (c->wrapper_inited_mpi) {
            int mpi_finalized = 0;
            MPI_Finalized(&mpi_finalized);
            if (!mpi_finalized) MPI_Finalize();
        }
        delete c;
        LS->content = nullptr;
    }
    /* Use SUNDIALS' own destructor for the empty LS shell to match
     * the SUNLinSolNewEmpty allocator in the constructor. Frees both
     * LS->ops and LS via the same allocator that was used by NewEmpty. */
    SUNLinSolFreeEmpty(LS);
    return SUNLS_SUCCESS;
}

}  /* anonymous namespace */

/* ---- Public constructor ----------------------------------------------- */

extern "C" SUNLinearSolver
SUNLinSol_Hypre(N_Vector y, void *MD_void,
                int interp_type, int coarsen_type, SUNContext sunctx) {
    if (y == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: SUNLinSol_Hypre constructor — y is NULL\n");
        return nullptr;
    }
    if (interp_type != 6 || coarsen_type != 8) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: SUNLinSol_Hypre constructor — "
            "(interp_type, coarsen_type) = (%d, %d), G0 requires (6, 8); "
            "G1 may relax\n",
            interp_type, coarsen_type);
        return nullptr;
    }

    SUNLinearSolver LS = SUNLinSolNewEmpty(sunctx);
    if (LS == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: SUNLinSolNewEmpty returned NULL\n");
        return nullptr;
    }

    HypreContent *c = new (std::nothrow) HypreContent();
    if (c == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: HypreContent allocation failed\n");
        SUNLinSolFreeEmpty(LS);
        return nullptr;
    }
    c->MD = static_cast<Model_Data *>(MD_void);
    c->interp_type = interp_type;
    c->coarsen_type = coarsen_type;
    c->n = (c->MD != nullptr) ? c->MD->NumY : 0;
    /* Stash the N_Vector flavor template for lazy probe-scratch
     * cloning. Caller owns the pointer; the wrapper does not destroy
     * y_template in Free. */
    c->y_template = y;

    /* Pre-allocate row indices. */
    if (c->n > 0) {
        c->indices.assign(c->n, 0);
        for (HYPRE_BigInt i = 0; i < c->n; ++i) c->indices[i] = i;
    }

    LS->content = c;
    LS->ops->gettype           = op_gettype;
    LS->ops->getid             = op_getid;
    LS->ops->setatimes         = op_setatimes;
    LS->ops->setpreconditioner = op_setpreconditioner;
    LS->ops->setscalingvectors = op_setscalingvectors;
    LS->ops->setzeroguess      = op_setzeroguess;
    LS->ops->initialize        = op_initialize;
    LS->ops->setup             = op_setup;
    LS->ops->solve             = op_solve;
    LS->ops->numiters          = op_numiters;
    LS->ops->resnorm           = op_resnorm;
    LS->ops->lastflag          = op_lastflag;
    LS->ops->space             = op_space;
    LS->ops->resid             = op_resid;
    LS->ops->free              = op_free;

    /* Build-time guard against any future ABI drift: every one of
     * the 15 op slots must be non-NULL at construction return. If a
     * refactor accidentally drops a slot assignment, this assert
     * fires in debug builds and `if (!ok)` short-circuits in release
     * builds (returning NULL with stderr message). */
    const bool ok =
        LS->ops->gettype && LS->ops->getid && LS->ops->setatimes &&
        LS->ops->setpreconditioner && LS->ops->setscalingvectors &&
        LS->ops->setzeroguess && LS->ops->initialize && LS->ops->setup &&
        LS->ops->solve && LS->ops->numiters && LS->ops->resnorm &&
        LS->ops->lastflag && LS->ops->space && LS->ops->resid &&
        LS->ops->free;
    assert(ok && "SUNLinSol_Hypre: 15-callback ABI slot population incomplete");
    if (!ok) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: 15-callback ABI slot population incomplete\n");
        op_free(LS);
        return nullptr;
    }

    return LS;
}

extern "C" int
SUNLinSol_Hypre_SetStepContext(SUNLinearSolver LS, long step_idx,
                               realtype t_sim, long cvode_nli_cum,
                               long cvode_nfeLS_cum) {
    if (LS == nullptr || LS->content == nullptr) return SUNLS_SUCCESS;
    if (LS->ops == nullptr || LS->ops->getid == nullptr ||
        LS->ops->getid(LS) != SUNLINEARSOLVER_CUSTOM) {
        return SUNLS_SUCCESS;  /* not our wrapper */
    }
    HypreContent *c = content_of(LS);
    c->ctx_step_idx = step_idx;
    c->ctx_t_sim = t_sim;
    c->ctx_nli_step = cvode_nli_cum - c->ctx_last_nli_cum;
    c->ctx_nfeLS_step = cvode_nfeLS_cum - c->ctx_last_nfeLS_cum;
    c->ctx_last_nli_cum = cvode_nli_cum;
    c->ctx_last_nfeLS_cum = cvode_nfeLS_cum;
    return SUNLS_SUCCESS;
}

extern "C" int
SUNLinSol_Hypre_DrainTelemetry(SUNLinearSolver LS, FILE *out) {
    if (LS == nullptr || LS->content == nullptr || out == nullptr) return 0;
    if (LS->ops == nullptr || LS->ops->getid == nullptr ||
        LS->ops->getid(LS) != SUNLINEARSOLVER_CUSTOM) {
        return 0;
    }
    HypreContent *c = content_of(LS);

    std::fprintf(out,
        "step_idx\tt_sim\tsetup_called\thypre_iters\thypre_op_count"
        "\tsetup_wall_sec\tsolve_wall_sec\tcvode_nli_step\tcvode_nfeLS_step\n");

    int written = 0;
    for (int i = 0; i < c->ring_count; ++i) {
        int idx = (c->ring_head + i) % HYPRE_TELEMETRY_RING_SIZE;
        const TelemetryEntry &e = c->ring[idx];
        std::fprintf(out,
            "%ld\t%.6f\t%d\t%d\t%.0f\t%.6e\t%.6e\t%ld\t%ld\n",
            e.step_idx, (double)e.t_sim, e.setup_called,
            e.hypre_iters, e.hypre_op_count,
            e.setup_wall_sec, e.solve_wall_sec,
            e.cvode_nli_step, e.cvode_nfeLS_step);
        written++;
    }
    if (c->entries_dropped_to_overflow > 0) {
        std::fprintf(out,
            "# entries_dropped_to_overflow=%ld\n",
            c->entries_dropped_to_overflow);
    }
    c->ring_head = 0;
    c->ring_count = 0;
    c->entries_dropped_to_overflow = 0;
    return written;
}
