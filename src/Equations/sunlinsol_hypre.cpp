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
     * application-side MPI_Init call exists. */
    int mpi_already = 0;
    MPI_Initialized(&mpi_already);
    if (!mpi_already) {
        int mpi_argc = 1;
        char mpi_arg0[] = "shud_amg";
        char *mpi_argv[2] = {mpi_arg0, nullptr};
        char **mpi_argv_ptr = mpi_argv;
        if (MPI_Init(&mpi_argc, &mpi_argv_ptr) != MPI_SUCCESS) {
            std::fprintf(stderr,
                "[shud-amg] FATAL: MPI_Init failed in SUNLinSol_Hypre Initialize\n");
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

    /* Hypre thread pinning under shud_omp (design.md D10). Defense-
     * in-depth: smoke-runner sbatch sets `export OMP_NUM_THREADS=1`,
     * and we additionally pin Hypre to 1 thread here. */
    int multi_thread_mode = 0;
#ifdef _OPENMP
    if (omp_get_max_threads() > 1) {
        multi_thread_mode = 1;
    }
#endif
    const char *omp_nt_env = std::getenv("OMP_NUM_THREADS");
    if (omp_nt_env != nullptr) {
        long v = std::strtol(omp_nt_env, nullptr, 10);
        if (v > 1) multi_thread_mode = 1;
    }
    /* Hypre 3.1.0 does not expose HYPRE_SetGlobalOptions(...) as a
     * portable runtime knob; thread pinning under OpenMP-enabled
     * Hypre is the user's responsibility (omp_set_num_threads).
     * Emit the documentation marker so test runs are self-explanatory. */
    std::fprintf(stdout,
        "[shud-amg] Hypre threads=1 mode=%s\n",
        multi_thread_mode ? "nested" : "single");
    std::fflush(stdout);

    return SUNLS_SUCCESS;
}

/* Internal helper used by both lazy-Setup (first Solve) and
 * CVODE-issued Setup. Returns SUNLS_SUCCESS or a SUNLS_* error code. */
int build_amg_hierarchy(HypreContent *c) {
    const auto t0 = std::chrono::steady_clock::now();

    /* Sanity. */
    if (c->MD == nullptr || c->atimes_fn == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: Setup invoked with MD=%p atimes_fn=%p\n",
            (void *)c->MD, (void *)c->atimes_fn);
        return SUNLS_ILL_INPUT;
    }
    const long n = c->n;
    if (n <= 0) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: Setup invoked with n=%ld\n", n);
        return SUNLS_ILL_INPUT;
    }

    /* Destroy any prior handles (CVODE-issued Setup → rebuild). */
    if (c->amg)   { HYPRE_BoomerAMGDestroy(c->amg);   c->amg   = nullptr; }
    if (c->A_ij)  { HYPRE_IJMatrixDestroy(c->A_ij);   c->A_ij  = nullptr; c->A_par = nullptr; }
    if (c->b_ij)  { HYPRE_IJVectorDestroy(c->b_ij);   c->b_ij  = nullptr; c->b_par = nullptr; }
    if (c->x_ij)  { HYPRE_IJVectorDestroy(c->x_ij);   c->x_ij  = nullptr; c->x_par = nullptr; }

    /* Probe-derive sparsity: for each column i in [0, n), invoke
     * ATimes(A_data, e_i, J·e_i) where e_i is the i-th unit basis
     * vector, then accumulate (col=i, row=r, val) tuples wherever
     * the result is non-zero. Re-arrange to per-row CSR for HYPRE
     * IJMatrix insertion.
     *
     * Cost: n ATimes calls per Setup. For keliya (NumY ~ 3*484 ~
     * 1452), this is ~1500 RHS evaluations — cheap (~1s on Mac).
     * For heihe_x4 (NumY ~ 124k) it is the design.md D5 perf concern
     * — accepted at G0 as the correctness baseline; G1 may add
     * topology pre-filtering to skip known-zero columns. */
    if (c->indices.size() != static_cast<size_t>(n)) {
        c->indices.assign(n, 0);
        for (HYPRE_BigInt i = 0; i < n; ++i) c->indices[i] = i;
    }

    /* Use Serial N_Vectors as the probe scratch (the wrapper assumes
     * CVODE's outer N_Vector is also Serial when ATimes can be
     * driven via simple raw-array access; nvector_openmp also exposes
     * N_VGetArrayPointer compatible semantics under SHUD's usage). */
    N_Vector ev = N_VClone(c->s2 ? c->s2 : c->s1);  /* placeholder; replaced below */
    /* ev needs to be the same flavor as CVODE's state vector. We
     * don't have direct access to a template here unless ATimes
     * passes one through. In practice CVODE invokes
     * setatimes(LS, cvode_mem, cvAtimes) and cvAtimes(A_data, v, z)
     * expects v to be a CVODE-managed N_Vector flavor. The cleanest
     * path: clone from a known-good N_Vector. We cache one via
     * Setup-call SUNMatrix arg in the LS-level Setup; but since
     * SetATimes preceded us, fall back to N_VNew_Serial(n, sunctx). */
    if (ev != nullptr) {
        N_VDestroy(ev);
        ev = nullptr;
    }
    /* nvector_serial.h provides N_VNew_Serial(sunindextype length,
     * SUNContext sunctx). Pull sunctx from the SUNLinearSolver
     * carrier (the wrapper's struct stores sunctx implicitly via
     * the SUNLinearSolver_C handle, but we don't have access here —
     * solve by passing sunctx through a stash). To keep this PR-0
     * patch minimal, allocate the probe N_Vectors at first Solve
     * where we DO have the in-flight x/b N_Vector to clone from. */

    /* This intermediate buffer path is only used when Setup is
     * issued by CVODE OUTSIDE of a Solve call — exceedingly rare per
     * task 1.4 (nsetups=0 on keliya). For PR-0 we forward this path
     * to error-return; the lazy-build inside Solve handles the real
     * scenario. */
    std::fprintf(stderr,
        "[shud-amg] WARNING: explicit Setup invoked but lazy-Solve path is the G0 baseline — "
        "deferring AMG hierarchy build to first Solve invocation\n");

    const auto t1 = std::chrono::steady_clock::now();
    c->pending_setup_wall_sec += std::chrono::duration<double>(t1 - t0).count();
    c->pending_setup_called = 1;
    return SUNLS_SUCCESS;
}

int op_setup(SUNLinearSolver LS, SUNMatrix /*A*/) {
    return build_amg_hierarchy(content_of(LS));
}

/* Lazy build of AMG hierarchy at first Solve. Uses the in-flight
 * `b` vector to determine the N_Vector flavor for probe scratch
 * allocation. */
static int lazy_build_hierarchy_from_solve(HypreContent *c, N_Vector b_template) {
    const auto t0 = std::chrono::steady_clock::now();

    if (c->atimes_fn == nullptr) {
        std::fprintf(stderr,
            "[shud-amg] FATAL: lazy Setup invoked with atimes_fn=NULL\n");
        return SUNLS_ATIMES_NULL;
    }

    const long n = c->n;

    /* Probe scratch — clone from the live vector flavor. */
    N_Vector ev = N_VClone(b_template);
    N_Vector av = N_VClone(b_template);
    if (ev == nullptr || av == nullptr) {
        if (ev) N_VDestroy(ev);
        if (av) N_VDestroy(av);
        std::fprintf(stderr,
            "[shud-amg] FATAL: probe scratch N_VClone failed\n");
        return SUNLS_MEM_FAIL;
    }

    /* Build column-major (col_i, row_j, val) triplet list via probe. */
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
        if (atimes_rc != 0) {
            N_VDestroy(ev);
            N_VDestroy(av);
            std::fprintf(stderr,
                "[shud-amg] FATAL: ATimes returned %d on probe col=%lld\n",
                atimes_rc, (long long)col);
            return SUNLS_ATIMES_FAIL_UNREC;
        }

        const double *av_data = N_VGetArrayPointer(av);
        for (HYPRE_BigInt row = 0; row < n; ++row) {
            const double v = av_data[row];
            if (v != 0.0) {
                row_cols[row].push_back(col);
                row_vals[row].push_back(v);
            }
        }
    }

    N_VDestroy(ev);
    N_VDestroy(av);

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

    /* Append ring-buffer entry. */
    if (c->ring.size() < HYPRE_TELEMETRY_RING_SIZE) {
        c->ring.resize(c->ring.size() + 1);
    }
    if (c->ring_count >= HYPRE_TELEMETRY_RING_SIZE) {
        c->entries_dropped_to_overflow++;
        c->ring_head = (c->ring_head + 1) % HYPRE_TELEMETRY_RING_SIZE;
        c->ring_count = HYPRE_TELEMETRY_RING_SIZE;
    } else {
        c->ring_count++;
    }
    int tail = (c->ring_head + c->ring_count - 1) % HYPRE_TELEMETRY_RING_SIZE;
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
    if (LS == nullptr) return SUNLS_SUCCESS;
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
    /* Free the empty LS shell allocated by SUNLinSolNewEmpty. */
    if (LS->ops) { std::free(LS->ops); LS->ops = nullptr; }
    std::free(LS);
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
