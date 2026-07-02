/* =====================================================================
 * P12-nvec PR-N0 (#442) — env-gated NVector op-share profiler (impl).
 * See MD_nvec_prof.hpp for the design rationale, composition-order rule,
 * clone-propagation contract, and bitwise-neutrality guarantee.
 * ===================================================================== */

#include "MD_nvec_prof.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

/* ---------------------------------------------------------------------
 * Monotonic nanosecond clock. CLOCK_MONOTONIC where available (Linux
 * server + modern macOS); a steady fallback keeps the profiler building
 * everywhere. Only ever used inside a shim, i.e. only when the profiler
 * is ON, so it never touches the default (gate-off) hot path.
 * ------------------------------------------------------------------- */
static inline uint64_t nvp_now_ns(void) {
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
#else
    return (uint64_t)clock() * (1000000000ull / CLOCKS_PER_SEC);
#endif
}

/* ---------------------------------------------------------------------
 * Fixed op registry. Each WRAPPED op slot has: a stable name (== the CSV
 * op_name and the ops-struct field name), a fixed op_class, a captured
 * ORIGINAL function pointer, and per-op counters. The op_class column is
 * the source-committed mapping table the spec requires (Requirement
 * "nvec_prof.csv pinned schema").
 *
 * CLASSIFICATION (the fixed mapping table — audit deliverable):
 *   elementwise : output element z[i] is a fixed per-element expression
 *                 of inputs at the SAME index i (no cross-element
 *                 accumulation). Parallelizable bitwise across threads.
 *   reduction   : accumulates a scalar (or per-element test) ACROSS all
 *                 elements — the class whose OpenMP `reduction(...)` order
 *                 varies with thread count (the Config D failure mode).
 *   other       : structural / clone ops that carry op work but are
 *                 neither of the above (nvclone, nvcloneempty).
 * ------------------------------------------------------------------- */
enum NvpOp {
    /* --- elementwise (standard) --- */
    NVP_LINEARSUM = 0,
    NVP_CONST,
    NVP_PROD,
    NVP_DIV,
    NVP_SCALE,
    NVP_ABS,
    NVP_INV,
    NVP_ADDCONST,
    NVP_COMPARE,
    /* --- reduction (standard) --- */
    NVP_DOTPROD,
    NVP_MAXNORM,
    NVP_WRMSNORM,
    NVP_WRMSNORMMASK,
    NVP_MIN,
    NVP_WL2NORM,
    NVP_L1NORM,
    NVP_INVTEST,
    NVP_CONSTRMASK,
    NVP_MINQUOTIENT,
    /* --- reduction (local kernels, may be populated by OpenMP backend) --- */
    NVP_DOTPRODLOCAL,
    NVP_MAXNORMLOCAL,
    NVP_MINLOCAL,
    NVP_L1NORMLOCAL,
    NVP_INVTESTLOCAL,
    NVP_CONSTRMASKLOCAL,
    NVP_MINQUOTIENTLOCAL,
    NVP_WSQRSUMLOCAL,
    NVP_WSQRSUMMASKLOCAL,
    /* --- fused / vector-array (NULL by default; wrapped only if enabled) --- */
    NVP_LINEARCOMBINATION,
    NVP_SCALEADDMULTI,
    NVP_DOTPRODMULTI,
    NVP_LINEARSUMVECTORARRAY,
    NVP_SCALEVECTORARRAY,
    NVP_CONSTVECTORARRAY,
    NVP_WRMSNORMVECTORARRAY,
    NVP_WRMSNORMMASKVECTORARRAY,
    NVP_DOTPRODMULTILOCAL,
    /* --- structural --- */
    NVP_CLONE,
    NVP_CLONEEMPTY,
    NVP_COUNT
};

struct NvpSlot {
    const char *name;
    const char *op_class;
    void       *orig;        /* captured original function pointer   */
    uint64_t    calls;
    uint64_t    total_ns;
    bool        wrapped;     /* this slot was populated & wrapped     */
};

/* The names + classes are fixed at compile time; orig/calls/total_ns/
 * wrapped are filled at install(). Indexed by enum NvpOp. */
static NvpSlot g_slot[NVP_COUNT] = {
    { "nvlinearsum",   "elementwise", NULL, 0, 0, false },
    { "nvconst",       "elementwise", NULL, 0, 0, false },
    { "nvprod",        "elementwise", NULL, 0, 0, false },
    { "nvdiv",         "elementwise", NULL, 0, 0, false },
    { "nvscale",       "elementwise", NULL, 0, 0, false },
    { "nvabs",         "elementwise", NULL, 0, 0, false },
    { "nvinv",         "elementwise", NULL, 0, 0, false },
    { "nvaddconst",    "elementwise", NULL, 0, 0, false },
    { "nvcompare",     "elementwise", NULL, 0, 0, false },
    { "nvdotprod",     "reduction",   NULL, 0, 0, false },
    { "nvmaxnorm",     "reduction",   NULL, 0, 0, false },
    { "nvwrmsnorm",    "reduction",   NULL, 0, 0, false },
    { "nvwrmsnormmask","reduction",   NULL, 0, 0, false },
    { "nvmin",         "reduction",   NULL, 0, 0, false },
    { "nvwl2norm",     "reduction",   NULL, 0, 0, false },
    { "nvl1norm",      "reduction",   NULL, 0, 0, false },
    { "nvinvtest",     "reduction",   NULL, 0, 0, false },
    { "nvconstrmask",  "reduction",   NULL, 0, 0, false },
    { "nvminquotient", "reduction",   NULL, 0, 0, false },
    { "nvdotprodlocal",     "reduction", NULL, 0, 0, false },
    { "nvmaxnormlocal",     "reduction", NULL, 0, 0, false },
    { "nvminlocal",         "reduction", NULL, 0, 0, false },
    { "nvl1normlocal",      "reduction", NULL, 0, 0, false },
    { "nvinvtestlocal",     "reduction", NULL, 0, 0, false },
    { "nvconstrmasklocal",  "reduction", NULL, 0, 0, false },
    { "nvminquotientlocal", "reduction", NULL, 0, 0, false },
    { "nvwsqrsumlocal",     "reduction", NULL, 0, 0, false },
    { "nvwsqrsummasklocal", "reduction", NULL, 0, 0, false },
    { "nvlinearcombination",         "elementwise", NULL, 0, 0, false },
    { "nvscaleaddmulti",             "elementwise", NULL, 0, 0, false },
    { "nvdotprodmulti",              "reduction",   NULL, 0, 0, false },
    { "nvlinearsumvectorarray",      "elementwise", NULL, 0, 0, false },
    { "nvscalevectorarray",          "elementwise", NULL, 0, 0, false },
    { "nvconstvectorarray",          "elementwise", NULL, 0, 0, false },
    { "nvwrmsnormvectorarray",       "reduction",   NULL, 0, 0, false },
    { "nvwrmsnormmaskvectorarray",   "reduction",   NULL, 0, 0, false },
    { "nvdotprodmultilocal",         "reduction",   NULL, 0, 0, false },
    { "nvclone",       "other", NULL, 0, 0, false },
    { "nvcloneempty",  "other", NULL, 0, 0, false },
};

/* Gate state, read once at the first install() and cached. */
static bool g_prof_on       = false;
static bool g_prof_resolved = false;

bool nvec_prof_env_is_1(const char *name) {
    const char *v = getenv(name);
    return (v != NULL) && (strcmp(v, "1") == 0);
}

bool nvec_prof_is_on(void) {
    if (!g_prof_resolved) {
        g_prof_on = nvec_prof_env_is_1("SHUD_NVEC_PROF");
        g_prof_resolved = true;
    }
    return g_prof_on;
}

/* Record one shim invocation. Called with the elapsed ns already
 * measured around the delegate call. */
static inline void nvp_account(int op, uint64_t dt_ns) {
    g_slot[op].calls    += 1;
    g_slot[op].total_ns += dt_ns;
}

/* ---------------------------------------------------------------------
 * Shims. C function pointers cannot capture state, so every wrapped slot
 * gets its own named shim that reads its captured original from the fixed
 * registry, times the delegate call, and returns the delegate's result.
 * Grouped by ops-struct signature; one macro per signature keeps the
 * boilerplate honest (identical structure → identical timing semantics).
 *
 * Casting through the exact field type is required: the ops table stores
 * distinct function-pointer types per slot.
 * ------------------------------------------------------------------- */

/* void (realtype, N_Vector, realtype, N_Vector, N_Vector) — linearsum */
static void nvp_linearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z) {
    typedef void (*F)(realtype, N_Vector, realtype, N_Vector, N_Vector);
    uint64_t t0 = nvp_now_ns();
    ((F)g_slot[NVP_LINEARSUM].orig)(a, x, b, y, z);
    nvp_account(NVP_LINEARSUM, nvp_now_ns() - t0);
}

/* void (realtype, N_Vector) — const */
static void nvp_const(realtype c, N_Vector z) {
    typedef void (*F)(realtype, N_Vector);
    uint64_t t0 = nvp_now_ns();
    ((F)g_slot[NVP_CONST].orig)(c, z);
    nvp_account(NVP_CONST, nvp_now_ns() - t0);
}

/* void (N_Vector, N_Vector, N_Vector) — prod / div */
#define NVP_SHIM_VVV(fn, IDX)                                          \
    static void fn(N_Vector x, N_Vector y, N_Vector z) {              \
        typedef void (*F)(N_Vector, N_Vector, N_Vector);             \
        uint64_t t0 = nvp_now_ns();                                   \
        ((F)g_slot[IDX].orig)(x, y, z);                               \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
    }
NVP_SHIM_VVV(nvp_prod, NVP_PROD)
NVP_SHIM_VVV(nvp_div,  NVP_DIV)

/* void (realtype, N_Vector, N_Vector) — scale / compare */
#define NVP_SHIM_RVV(fn, IDX)                                          \
    static void fn(realtype c, N_Vector x, N_Vector z) {             \
        typedef void (*F)(realtype, N_Vector, N_Vector);            \
        uint64_t t0 = nvp_now_ns();                                   \
        ((F)g_slot[IDX].orig)(c, x, z);                               \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
    }
NVP_SHIM_RVV(nvp_scale,   NVP_SCALE)
NVP_SHIM_RVV(nvp_compare, NVP_COMPARE)

/* void (N_Vector, N_Vector) — abs / inv */
#define NVP_SHIM_VV(fn, IDX)                                           \
    static void fn(N_Vector x, N_Vector z) {                         \
        typedef void (*F)(N_Vector, N_Vector);                      \
        uint64_t t0 = nvp_now_ns();                                   \
        ((F)g_slot[IDX].orig)(x, z);                                  \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
    }
NVP_SHIM_VV(nvp_abs, NVP_ABS)
NVP_SHIM_VV(nvp_inv, NVP_INV)

/* void (N_Vector, realtype, N_Vector) — addconst */
static void nvp_addconst(N_Vector x, realtype b, N_Vector z) {
    typedef void (*F)(N_Vector, realtype, N_Vector);
    uint64_t t0 = nvp_now_ns();
    ((F)g_slot[NVP_ADDCONST].orig)(x, b, z);
    nvp_account(NVP_ADDCONST, nvp_now_ns() - t0);
}

/* realtype (N_Vector, N_Vector) — dotprod / wrmsnorm / wl2norm /
 * minquotient + local dotprod / minquotient / wsqrsum */
#define NVP_SHIM_R_VV(fn, IDX)                                         \
    static realtype fn(N_Vector x, N_Vector y) {                     \
        typedef realtype (*F)(N_Vector, N_Vector);                  \
        uint64_t t0 = nvp_now_ns();                                   \
        realtype r = ((F)g_slot[IDX].orig)(x, y);                     \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_R_VV(nvp_dotprod,          NVP_DOTPROD)
NVP_SHIM_R_VV(nvp_wrmsnorm,         NVP_WRMSNORM)
NVP_SHIM_R_VV(nvp_wl2norm,          NVP_WL2NORM)
NVP_SHIM_R_VV(nvp_minquotient,      NVP_MINQUOTIENT)
NVP_SHIM_R_VV(nvp_dotprodlocal,     NVP_DOTPRODLOCAL)
NVP_SHIM_R_VV(nvp_minquotientlocal, NVP_MINQUOTIENTLOCAL)
NVP_SHIM_R_VV(nvp_wsqrsumlocal,     NVP_WSQRSUMLOCAL)

/* realtype (N_Vector) — maxnorm / min / l1norm + local variants */
#define NVP_SHIM_R_V(fn, IDX)                                          \
    static realtype fn(N_Vector x) {                                 \
        typedef realtype (*F)(N_Vector);                            \
        uint64_t t0 = nvp_now_ns();                                   \
        realtype r = ((F)g_slot[IDX].orig)(x);                        \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_R_V(nvp_maxnorm,      NVP_MAXNORM)
NVP_SHIM_R_V(nvp_min,          NVP_MIN)
NVP_SHIM_R_V(nvp_l1norm,       NVP_L1NORM)
NVP_SHIM_R_V(nvp_maxnormlocal, NVP_MAXNORMLOCAL)
NVP_SHIM_R_V(nvp_minlocal,     NVP_MINLOCAL)
NVP_SHIM_R_V(nvp_l1normlocal,  NVP_L1NORMLOCAL)

/* realtype (N_Vector, N_Vector, N_Vector) — wrmsnormmask + wsqrsummasklocal */
#define NVP_SHIM_R_VVV(fn, IDX)                                        \
    static realtype fn(N_Vector x, N_Vector w, N_Vector id) {        \
        typedef realtype (*F)(N_Vector, N_Vector, N_Vector);        \
        uint64_t t0 = nvp_now_ns();                                   \
        realtype r = ((F)g_slot[IDX].orig)(x, w, id);                 \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_R_VVV(nvp_wrmsnormmask,     NVP_WRMSNORMMASK)
NVP_SHIM_R_VVV(nvp_wsqrsummasklocal, NVP_WSQRSUMMASKLOCAL)

/* booleantype (N_Vector, N_Vector) — invtest + invtestlocal */
#define NVP_SHIM_B_VV(fn, IDX)                                         \
    static booleantype fn(N_Vector x, N_Vector z) {                  \
        typedef booleantype (*F)(N_Vector, N_Vector);               \
        uint64_t t0 = nvp_now_ns();                                   \
        booleantype r = ((F)g_slot[IDX].orig)(x, z);                  \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_B_VV(nvp_invtest,      NVP_INVTEST)
NVP_SHIM_B_VV(nvp_invtestlocal, NVP_INVTESTLOCAL)

/* booleantype (N_Vector, N_Vector, N_Vector) — constrmask + constrmasklocal */
#define NVP_SHIM_B_VVV(fn, IDX)                                        \
    static booleantype fn(N_Vector c, N_Vector x, N_Vector m) {      \
        typedef booleantype (*F)(N_Vector, N_Vector, N_Vector);     \
        uint64_t t0 = nvp_now_ns();                                   \
        booleantype r = ((F)g_slot[IDX].orig)(c, x, m);               \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_B_VVV(nvp_constrmask,      NVP_CONSTRMASK)
NVP_SHIM_B_VVV(nvp_constrmasklocal, NVP_CONSTRMASKLOCAL)

/* N_Vector (N_Vector) — clone / cloneempty. Return value is a fresh
 * vector whose ops table was copied from the (already-wrapped) source, so
 * the clone inherits the shims automatically — nothing extra to do here
 * beyond timing + delegation. */
#define NVP_SHIM_NV_V(fn, IDX)                                         \
    static N_Vector fn(N_Vector w) {                                 \
        typedef N_Vector (*F)(N_Vector);                            \
        uint64_t t0 = nvp_now_ns();                                   \
        N_Vector r = ((F)g_slot[IDX].orig)(w);                        \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_NV_V(nvp_clone,      NVP_CLONE)
NVP_SHIM_NV_V(nvp_cloneempty, NVP_CLONEEMPTY)

/* Fused / vector-array shims. These are NULL by default in the OpenMP and
 * Serial backends (CVODE only populates them via N_VEnable*Ops, which
 * SHUD does not call), so in practice they are never installed. They are
 * wrapped defensively IF a future config enables them, keeping the "wrap
 * every populated op" contract complete. */
#define NVP_SHIM_FUSED_LC(fn, IDX)                                    \
    static int fn(int n, realtype *c, N_Vector *X, N_Vector z) {     \
        typedef int (*F)(int, realtype *, N_Vector *, N_Vector);    \
        uint64_t t0 = nvp_now_ns();                                   \
        int r = ((F)g_slot[IDX].orig)(n, c, X, z);                    \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_FUSED_LC(nvp_linearcombination, NVP_LINEARCOMBINATION)

static int nvp_scaleaddmulti(int n, realtype *a, N_Vector x, N_Vector *Y, N_Vector *Z) {
    typedef int (*F)(int, realtype *, N_Vector, N_Vector *, N_Vector *);
    uint64_t t0 = nvp_now_ns();
    int r = ((F)g_slot[NVP_SCALEADDMULTI].orig)(n, a, x, Y, Z);
    nvp_account(NVP_SCALEADDMULTI, nvp_now_ns() - t0);
    return r;
}

/* int (int, N_Vector, N_Vector*, realtype*) — dotprodmulti + dotprodmultilocal */
#define NVP_SHIM_DPM(fn, IDX)                                          \
    static int fn(int n, N_Vector x, N_Vector *Y, realtype *d) {     \
        typedef int (*F)(int, N_Vector, N_Vector *, realtype *);    \
        uint64_t t0 = nvp_now_ns();                                   \
        int r = ((F)g_slot[IDX].orig)(n, x, Y, d);                    \
        nvp_account(IDX, nvp_now_ns() - t0);                          \
        return r;                                                     \
    }
NVP_SHIM_DPM(nvp_dotprodmulti,      NVP_DOTPRODMULTI)
NVP_SHIM_DPM(nvp_dotprodmultilocal, NVP_DOTPRODMULTILOCAL)

static int nvp_linearsumvectorarray(int n, realtype a, N_Vector *X, realtype b, N_Vector *Y, N_Vector *Z) {
    typedef int (*F)(int, realtype, N_Vector *, realtype, N_Vector *, N_Vector *);
    uint64_t t0 = nvp_now_ns();
    int r = ((F)g_slot[NVP_LINEARSUMVECTORARRAY].orig)(n, a, X, b, Y, Z);
    nvp_account(NVP_LINEARSUMVECTORARRAY, nvp_now_ns() - t0);
    return r;
}

static int nvp_scalevectorarray(int n, realtype *c, N_Vector *X, N_Vector *Z) {
    typedef int (*F)(int, realtype *, N_Vector *, N_Vector *);
    uint64_t t0 = nvp_now_ns();
    int r = ((F)g_slot[NVP_SCALEVECTORARRAY].orig)(n, c, X, Z);
    nvp_account(NVP_SCALEVECTORARRAY, nvp_now_ns() - t0);
    return r;
}

static int nvp_constvectorarray(int n, realtype c, N_Vector *Z) {
    typedef int (*F)(int, realtype, N_Vector *);
    uint64_t t0 = nvp_now_ns();
    int r = ((F)g_slot[NVP_CONSTVECTORARRAY].orig)(n, c, Z);
    nvp_account(NVP_CONSTVECTORARRAY, nvp_now_ns() - t0);
    return r;
}

static int nvp_wrmsnormvectorarray(int n, N_Vector *X, N_Vector *W, realtype *nrm) {
    typedef int (*F)(int, N_Vector *, N_Vector *, realtype *);
    uint64_t t0 = nvp_now_ns();
    int r = ((F)g_slot[NVP_WRMSNORMVECTORARRAY].orig)(n, X, W, nrm);
    nvp_account(NVP_WRMSNORMVECTORARRAY, nvp_now_ns() - t0);
    return r;
}

static int nvp_wrmsnormmaskvectorarray(int n, N_Vector *X, N_Vector *W, N_Vector id, realtype *nrm) {
    typedef int (*F)(int, N_Vector *, N_Vector *, N_Vector, realtype *);
    uint64_t t0 = nvp_now_ns();
    int r = ((F)g_slot[NVP_WRMSNORMMASKVECTORARRAY].orig)(n, X, W, id, nrm);
    nvp_account(NVP_WRMSNORMMASKVECTORARRAY, nvp_now_ns() - t0);
    return r;
}

/* ---------------------------------------------------------------------
 * install(): capture-and-replace one slot iff it is populated. `field`
 * is the live ops-table entry (lvalue); `shim` is the matching typed
 * shim. If the same original was already captured (du after udata, or a
 * clone re-entering), the capture is idempotent — we detect an already-
 * installed shim by pointer identity and skip re-capturing (which would
 * otherwise store the shim as its own "original" and infinite-loop).
 * ------------------------------------------------------------------- */
#define NVP_WRAP(FIELD, IDX, SHIM)                                     \
    do {                                                               \
        if ((FIELD) != NULL && (void *)(FIELD) != (void *)(SHIM)) {   \
            g_slot[IDX].orig    = (void *)(FIELD);                     \
            g_slot[IDX].wrapped = true;                                \
            (FIELD) = (SHIM);                                          \
        } else if ((void *)(FIELD) == (void *)(SHIM)) {               \
            g_slot[IDX].wrapped = true; /* already wrapped (idempotent) */ \
        }                                                              \
    } while (0)

void nvec_prof_install(N_Vector v, const char *backend) {
    (void)backend; /* backend recorded at dump time */
    if (!nvec_prof_is_on()) return;
    if (v == NULL || v->ops == NULL) return;

    N_Vector_Ops o = v->ops;

    /* elementwise */
    NVP_WRAP(o->nvlinearsum, NVP_LINEARSUM, nvp_linearsum);
    NVP_WRAP(o->nvconst,     NVP_CONST,     nvp_const);
    NVP_WRAP(o->nvprod,      NVP_PROD,      nvp_prod);
    NVP_WRAP(o->nvdiv,       NVP_DIV,       nvp_div);
    NVP_WRAP(o->nvscale,     NVP_SCALE,     nvp_scale);
    NVP_WRAP(o->nvabs,       NVP_ABS,       nvp_abs);
    NVP_WRAP(o->nvinv,       NVP_INV,       nvp_inv);
    NVP_WRAP(o->nvaddconst,  NVP_ADDCONST,  nvp_addconst);
    NVP_WRAP(o->nvcompare,   NVP_COMPARE,   nvp_compare);
    /* reduction (standard) */
    NVP_WRAP(o->nvdotprod,     NVP_DOTPROD,      nvp_dotprod);
    NVP_WRAP(o->nvmaxnorm,     NVP_MAXNORM,      nvp_maxnorm);
    NVP_WRAP(o->nvwrmsnorm,    NVP_WRMSNORM,     nvp_wrmsnorm);
    NVP_WRAP(o->nvwrmsnormmask,NVP_WRMSNORMMASK, nvp_wrmsnormmask);
    NVP_WRAP(o->nvmin,         NVP_MIN,          nvp_min);
    NVP_WRAP(o->nvwl2norm,     NVP_WL2NORM,      nvp_wl2norm);
    NVP_WRAP(o->nvl1norm,      NVP_L1NORM,       nvp_l1norm);
    NVP_WRAP(o->nvinvtest,     NVP_INVTEST,      nvp_invtest);
    NVP_WRAP(o->nvconstrmask,  NVP_CONSTRMASK,   nvp_constrmask);
    NVP_WRAP(o->nvminquotient, NVP_MINQUOTIENT,  nvp_minquotient);
    /* reduction (local kernels) */
    NVP_WRAP(o->nvdotprodlocal,     NVP_DOTPRODLOCAL,     nvp_dotprodlocal);
    NVP_WRAP(o->nvmaxnormlocal,     NVP_MAXNORMLOCAL,     nvp_maxnormlocal);
    NVP_WRAP(o->nvminlocal,         NVP_MINLOCAL,         nvp_minlocal);
    NVP_WRAP(o->nvl1normlocal,      NVP_L1NORMLOCAL,      nvp_l1normlocal);
    NVP_WRAP(o->nvinvtestlocal,     NVP_INVTESTLOCAL,     nvp_invtestlocal);
    NVP_WRAP(o->nvconstrmasklocal,  NVP_CONSTRMASKLOCAL,  nvp_constrmasklocal);
    NVP_WRAP(o->nvminquotientlocal, NVP_MINQUOTIENTLOCAL, nvp_minquotientlocal);
    NVP_WRAP(o->nvwsqrsumlocal,     NVP_WSQRSUMLOCAL,     nvp_wsqrsumlocal);
    NVP_WRAP(o->nvwsqrsummasklocal, NVP_WSQRSUMMASKLOCAL, nvp_wsqrsummasklocal);
    /* fused / vector-array (NULL by default; wrapped only if populated) */
    NVP_WRAP(o->nvlinearcombination,       NVP_LINEARCOMBINATION,       nvp_linearcombination);
    NVP_WRAP(o->nvscaleaddmulti,           NVP_SCALEADDMULTI,           nvp_scaleaddmulti);
    NVP_WRAP(o->nvdotprodmulti,            NVP_DOTPRODMULTI,            nvp_dotprodmulti);
    NVP_WRAP(o->nvlinearsumvectorarray,    NVP_LINEARSUMVECTORARRAY,    nvp_linearsumvectorarray);
    NVP_WRAP(o->nvscalevectorarray,        NVP_SCALEVECTORARRAY,        nvp_scalevectorarray);
    NVP_WRAP(o->nvconstvectorarray,        NVP_CONSTVECTORARRAY,        nvp_constvectorarray);
    NVP_WRAP(o->nvwrmsnormvectorarray,     NVP_WRMSNORMVECTORARRAY,     nvp_wrmsnormvectorarray);
    NVP_WRAP(o->nvwrmsnormmaskvectorarray, NVP_WRMSNORMMASKVECTORARRAY, nvp_wrmsnormmaskvectorarray);
    NVP_WRAP(o->nvdotprodmultilocal,       NVP_DOTPRODMULTILOCAL,       nvp_dotprodmultilocal);
    /* structural */
    NVP_WRAP(o->nvclone,      NVP_CLONE,      nvp_clone);
    NVP_WRAP(o->nvcloneempty, NVP_CLONEEMPTY, nvp_cloneempty);
}

/* Count how many distinct slots point at OUR shim on this table (the
 * "wrapped-entry debug count" the spec cross-checks against). */
static int nvp_count_shims_on(N_Vector v) {
    if (v == NULL || v->ops == NULL) return 0;
    N_Vector_Ops o = v->ops;
    int n = 0;
    if ((void *)o->nvlinearsum    == (void *)nvp_linearsum)    n++;
    if ((void *)o->nvconst        == (void *)nvp_const)        n++;
    if ((void *)o->nvprod         == (void *)nvp_prod)         n++;
    if ((void *)o->nvdiv          == (void *)nvp_div)          n++;
    if ((void *)o->nvscale        == (void *)nvp_scale)        n++;
    if ((void *)o->nvabs          == (void *)nvp_abs)          n++;
    if ((void *)o->nvinv          == (void *)nvp_inv)          n++;
    if ((void *)o->nvaddconst     == (void *)nvp_addconst)     n++;
    if ((void *)o->nvcompare      == (void *)nvp_compare)      n++;
    if ((void *)o->nvdotprod      == (void *)nvp_dotprod)      n++;
    if ((void *)o->nvmaxnorm      == (void *)nvp_maxnorm)      n++;
    if ((void *)o->nvwrmsnorm     == (void *)nvp_wrmsnorm)     n++;
    if ((void *)o->nvwrmsnormmask == (void *)nvp_wrmsnormmask) n++;
    if ((void *)o->nvmin          == (void *)nvp_min)          n++;
    if ((void *)o->nvwl2norm      == (void *)nvp_wl2norm)      n++;
    if ((void *)o->nvl1norm       == (void *)nvp_l1norm)       n++;
    if ((void *)o->nvinvtest      == (void *)nvp_invtest)      n++;
    if ((void *)o->nvconstrmask   == (void *)nvp_constrmask)   n++;
    if ((void *)o->nvminquotient  == (void *)nvp_minquotient)  n++;
    if ((void *)o->nvdotprodlocal     == (void *)nvp_dotprodlocal)     n++;
    if ((void *)o->nvmaxnormlocal     == (void *)nvp_maxnormlocal)     n++;
    if ((void *)o->nvminlocal         == (void *)nvp_minlocal)         n++;
    if ((void *)o->nvl1normlocal      == (void *)nvp_l1normlocal)      n++;
    if ((void *)o->nvinvtestlocal     == (void *)nvp_invtestlocal)     n++;
    if ((void *)o->nvconstrmasklocal  == (void *)nvp_constrmasklocal)  n++;
    if ((void *)o->nvminquotientlocal == (void *)nvp_minquotientlocal) n++;
    if ((void *)o->nvwsqrsumlocal     == (void *)nvp_wsqrsumlocal)     n++;
    if ((void *)o->nvwsqrsummasklocal == (void *)nvp_wsqrsummasklocal) n++;
    /* fused / vector-array (populated on some backends; NULL by default) */
    if ((void *)o->nvlinearcombination       == (void *)nvp_linearcombination)       n++;
    if ((void *)o->nvscaleaddmulti           == (void *)nvp_scaleaddmulti)           n++;
    if ((void *)o->nvdotprodmulti            == (void *)nvp_dotprodmulti)            n++;
    if ((void *)o->nvlinearsumvectorarray    == (void *)nvp_linearsumvectorarray)    n++;
    if ((void *)o->nvscalevectorarray        == (void *)nvp_scalevectorarray)        n++;
    if ((void *)o->nvconstvectorarray        == (void *)nvp_constvectorarray)        n++;
    if ((void *)o->nvwrmsnormvectorarray     == (void *)nvp_wrmsnormvectorarray)     n++;
    if ((void *)o->nvwrmsnormmaskvectorarray == (void *)nvp_wrmsnormmaskvectorarray) n++;
    if ((void *)o->nvdotprodmultilocal       == (void *)nvp_dotprodmultilocal)       n++;
    if ((void *)o->nvclone      == (void *)nvp_clone)      n++;
    if ((void *)o->nvcloneempty == (void *)nvp_cloneempty) n++;
    return n;
}

bool nvec_prof_clone_carries_shims(N_Vector v) {
    if (!nvec_prof_is_on()) return true;
    if (v == NULL) return false;

    int base = nvp_count_shims_on(v);

    N_Vector c1 = N_VClone(v);
    N_Vector c2 = N_VCloneEmpty(v);
    int n1 = nvp_count_shims_on(c1);
    int n2 = nvp_count_shims_on(c2);
    if (c1 != NULL) N_VDestroy(c1);
    if (c2 != NULL) N_VDestroy(c2);

    bool pass = (base > 0) && (n1 == base) && (n2 == base);
    fprintf(stdout,
            "[NVEC_PROF] clone-propagation smoke assert: base=%d "
            "N_VClone=%d N_VCloneEmpty=%d -> %s\n",
            base, n1, n2, pass ? "PASS" : "FAIL");
    return pass;
}

void nvec_prof_dump(const char *project_name, int NY, int nthreads,
                    const char *backend, const char *outpath) {
    if (!nvec_prof_is_on()) return;

    char path[4096];
    snprintf(path, sizeof(path), "%s/nvec_prof.csv", outpath);
    FILE *fp = fopen(path, "w");
    if (fp == NULL) {
        fprintf(stderr, "[NVEC_PROF] WARNING: cannot open %s for write; "
                        "profile dump skipped.\n", path);
        return;
    }

    /* Header lines (pinned schema): project/NY/nthreads/backend. */
    fprintf(fp, "# project_name=%s\n", project_name ? project_name : "");
    fprintf(fp, "# NY=%d\n", NY);
    fprintf(fp, "# nthreads=%d\n", nthreads);
    fprintf(fp, "# backend=%s\n", backend ? backend : "");
    fprintf(fp, "op_name,op_class,calls,total_ns\n");

    int wrapped_entries = 0;
    int invoked_entries = 0;
    for (int i = 0; i < NVP_COUNT; i++) {
        if (!g_slot[i].wrapped) continue; /* only rows for wrapped ops */
        wrapped_entries++;
        if (g_slot[i].calls > 0) invoked_entries++;
        fprintf(fp, "%s,%s,%llu,%llu\n",
                g_slot[i].name, g_slot[i].op_class,
                (unsigned long long)g_slot[i].calls,
                (unsigned long long)g_slot[i].total_ns);
    }
    fclose(fp);

    /* Cross-check line to stdout (spec "no unwrapped-op leakage"): the
     * number of wrapped table entries and how many CVODE actually
     * invoked. A reviewer compares wrapped_entries against the shim count
     * on the live table. */
    fprintf(stdout,
            "[NVEC_PROF] dump %s : wrapped_entries=%d invoked(calls>0)=%d "
            "backend=%s NY=%d nthreads=%d\n",
            path, wrapped_entries, invoked_entries,
            backend ? backend : "", NY, nthreads);
}
