/* tests/s1d_configd_nvec_smoke.cpp -- S1d.2 (openMP #49).
 *
 * Runtime smoke for Config D
 * (`SHUD_ENABLE_OPENMP_RHS=1 + SHUD_USE_OPENMP_NVECTOR=1`).
 *
 * What it asserts (the only thing that matters here):
 *   - `N_VNew_OpenMP` actually constructs an OpenMP-backed N_Vector
 *     (`N_VGetVectorID(v) == SUNDIALS_NVEC_OPENMP`).
 *   - The generic `N_VDestroy(v)` correctly dispatches via the
 *     vector's `ops->nvdestroy` slot. This validates the §4.19 fix
 *     landed in #48: pre-#48 SHUD called `N_VDestroy_Serial(v)`
 *     unconditionally, which on an OpenMP-backed vector triggers a
 *     type-tag mismatch and is undefined behavior. The fix routes
 *     every destroy through the generic entry; this test exercises
 *     that path on the only backend where it matters.
 *
 * What it does NOT exercise (out of scope -- covered elsewhere or
 * not yet implemented):
 *   - RHS evaluation under OpenMP NVector backend (S2 scope).
 *   - Multiple-threads correctness / scalability (P-phase scope).
 *   - SHUD framework integration (separate from this standalone
 *     smoke -- the Makefile target deliberately does NOT link
 *     `$(SHUD_SRC_NOMAIN)` because we are only validating SUNDIALS
 *     OpenMP NVector + the generic destroy contract).
 *
 * Build / run via `make smoke_configd` (Makefile target).
 */

#include <sundials/sundials_context.h>
#include <sundials/sundials_nvector.h>
#include <nvector/nvector_openmp.h>

#include <cstdio>
#include <cstdlib>

int main(void) {
    /* SUNDIALS 6.0 mandates an explicit context for every N_Vector
     * allocator. NULL comm is legal in serial / OpenMP-only builds
     * (no MPI). */
    SUNContext ctx = nullptr;
    int flag = SUNContext_Create(NULL, &ctx);
    if (flag != 0 || ctx == nullptr) {
        std::fprintf(stderr,
                     "FAIL: SUNContext_Create returned %d (ctx=%p)\n",
                     flag, (void *)ctx);
        return 1;
    }

    /* Length 10 + 4 threads is arbitrary; the test only depends on
     * the type tag, not the contents or thread count. */
    N_Vector v = N_VNew_OpenMP(10, 4, ctx);
    if (v == nullptr) {
        std::fprintf(stderr, "FAIL: N_VNew_OpenMP returned NULL\n");
        SUNContext_Free(&ctx);
        return 1;
    }

    N_Vector_ID id = N_VGetVectorID(v);
    if (id != SUNDIALS_NVEC_OPENMP) {
        std::fprintf(stderr,
                     "FAIL: N_VGetVectorID returned %d, expected %d "
                     "(SUNDIALS_NVEC_OPENMP)\n",
                     (int)id, (int)SUNDIALS_NVEC_OPENMP);
        N_VDestroy(v);
        SUNContext_Free(&ctx);
        return 1;
    }

    /* Generic destroy: the §4.19 contract. If this dispatches to the
     * wrong backend's destroy (e.g. the pre-#48 `_Serial` flavor),
     * AddressSanitizer / valgrind would catch it; in a release build
     * we at least exercise the code path that was previously UB. */
    N_VDestroy(v);
    SUNContext_Free(&ctx);

    std::printf("OK: N_VGetVectorID == SUNDIALS_NVEC_OPENMP\n");
    return 0;
}
