/* tests/s1d_configd_nvec_smoke.cpp
 *
 * Runtime smoke test for the plain SUNDIALS OpenMP N_Vector, i.e. the
 * build variant `SHUD_ENABLE_OPENMP_RHS=1 + SHUD_USE_OPENMP_NVECTOR=1`
 * without the SHUD serial-sum overrides (`SHUD_NVEC_HYBRID`). That
 * variant gives thread-count-dependent results and is blocked by the
 * Makefile for the model binary; this probe only checks the vector
 * backend itself.
 *
 * What it asserts:
 *   - `N_VNew_OpenMP` actually constructs an OpenMP-backed N_Vector
 *     (`N_VGetVectorID(v) == SUNDIALS_NVEC_OPENMP`).
 *   - The generic `N_VDestroy(v)` correctly dispatches via the
 *     vector's `ops->nvdestroy` slot. SHUD destroys every vector
 *     through this generic entry, because calling
 *     `N_VDestroy_Serial(v)` on an OpenMP-backed vector is a type-tag
 *     mismatch and undefined behavior. This test exercises that path
 *     on the only backend where it matters.
 *
 * What it does NOT exercise:
 *   - RHS evaluation under the OpenMP N_Vector backend.
 *   - Multi-thread correctness / scalability.
 *   - SHUD framework integration (the Makefile target deliberately
 *     does NOT link `$(SHUD_SRC_NOMAIN)`; only the SUNDIALS OpenMP
 *     N_Vector and the generic destroy are validated).
 *
 * Build and run: `make smoke_configd`.
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

    /* Generic destroy. If this dispatched to the wrong backend's
     * destroy (e.g. the `_Serial` flavor), AddressSanitizer / valgrind
     * would catch it; in a release build we at least exercise the
     * code path. */
    N_VDestroy(v);
    SUNContext_Free(&ctx);

    std::printf("OK: N_VGetVectorID == SUNDIALS_NVEC_OPENMP\n");
    return 0;
}
