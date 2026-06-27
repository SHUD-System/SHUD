/* MD_precond_identity.cpp — P8-precond-0 identity preconditioner body.
 *
 * Implements the two CVLS preconditioner callbacks declared in
 * MD_precond_identity.h. Both bodies are intentionally minimal:
 *
 *   PSetupIdentity — declares the Jacobian as NOT recomputed (set
 *     *jcurPtr = SUNFALSE). Bracketed with an RAII shud_profile::Timer
 *     bucket "t_precond_setup" so soft-gate-6 (setup-overhead evidence
 *     per spec p8precond-zero-identity-spike Scenario L108-113) surfaces
 *     via timer.cpp catch-all auto-emit (extras:t_precond_setup).
 *
 *   PSolveIdentity — solves P z = r with P = I, i.e. copies r into z
 *     via N_VScale(1.0, r, z). No timer (cheap operation; spec only
 *     requires setup-side bucket).
 *
 * Both return 0 (CVLS success). The signatures match SUNDIALS 6.0.0
 * CVLsPrecSetupFn (cvode/cvode_ls.h:57) + CVLsPrecSolveFn
 * (cvode/cvode_ls.h:61) verbatim.
 *
 * Part of p8pre-spike Step 2 P8-precond-0 (DankerMu/SHUD-OpenMP #338,
 * #345). */
#include "MD_precond_identity.h"

#include <sundials/sundials_types.h>   /* realtype, booleantype, SUNFALSE, SUN_RCONST */
#include <sundials/sundials_nvector.h> /* N_VScale */
#include <nvector/nvector_serial.h>    /* N_Vector (concrete serial) */

/* Timer include + RAII instance are #ifdef-gated to match the project
 * convention (see SHUD/src/Model/f.cpp). Under SHUD_ENABLE_PROFILE=0
 * the include path does NOT carry tools/profile (Makefile L288), so an
 * unguarded include would break PROFILE=0 builds. The header itself
 * provides a no-op Timer stub for safety, but the guard is the
 * canonical pattern. */
#ifdef SHUD_ENABLE_PROFILE
#include "timer.h"
#endif

extern "C" {

int PSetupIdentity(realtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                   booleantype jok, booleantype *jcurPtr,
                   realtype /*gamma*/, void * /*user_data*/) {
#ifdef SHUD_ENABLE_PROFILE
    /* Soft-gate-6 setup-overhead evidence. timer.cpp's emit_extra
     * catch-all surfaces this bucket under `extras:` in profile yaml
     * because it's NOT in kKnownRawOrCanonical[] (timer.cpp L188-191). */
    shud_profile::Timer _t("t_precond_setup");
#endif
    /* Mirror the SUNDIALS canonical Precond convention (see SUNDIALS
     * 6.0.0 example cvDiurnal_kry.c) so the CVLS internal accounting
     * counters (CVodeGetNumPrecEvals -> npe) increment as expected:
     *   jok == SUNTRUE  -> "Jacobian still good": *jcurPtr = SUNFALSE
     *   jok == SUNFALSE -> "rebuild Jacobian":    *jcurPtr = SUNTRUE
     *
     * P = I either way (P^{-1} z = r implemented via N_VScale in
     * PSolveIdentity), so the *jcurPtr value has no numerical effect
     * on the iteration. With this canonical pattern CVodeGetNumPrecEvals
     * returns the expected nonzero npe counter (spec gate 3
     * "nps and npe accumulation" — p8precond-zero-identity-spike). */
    (void)jok;
    *jcurPtr = jok ? SUNFALSE : SUNTRUE;
    return 0;
}

int PSolveIdentity(realtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                   N_Vector r, N_Vector z, realtype /*gamma*/,
                   realtype /*delta*/, int /*lr*/, void * /*user_data*/) {
    /* P^{-1} r = I r = r → z := 1.0 * r. */
    N_VScale(SUN_RCONST(1.0), r, z);
    return 0;
}

} /* extern "C" */
