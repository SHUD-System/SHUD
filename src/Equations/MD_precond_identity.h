/* MD_precond_identity.h — P8-precond-0 identity preconditioner stub.
 *
 * Declares PSetupIdentity / PSolveIdentity satisfying SUNDIALS 6.0.0
 * CVLsPrecSetupFn (cvode/cvode_ls.h:57) + CVLsPrecSolveFn
 * (cvode/cvode_ls.h:61). The pair installs P = I via
 * CVodeSetPreconditioner so SPGMR (PREC_LEFT) wires the preconditioner
 * call path without changing the iteration numerics (P^{-1} = I).
 *
 * Purpose: P8-precond-0 zero-overhead identity spike — exercises CVLS
 * preconditioner registration + nps/npe stat accumulation while
 * preserving B1b bitwise neutrality. Future P8 stages replace the
 * identity body with a real preconditioner (block-Jacobi / KLU).
 *
 * Part of p8pre-spike Step 2 P8-precond-0 (DankerMu/SHUD-OpenMP #338,
 * #345). */
#ifndef MD_PRECOND_IDENTITY_H
#define MD_PRECOND_IDENTITY_H

#include <sundials/sundials_types.h>     /* realtype, booleantype */
#include <nvector/nvector_serial.h>      /* N_Vector */

#ifdef __cplusplus
extern "C" {
#endif

/* CVLsPrecSetupFn signature — SUNDIALS 6.0.0 cvode_ls.h:57. */
int PSetupIdentity(realtype t, N_Vector y, N_Vector fy,
                   booleantype jok, booleantype *jcurPtr,
                   realtype gamma, void *user_data);

/* CVLsPrecSolveFn signature — SUNDIALS 6.0.0 cvode_ls.h:61. */
int PSolveIdentity(realtype t, N_Vector y, N_Vector fy,
                   N_Vector r, N_Vector z, realtype gamma,
                   realtype delta, int lr, void *user_data);

#ifdef __cplusplus
}
#endif

#endif /* MD_PRECOND_IDENTITY_H */
