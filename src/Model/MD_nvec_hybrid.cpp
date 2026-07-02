/* =====================================================================
 * P12-nvec PR-N1 (#443) — Config E hybrid NVector reduction overrides.
 * See MD_nvec_hybrid.hpp for the design rationale, hard rules, aliasing
 * note, and bitwise-equivalence guarantee.
 *
 * The whole translation unit is compiled ONLY under SHUD_NVEC_HYBRID, so
 * a default (Config C) build sees an empty TU → preprocessor-identical to
 * pre-change. The install/assert symbols still exist for the linker via
 * the no-op fallback at the bottom of this file so shud.cpp can call them
 * unconditionally under the `#ifdef SHUD_USE_OPENMP_NVECTOR` block.
 * ===================================================================== */

#include "MD_nvec_hybrid.hpp"

#ifdef SHUD_NVEC_HYBRID

#include <stdio.h>
#include <sundials/sundials_math.h>   /* SUNRabs, SUNRsqrt, SUNSQR, BIG_REAL */
#include <sundials/sundials_types.h>  /* realtype, sunindextype, booleantype */

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define HALF   RCONST(0.5)
#define ONEPT5 RCONST(1.5)

/* ---------------------------------------------------------------------
 * SHUD_NVEC_NOOPT — pin each reduction body to -O0 codegen.
 *
 * WHY THIS EXISTS (the G-E1 bitwise gate depends on it):
 * The vendored SUNDIALS is built with CMAKE_BUILD_TYPE="" (empty) →
 * effectively -O0; its serial reductions (N_VWrmsNorm_Serial, ...) are
 * the reference Config C runs against. SHUD itself compiles at
 * `-O2 -ffp-contract=off` (B0 IEEE-754 lockdown). Compiling the *same*
 * scalar serial loop at -O2 lets the optimizer reassociate / reschedule
 * the accumulation (SLP, strength reduction) so its rounding drifts from
 * the -O0 library fold by ~1 ULP. Measured on keliya solver data:
 *   override(-O2)=2.1536859975410953038e-05
 *   libSerial(-O0)=2.153685997541094965e-05   (eq=0, 1 ULP)
 * cascading to 96.7% of keliya.rivqdown.dat values (first divergence at
 * byte offset 3728) → G-E1 vs-Config-C would FAIL. A 3000-dataset unit
 * sweep confirms: -O2 body diverges from the library on 363/3000 inputs;
 * the -O0-pinned body diverges on 0/3000. Pinning the reduction bodies
 * to -O0 makes the SHUD-compiled scalar fold bit-match the -O0 library.
 *
 * This is NOT a spec-approach change: the bodies remain plain serial
 * loops over the generic API (N_VGetArrayPointer / N_VGetLength). The
 * attribute only removes the optimizer's freedom to reassociate FP ops.
 *
 * DOCUMENTED ASSUMPTION / fragility: correctness of the *bitwise-to-C*
 * guarantee rests on the SUNDIALS library being -O0-equivalent scalar.
 * If `./configure` is ever changed to build SUNDIALS at -O3
 * (CMAKE_BUILD_TYPE=Release → CMAKE_C_FLAGS_RELEASE=-O3), Config C's
 * reference reductions change and this -O0 pin would no longer match;
 * the G-E1 gate would catch it. The construction-guaranteed alternative
 * (copy into a Serial N_Vector and call the library N_V*_Serial) is
 * documented in the PR-N1 evidence as the robust fallback.
 * ------------------------------------------------------------------- */
#if defined(__clang__)
#  define SHUD_NVEC_NOOPT __attribute__((optnone))
#elif defined(__GNUC__)
#  define SHUD_NVEC_NOOPT __attribute__((optimize("O0")))
#else
#  define SHUD_NVEC_NOOPT
#endif

/* ---------------------------------------------------------------------
 * Serial reduction overrides — generic API ONLY.
 *
 * Every body below is a line-for-line port of the matching
 * `cvode-6.0.0/src/nvector/serial/nvector_serial.c` function, with
 * `NV_LENGTH_S`/`NV_DATA_S` replaced by the backend-agnostic
 * `N_VGetLength`/`N_VGetArrayPointer`. Reproducing the serial body (not
 * the OpenMP body) guarantees the exact left-to-right accumulation order,
 * so results are bitwise-identical to the Serial NVector backend AND
 * across every thread count (there is no parallel accumulation left).
 * Each reduction that folds FP values carries SHUD_NVEC_NOOPT (see above)
 * so the -O2 optimizer cannot reassociate the fold away from the -O0
 * library order. Pure integer/branch bodies (invtest, constrmask) do not
 * strictly need it but carry it for uniformity + defence against future
 * FP creep.
 * ------------------------------------------------------------------- */

/* a = sum(x[i]*y[i]) — mirrors N_VDotProd_Serial */
SHUD_NVEC_NOOPT static realtype shud_dotprod(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype sum = ZERO;
  realtype *xd = N_VGetArrayPointer(x);
  realtype *yd = N_VGetArrayPointer(y);
  N = N_VGetLength(x);
  for (i = 0; i < N; i++)
    sum += xd[i] * yd[i];
  return (sum);
}

/* max|x[i]| — mirrors N_VMaxNorm_Serial */
SHUD_NVEC_NOOPT static realtype shud_maxnorm(N_Vector x)
{
  sunindextype i, N;
  realtype max = ZERO;
  realtype *xd = N_VGetArrayPointer(x);
  N = N_VGetLength(x);
  for (i = 0; i < N; i++) {
    if (SUNRabs(xd[i]) > max) max = SUNRabs(xd[i]);
  }
  return (max);
}

/* sum(SUNSQR(x[i]*w[i])) — mirrors N_VWSqrSumLocal_Serial (backs wrmsnorm) */
SHUD_NVEC_NOOPT static realtype shud_wsqrsum(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum = ZERO;
  realtype prodi;
  realtype *xd = N_VGetArrayPointer(x);
  realtype *wd = N_VGetArrayPointer(w);
  N = N_VGetLength(x);
  for (i = 0; i < N; i++) {
    prodi = xd[i] * wd[i];
    sum += SUNSQR(prodi);
  }
  return (sum);
}

/* masked sum(SUNSQR(x[i]*w[i])) — mirrors N_VWSqrSumMaskLocal_Serial */
SHUD_NVEC_NOOPT static realtype shud_wsqrsummask(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  realtype sum = ZERO;
  realtype prodi;
  realtype *xd  = N_VGetArrayPointer(x);
  realtype *wd  = N_VGetArrayPointer(w);
  realtype *idd = N_VGetArrayPointer(id);
  N = N_VGetLength(x);
  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i] * wd[i];
      sum += SUNSQR(prodi);
    }
  }
  return (sum);
}

/* sqrt(wsqrsum(x,w)/N) — mirrors N_VWrmsNorm_Serial */
SHUD_NVEC_NOOPT static realtype shud_wrmsnorm(N_Vector x, N_Vector w)
{
  return (SUNRsqrt(shud_wsqrsum(x, w) / (N_VGetLength(x))));
}

/* sqrt(wsqrsummask(x,w,id)/N) — mirrors N_VWrmsNormMask_Serial */
SHUD_NVEC_NOOPT static realtype shud_wrmsnormmask(N_Vector x, N_Vector w, N_Vector id)
{
  return (SUNRsqrt(shud_wsqrsummask(x, w, id) / (N_VGetLength(x))));
}

/* min(x[i]) — mirrors N_VMin_Serial (seed x[0], loop 1..N) */
SHUD_NVEC_NOOPT static realtype shud_min(N_Vector x)
{
  sunindextype i, N;
  realtype *xd = N_VGetArrayPointer(x);
  realtype min;
  N = N_VGetLength(x);
  min = xd[0];
  for (i = 1; i < N; i++) {
    if (xd[i] < min) min = xd[i];
  }
  return (min);
}

/* sqrt(sum(SUNSQR(x[i]*w[i]))) — mirrors N_VWL2Norm_Serial */
SHUD_NVEC_NOOPT static realtype shud_wl2norm(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum = ZERO;
  realtype prodi;
  realtype *xd = N_VGetArrayPointer(x);
  realtype *wd = N_VGetArrayPointer(w);
  N = N_VGetLength(x);
  for (i = 0; i < N; i++) {
    prodi = xd[i] * wd[i];
    sum += SUNSQR(prodi);
  }
  return (SUNRsqrt(sum));
}

/* sum|x[i]| — mirrors N_VL1Norm_Serial */
SHUD_NVEC_NOOPT static realtype shud_l1norm(N_Vector x)
{
  sunindextype i, N;
  realtype sum = ZERO;
  realtype *xd = N_VGetArrayPointer(x);
  N = N_VGetLength(x);
  for (i = 0; i < N; i++)
    sum += SUNRabs(xd[i]);
  return (sum);
}

/* z[i]=1/x[i], returns SUNFALSE if any x[i]==0 — mirrors N_VInvTest_Serial */
SHUD_NVEC_NOOPT static booleantype shud_invtest(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd = N_VGetArrayPointer(x);
  realtype *zd = N_VGetArrayPointer(z);
  booleantype no_zero_found = SUNTRUE;
  N = N_VGetLength(x);
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO)
      no_zero_found = SUNFALSE;
    else
      zd[i] = ONE / xd[i];
  }
  return no_zero_found;
}

/* constraint mask — mirrors N_VConstrMask_Serial */
SHUD_NVEC_NOOPT static booleantype shud_constrmask(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  booleantype test;
  realtype *cd = N_VGetArrayPointer(c);
  realtype *xd = N_VGetArrayPointer(x);
  realtype *md = N_VGetArrayPointer(m);
  N = N_VGetLength(x);

  temp = ZERO;
  for (i = 0; i < N; i++) {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO) continue;

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i] * cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF   && xd[i] * cd[i] <  ZERO);
    if (test) temp = md[i] = ONE;
  }

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

/* min(num[i]/denom[i]) over denom!=0 — mirrors N_VMinQuotient_Serial */
SHUD_NVEC_NOOPT static realtype shud_minquotient(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce = SUNTRUE;
  sunindextype i, N;
  realtype *nd = N_VGetArrayPointer(num);
  realtype *dd = N_VGetArrayPointer(denom);
  realtype min = BIG_REAL;
  N = N_VGetLength(num);

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce)
        min = SUNMIN(min, nd[i] / dd[i]);
      else {
        min = nd[i] / dd[i];
        notEvenOnce = SUNFALSE;
      }
    }
  }
  return (min);
}

/* single-buffer multi dot product — mirrors N_VDotProdMulti_Serial:
 * dotprods[i] = sum_j x[j]*Y[i][j], each accumulated serially. */
SHUD_NVEC_NOOPT static int shud_dotprodmultilocal(int nvec, N_Vector x, N_Vector *Y, realtype *dotprods)
{
  sunindextype j, N;
  int i;
  realtype *xd = N_VGetArrayPointer(x);
  realtype *yd;
  N = N_VGetLength(x);

  if (nvec < 1) return (-1);

  for (i = 0; i < nvec; i++) {
    dotprods[i] = ZERO;
    yd = N_VGetArrayPointer(Y[i]);
    for (j = 0; j < N; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------
 * install(): overwrite each populated reduction slot (and its aliased
 * `*local` sibling) with the serial override. Element-wise slots are left
 * untouched (stay stock OpenMP). We overwrite unconditionally — the
 * OpenMP backend always populates these slots (see the PR-N1 runtime slot
 * probe + docs/p12-nvec/nvec_reduction_audit.md) — but guard on non-NULL
 * defensively so a future ops-table change cannot install a bogus pointer
 * into a NULL slot. Writing the standard AND the `*local` slot closes the
 * aliasing hazard: the stock pointers are identical, so leaving one slot
 * would keep a stock parallel body reachable.
 * ------------------------------------------------------------------- */
#define HYB_SET(FIELD, FN) do { if ((FIELD) != NULL) (FIELD) = (FN); } while (0)

void nvec_hybrid_install(N_Vector v)
{
  if (v == NULL || v->ops == NULL) return;
  N_Vector_Ops o = v->ops;

  /* standard reductions */
  HYB_SET(o->nvdotprod,      shud_dotprod);
  HYB_SET(o->nvmaxnorm,      shud_maxnorm);
  HYB_SET(o->nvwrmsnorm,     shud_wrmsnorm);
  HYB_SET(o->nvwrmsnormmask, shud_wrmsnormmask);
  HYB_SET(o->nvmin,          shud_min);
  HYB_SET(o->nvwl2norm,      shud_wl2norm);
  HYB_SET(o->nvl1norm,       shud_l1norm);
  HYB_SET(o->nvinvtest,      shud_invtest);
  HYB_SET(o->nvconstrmask,   shud_constrmask);
  HYB_SET(o->nvminquotient,  shud_minquotient);

  /* local reduction kernels (aliased to the standard pointer on the OpenMP
   * backend — must be overridden too, else the alias keeps the stock body) */
  HYB_SET(o->nvdotprodlocal,     shud_dotprod);
  HYB_SET(o->nvmaxnormlocal,     shud_maxnorm);
  HYB_SET(o->nvminlocal,         shud_min);
  HYB_SET(o->nvl1normlocal,      shud_l1norm);
  HYB_SET(o->nvinvtestlocal,     shud_invtest);
  HYB_SET(o->nvconstrmasklocal,  shud_constrmask);
  HYB_SET(o->nvminquotientlocal, shud_minquotient);
  HYB_SET(o->nvwsqrsumlocal,     shud_wsqrsum);
  HYB_SET(o->nvwsqrsummasklocal, shud_wsqrsummask);

  /* single-buffer reduction (populated by the OpenMP backend) */
  HYB_SET(o->nvdotprodmultilocal, shud_dotprodmultilocal);

  /* fused / vector-array reductions are NULL by default (SHUD never calls
   * N_VEnable*Ops) → nothing to override; HYB_SET no-ops on the NULL slots
   * if a future config enables them it must extend this list + the audit. */

  fprintf(stdout,
          "[NVEC_HYBRID] serial reduction overrides installed on ops table "
          "(dotprod/maxnorm/wrmsnorm[mask]/min/wl2norm/l1norm/invtest/"
          "constrmask/minquotient + aliased *local + wsqrsum[mask]local + "
          "dotprodmultilocal); element-wise ops stay stock-OpenMP.\n");
}

/* ---------------------------------------------------------------------
 * Clone-propagation smoke assert. Verify that N_VClone / N_VCloneEmpty
 * carry the override pointers (the generic clone copies the ops table via
 * N_VCopyOps), and that those pointers DIFFER from a fresh stock vector's
 * reduction pointers. Checks the 10 standard reduction slots.
 * ------------------------------------------------------------------- */
static int hyb_count_overrides_on(N_Vector v)
{
  if (v == NULL || v->ops == NULL) return 0;
  N_Vector_Ops o = v->ops;
  int n = 0;
  if ((void *)o->nvdotprod      == (void *)shud_dotprod)      n++;
  if ((void *)o->nvmaxnorm      == (void *)shud_maxnorm)      n++;
  if ((void *)o->nvwrmsnorm     == (void *)shud_wrmsnorm)     n++;
  if ((void *)o->nvwrmsnormmask == (void *)shud_wrmsnormmask) n++;
  if ((void *)o->nvmin          == (void *)shud_min)          n++;
  if ((void *)o->nvwl2norm      == (void *)shud_wl2norm)      n++;
  if ((void *)o->nvl1norm       == (void *)shud_l1norm)       n++;
  if ((void *)o->nvinvtest      == (void *)shud_invtest)      n++;
  if ((void *)o->nvconstrmask   == (void *)shud_constrmask)   n++;
  if ((void *)o->nvminquotient  == (void *)shud_minquotient)  n++;
  return n;
}

bool nvec_hybrid_clone_carries_overrides(N_Vector v)
{
  if (v == NULL) return false;

  int base = hyb_count_overrides_on(v);

  N_Vector c1 = N_VClone(v);
  N_Vector c2 = N_VCloneEmpty(v);
  int n1 = hyb_count_overrides_on(c1);
  int n2 = hyb_count_overrides_on(c2);

  /* dotprod slot on the clone must be the override, not the stock pointer.
   * Read the stock pointer from a fresh clone's table BEFORE we would have
   * touched it — but the clone already carries the override, so compare the
   * clone's slot against the KNOWN stock body is not directly available
   * here (stock symbols are internal to libsundials_nvecopenmp). Instead we
   * assert (a) the clone matches the override address, and (b) the override
   * address differs from what an un-overridden sibling slot would hold —
   * covered by base==10 AND n1==base AND n2==base with distinct addresses. */
  bool pass = (base == 10) && (n1 == base) && (n2 == base);

  if (c1 != NULL) N_VDestroy(c1);
  if (c2 != NULL) N_VDestroy(c2);

  fprintf(stdout,
          "[NVEC_HYBRID] clone-propagation smoke assert: base=%d "
          "N_VClone=%d N_VCloneEmpty=%d -> %s\n",
          base, n1, n2, pass ? "PASS" : "FAIL");
  return pass;
}

/* Is `fp` one of the serial reduction override addresses? Compares against
 * every function install() writes into the ops table (standard + aliased
 * *local + wsqrsum[mask]local + dotprodmultilocal). Cast through void* to
 * compare pointer identity across the distinct ops-slot signatures. */
bool nvec_hybrid_addr_is_override(void *fp)
{
  return    fp == (void *)shud_dotprod
         || fp == (void *)shud_maxnorm
         || fp == (void *)shud_wrmsnorm
         || fp == (void *)shud_wrmsnormmask
         || fp == (void *)shud_min
         || fp == (void *)shud_wl2norm
         || fp == (void *)shud_l1norm
         || fp == (void *)shud_invtest
         || fp == (void *)shud_constrmask
         || fp == (void *)shud_minquotient
         || fp == (void *)shud_wsqrsum
         || fp == (void *)shud_wsqrsummask
         || fp == (void *)shud_dotprodmultilocal;
}

#else /* !SHUD_NVEC_HYBRID — no-op fallbacks so shud.cpp links unconditionally */

#include <sundials/sundials_nvector.h>

void nvec_hybrid_install(N_Vector /*v*/) { }
bool nvec_hybrid_clone_carries_overrides(N_Vector /*v*/) { return true; }
bool nvec_hybrid_addr_is_override(void * /*fp*/) { return false; }

#endif /* SHUD_NVEC_HYBRID */
