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
 * SHUD_NVEC_NOOPT — force each reduction body to a scalar, FMA-preserving
 * fold that bit-matches the vendored SUNDIALS *_Serial reductions.
 *
 * WHY THIS EXISTS (the G-E1 bitwise gate needs it on Apple clang / ARM):
 * Config C's reference reductions are the SUNDIALS library serial
 * functions (N_VDotProd_Serial / N_VWSqrSumLocal_Serial / ...). Whether a
 * SHUD-compiled generic-API serial loop bit-matches them is PLATFORM-
 * DEPENDENT, because the divergence is driven by two codegen choices that
 * SHUD's `-ffp-contract=off` (B0 IEEE-754 lockdown) forces on us:
 *
 *   (1) FMA contraction. The vendored library carries no -ffp-contract
 *       flag, so its default contraction produces a fused multiply-add for
 *       `sum += x[i]*y[i]` on targets that have FMA. SHUD's
 *       -ffp-contract=off STRIPS that FMA from our override, changing the
 *       per-iteration rounding (x*y rounded, then +sum rounded — two
 *       roundings vs the library's single fma rounding).
 *   (2) Auto-vectorization. At -O2 the non-FMA reduction is a candidate
 *       for the loop / SLP vectorizer, which folds partial sums across
 *       SIMD lanes → a different (pairwise/tree) summation order than the
 *       library's sequential scalar loop.
 *
 * Measured (unit sweep, keliya/heihe-shaped data, override vs library):
 *   - Apple clang / ARM (Mac): the vendored lib emits scalar `fmadd`;
 *     our -ffp-contract=off override emits non-FMA and (at -O2) vectorizes
 *     → dot diverges 4785/5000, wsqrsum 1819/5000 (~1 ULP). BOTH knobs
 *     matter: -ffp-contract=on ALONE still diverges (still vectorized,
 *     1823/5000); scalar+FMA together (or optnone) → 0/5000.
 *   - x86_64 / gcc (server + CI): the gcc-built lib emits scalar non-FMA
 *     (`mulsd`+`addsd`); gcc at -O2 -ffp-contract=off produces the SAME
 *     scalar non-FMA fold → the plain override ALREADY matches, 0/5000,
 *     with or without this attribute (verified by disassembly + sweep on
 *     the server gcc-13 toolchain, .review-evidence/.../gcc_spot_leg/).
 *
 * THE FIX. `optnone` on clang does BOTH: it disables vectorization AND
 * discards the function-level -ffp-contract=off (restoring the default
 * contraction → FMA), so the override reproduces the library's scalar-FMA
 * fold exactly (0/5000). On gcc `optimize("O0")` disables vectorization
 * and keeps default contraction; it is not *needed* there (the plain
 * override already matches) but is harmless and keeps one portable knob.
 * So SHUD_NVEC_NOOPT is REQUIRED on clang/ARM and correctness-INERT on
 * gcc/x86 — safe on both, validated on both.
 *
 * This is NOT a spec-approach change: the bodies remain plain serial loops
 * over the generic API (N_VGetArrayPointer / N_VGetLength), with NO
 * N_V*_Serial call and NO content-struct macro. The attribute only
 * constrains codegen (scalar + default contraction), not the source logic.
 *
 * FRAGILITY (the honest statement): correctness of the bitwise-to-Config-C
 * guarantee rests on a per-compiler codegen coincidence — on clang the
 * `optnone`-restores-contraction side effect, on gcc the fact that -O2
 * -ffp-contract=off already yields the library's scalar non-FMA fold. It
 * does NOT rest on the library being rebuilt at a different -O level. If a
 * future toolchain changes either the library's contraction default or
 * clang's optnone contraction behaviour, the G-E1 SHA gate (Mac clang +
 * the PR-N2 server gcc matrix + CI GCC) is the backstop that catches the
 * regression. The copy-into-Serial + library-call alternative is
 * spec-PROHIBITED here (no N_V*_Serial in overrides), so this codegen pin
 * is the spec-compliant path.
 * ------------------------------------------------------------------- */
#if defined(__clang__)
/* clang: optnone disables vectorization AND drops the fn-level
 * -ffp-contract=off, restoring the default contraction (FMA) so the fold
 * matches the vendored *_Serial reductions on ARM. */
#  define SHUD_NVEC_NOOPT __attribute__((optnone))
#elif defined(__GNUC__)
/* gcc: -O0 + explicit no-tree-vectorize → scalar sequential fold with
 * default contraction. Not required for correctness on x86 (the plain
 * override already matches the gcc-built lib) but pins the intent and
 * guards a future gcc that might vectorize this loop at -O2. */
#  define SHUD_NVEC_NOOPT __attribute__((optimize("O0", "no-tree-vectorize")))
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
 * the OpenMP body) gives the same sequential left-to-right accumulation as
 * the library *_Serial reductions, so results are bitwise-identical to the
 * Serial NVector backend (Config C) AND across every thread count (there is
 * no parallel accumulation left). Each reduction that folds FP values
 * carries SHUD_NVEC_NOOPT (see above) so the compiler keeps a scalar,
 * default-contraction (FMA-preserving) fold that matches the vendored
 * library codegen — on clang/ARM this is REQUIRED (else -ffp-contract=off
 * strips FMA and -O2 vectorizes the loop → ~1 ULP drift); on gcc/x86 it is
 * inert (the plain override already matches). Pure integer/branch bodies
 * (invtest, constrmask) carry no FP fold but keep the attribute for
 * uniformity + defence against future FP creep.
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
