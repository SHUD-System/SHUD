# SUNDIALS 6.0.0 OpenMP NVector — reduction audit (P12-nvec PR-N1, task 2.1)

Config E (`SHUD_NVEC_HYBRID=1`) keeps every **element-wise** op stock-OpenMP
(per-element FP order is index-fixed → bitwise across thread counts and vs the
Serial backend) and replaces every op whose stock OpenMP implementation
**accumulates across elements** (a `reduction(...)` clause, an `omp critical`
combine, or a cross-thread race write) with a SHUD-owned serial generic-API
loop. This table enumerates the FULL populated ops surface so a reviewer can
confirm no accumulating op escapes the override.

## Audit source (SOURCE TREE, not `InstallSundials/`)

- Element/reduction bodies: `cvode-6.0.0/src/nvector/openmp/nvector_openmp.c`
  — `grep -n 'reduction('` → **5 sites**: L718 `N_VDotProd_OpenMP`,
  L837 `N_VWL2Norm_OpenMP`, L863 `N_VL1Norm_OpenMP`,
  L1032 `N_VWSqrSumLocal_OpenMP`, L1060 `N_VWSqrSumMaskLocal_OpenMP`.
  Four further accumulating ops combine partials with `#pragma omp critical`
  instead of a `reduction(` clause — also thread-count-dependent and audited
  as reductions: `N_VMaxNorm_OpenMP` (L751), `N_VMin_OpenMP` (L808),
  `N_VMinQuotient_OpenMP` (L1003), and the fused `N_VDotProdMulti_OpenMP`
  (L1273). `N_VConstrMask_OpenMP` / `N_VInvTest_OpenMP` write a shared flag
  under a race but return an **order-independent boolean**; they are
  nonetheless overridden (serial) so the whole reduction class is
  SHUD-owned and no stock parallel body is reachable.
- Fused / vector-array defaults: `cvode-6.0.0/src/sundials/sundials_nvector.c`
  L103-134 — the generic empty ops set ALL fused + vector-array + local
  ops to `NULL`; `N_VNewEmpty_OpenMP` (nvector_openmp.c L119) then leaves the
  fused + vector-array slots NULL ("fused and vector array operations are
  disabled (NULL) by default") and populates ONLY the local reduction kernels
  + the single-buffer `nvdotprodmultilocal`.
- `InstallSundials/include/nvector/nvector_openmp.h` is declaration-only (no
  loop bodies) — grepping it alone yields a false-empty audit; it is NOT the
  audit source.

## Runtime slot probe (confirms what is actually populated + the aliasing)

`N_VNew_OpenMP(16, 2)` ops-table dump (scratch probe, PR-N1):

- **Populated reductions (20 slots)** — the standard slot and its `*local`
  sibling are the SAME function pointer (SUNDIALS aliases them):
  `nvdotprod`==`nvdotprodlocal`, `nvmaxnorm`==`nvmaxnormlocal`,
  `nvmin`==`nvminlocal`, `nvl1norm`==`nvl1normlocal`,
  `nvinvtest`==`nvinvtestlocal`, `nvconstrmask`==`nvconstrmasklocal`,
  `nvminquotient`==`nvminquotientlocal` (each pair identical address);
  plus `nvwrmsnorm`, `nvwrmsnormmask`, `nvwl2norm`, `nvwsqrsumlocal`,
  `nvwsqrsummasklocal`, `nvdotprodmultilocal`.
- **Populated element-wise (9 slots, stay stock)**: `nvlinearsum`, `nvconst`,
  `nvprod`, `nvdiv`, `nvscale`, `nvabs`, `nvinv`, `nvaddconst`, `nvcompare`.
- **NULL by default (not populated → not overridden)**: `nvdotprodmulti`,
  `nvlinearcombination`, `nvscaleaddmulti`, `nvlinearsumvectorarray`,
  `nvscalevectorarray`, `nvconstvectorarray`, `nvwrmsnormvectorarray`,
  `nvwrmsnormmaskvectorarray`, `nvlinearcombinationvectorarray`,
  `nvscaleaddmultivectorarray`, `nvdotprodmultiallreduce`, plus all XBraid
  buf ops and the `nvwsqrsum*` helpers exposed only through the local slots.

**Aliasing consequence (reviewer note, PR-N0):** because `nvdotprod` and
`nvdotprodlocal` are the same pointer, overriding only one slot would leave the
other pointing at the stock parallel body. `MD_nvec_hybrid` therefore writes the
serial override into BOTH the standard slot AND its `*local` sibling for every
aliased reduction — the override is idempotent per slot (it writes a SHUD
address, never reads the stock one), so no sibling is corrupted.

## Audit table

| # | op (ops slot) | stock OpenMP parallel? | `reduction(` pragma? | cross-element accumulation? | overridden (Config E)? | rationale |
|---|---|---|---|---|---|---|
| 1 | `nvlinearsum` | yes (`parallel for`) | no | no | **no** | element-wise `z[i]=a·x[i]+b·y[i]`; index-fixed FP order |
| 2 | `nvconst` | yes | no | no | **no** | element-wise `z[i]=c` |
| 3 | `nvprod` | yes | no | no | **no** | element-wise `z[i]=x[i]·y[i]` |
| 4 | `nvdiv` | yes | no | no | **no** | element-wise `z[i]=x[i]/y[i]` |
| 5 | `nvscale` | yes | no | no | **no** | element-wise `z[i]=c·x[i]` |
| 6 | `nvabs` | yes | no | no | **no** | element-wise `z[i]=|x[i]|` |
| 7 | `nvinv` | yes | no | no | **no** | element-wise `z[i]=1/x[i]` |
| 8 | `nvaddconst` | yes | no | no | **no** | element-wise `z[i]=x[i]+b` |
| 9 | `nvcompare` | yes | no | no | **no** | element-wise `z[i]=(|x[i]|≥c)?1:0` |
| 10 | `nvdotprod` | yes | **yes** (L718) | yes (Σ x·y) | **yes** | `reduction(+:sum)` → thread-count-dependent order |
| 11 | `nvdotprodlocal` | yes | **yes** (alias of #10) | yes | **yes** | SAME pointer as #10; override both slots |
| 12 | `nvmaxnorm` | yes | no (`critical`) | yes (max reduce, L751) | **yes** | per-thread `tmax` combined under `critical`; override for class-completeness (max is value-order-independent but kept serial) |
| 13 | `nvmaxnormlocal` | yes | no (alias of #12) | yes | **yes** | SAME pointer as #12 |
| 14 | `nvwrmsnorm` | yes | via `nvwsqrsumlocal` | yes (Σ (x·w)²) | **yes** | calls stock parallel `N_VWSqrSumLocal_OpenMP`; serial override computes `sqrt(Σ(x·w)²/N)` directly |
| 15 | `nvwrmsnormmask` | yes | via `nvwsqrsummasklocal` | yes | **yes** | masked variant of #14 |
| 16 | `nvmin` | yes | no (`critical`) | yes (min reduce, L808) | **yes** | per-thread `tmin` under `critical`; serial `min=x[0]; i=1..N` |
| 17 | `nvminlocal` | yes | no (alias of #16) | yes | **yes** | SAME pointer as #16 |
| 18 | `nvwl2norm` | yes | **yes** (L837) | yes (Σ (x·w)²) | **yes** | `reduction(+:sum)`; serial `sqrt(Σ(x·w)²)` |
| 19 | `nvl1norm` | yes | **yes** (L863) | yes (Σ |x|) | **yes** | `reduction(+:sum)`; serial `Σ|x[i]|` |
| 20 | `nvl1normlocal` | yes | **yes** (alias of #19) | yes | **yes** | SAME pointer as #19 |
| 21 | `nvinvtest` | yes | no (race on `val`) | flag (order-independent bool) | **yes** | `z[i]=1/x[i]`, returns "no zero found"; overridden serial for class-completeness |
| 22 | `nvinvtestlocal` | yes | no (alias of #21) | flag | **yes** | SAME pointer as #21 |
| 23 | `nvconstrmask` | yes | no (race on `temp`) | flag (order-independent bool) | **yes** | constraint mask + "any violated" flag; overridden serial for class-completeness |
| 24 | `nvconstrmasklocal` | yes | no (alias of #23) | flag | **yes** | SAME pointer as #23 |
| 25 | `nvminquotient` | yes | no (`critical`) | yes (min reduce, L1003) | **yes** | per-thread `tmin` under `critical`; serial min of `num/denom` over `denom≠0` |
| 26 | `nvminquotientlocal` | yes | no (alias of #25) | yes | **yes** | SAME pointer as #25 |
| 27 | `nvwsqrsumlocal` | yes | **yes** (L1032) | yes (Σ (x·w)²) | **yes** | `reduction(+:sum)`; serial `Σ(x·w)²` (backing kernel for #14) |
| 28 | `nvwsqrsummasklocal` | yes | **yes** (L1060) | yes | **yes** | `reduction(+:sum)`; masked (backing kernel for #15) |
| 29 | `nvdotprodmultilocal` | yes | no (`critical`) | yes (Σ per vec, L1273) | **yes** | single-buffer multi-dot; serial per-vector `Σ x·Y[k]` |
| — | `nvdotprodmulti` | (NULL) | — | — | n/a | not populated by `N_VNewEmpty_OpenMP` |
| — | `nvlinearcombination` | (NULL) | — | — | n/a | not populated (fused disabled by default) |
| — | `nvscaleaddmulti` | (NULL) | — | — | n/a | not populated |
| — | `nvlinearsumvectorarray` | (NULL) | — | — | n/a | not populated (vector-array disabled by default) |
| — | `nvscalevectorarray` | (NULL) | — | — | n/a | not populated |
| — | `nvconstvectorarray` | (NULL) | — | — | n/a | not populated |
| — | `nvwrmsnormvectorarray` | (NULL) | — | — | n/a | not populated |
| — | `nvwrmsnormmaskvectorarray` | (NULL) | — | — | n/a | not populated |

**Summary:** 9 element-wise ops stay stock; **20 populated reduction slots**
(10 standard + 9 aliased `*local` + `nvdotprodmultilocal`) are overridden with
serial generic-API loops. Every one of the 5 `reduction(` source hits maps to an
overridden slot (dotprod, wl2norm, l1norm, wsqrsumlocal, wsqrsummasklocal), and
every `critical`-combined accumulator (maxnorm, min, minquotient,
dotprodmultilocal) is overridden too. The fused/vector-array ops are NULL by
default (SHUD never calls `N_VEnable*`), so there is no reachable stock parallel
reduction after installation. The G-E1 keliya cross-thread SHA gate is the
behavioral backstop.

**Serial-equivalence note:** each override mirrors the corresponding
`nvector_serial.c` body exactly (same left-to-right accumulation, same
`min=x[0]` seed, same `SUNSQR(x·w)` per-element square, same `SUNRabs`), so
Config E is bitwise-identical to the Serial NVector backend (Config C) as well
as across thread counts.

**Codegen pin (`SHUD_NVEC_NOOPT`) — platform-dependent, needed on Apple
clang / ARM:** mirroring the serial *source* is necessary but not, on every
toolchain, sufficient. Config C's reference reductions are the **SUNDIALS
library** serial functions (`N_VDotProd_Serial`, `N_VWSqrSumLocal_Serial`, …).
Whether a SHUD-compiled generic-API serial loop bit-matches them depends on the
compiler, because SHUD's `-ffp-contract=off` (B0 IEEE-754 lockdown) forces two
codegen choices on the override that the vendored library does not make:

1. **FMA contraction.** The vendored library is compiled with **no
   `-ffp-contract` flag** (`InstallSundials/.../flags.make` carries no contract
   flag; `CMAKE_BUILD_TYPE=""`), so its *default* contraction fuses
   `sum += x[i]*y[i]` into a single-rounding `fma` on targets that have FMA.
   `-ffp-contract=off` **strips** that FMA from our override → two roundings
   instead of one.
2. **Auto-vectorization.** At `-O2` the non-FMA reduction is a vectorizer
   candidate; the loop/SLP vectorizer folds partial sums across SIMD lanes → a
   different (pairwise/tree) order than the library's sequential scalar loop.

The mechanism was NOT reassociation of the scalar loop (IEEE FP-add reassoc is
illegal at plain `-O2` without `-ffast-math`, and neither clang nor gcc does
it — a plain `-O2 -ffp-contract=off` loop matches an `-O0`-same-flags loop 0/N).
It is the FMA-state + vectorization difference above, and it is
**platform-specific**:

| toolchain | library codegen | plain override (`-O2 -ffp-contract=off`) | match? |
|-----------|-----------------|-------------------------------------------|--------|
| Apple clang / ARM (Mac) | scalar **`fmadd`** | non-FMA, `-O2`-vectorized | **NO** — dot 4785/5000, wsqrsum 1819/5000 diverge (~1 ULP) |
| x86_64 / gcc (server + CI) | scalar non-FMA (`mulsd`+`addsd`) | scalar non-FMA (same) | **YES** — 0/5000, attribute not needed |

On clang, `-ffp-contract=on` **alone** does not fix it (the loop still
vectorizes → 1823/5000); only scalar **and** FMA together match. **Fix:** each
FP-folding override carries `SHUD_NVEC_NOOPT`. On clang it expands to
`__attribute__((optnone))`, which disables vectorization AND discards the
function-level `-ffp-contract=off` (restoring default contraction → FMA), so the
override reproduces the library's scalar-FMA fold exactly (**0/5000**). On gcc
it expands to `__attribute__((optimize("O0","no-tree-vectorize")))` — scalar +
default contraction; not required for correctness on x86 (the plain override
already matches, **0/5000** measured on server gcc-13) but pins the intent and
guards a future gcc that might vectorize this loop.

Spec mechanism is intact: the bodies are plain serial loops over the generic
API (`N_VGetArrayPointer` / `N_VGetLength`), NO `N_V*_Serial` call, NO
content-struct macro; the attribute only constrains codegen (scalar + default
contraction), not the source logic. **Honest fragility:** the bitwise-to-C
guarantee rests on a per-compiler codegen coincidence (clang's optnone-restores-
contraction side effect; gcc's `-O2 -ffp-contract=off` already yielding the
library fold) — NOT on the library's `-O` level. The G-E1 SHA gate (Mac clang +
PR-N2 server gcc matrix + CI GCC) is the backstop for any future toolchain
change. Evidence: `.review-evidence/p12-nvec/pr-n1/optnone_fold_order.evidence.txt`
(per-variant divergence + FMA disassembly, both toolchains) and
`.../gcc_spot_leg/` (server gcc G-E1 heihe spot verdict).
