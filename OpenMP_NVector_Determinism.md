# Determinism of the OpenMP vector layer

This note is for maintainers. It explains why the default `make shud_omp`
build gives bit-identical output at every thread count and bit-identical
output to the serial-vector build (`SHUD_USE_OPENMP_NVECTOR=0`), and which
assumption about the compiler this depends on. Read it before changing the
compiler, the compiler flags, or the SUNDIALS version.

The implementation is in `src/Model/MD_nvec_hybrid.cpp`.

## Principle

CVODE works on vectors through a table of function pointers (the `N_Vector`
operations table). SUNDIALS ships an OpenMP implementation of that table.
Its operations fall into two classes:

- **Element-wise** operations (`z[i] = a*x[i] + b*y[i]`, scaling, …). Each
  output entry depends only on the input entries with the same index, so
  the result does not depend on how the loop is split between threads.
- **Accumulating** operations (dot product, norms, minimum, …). The stock
  OpenMP versions combine per-thread partial results, so the order of the
  floating-point additions — and therefore the last bits of the result —
  changes with the thread count. CVODE uses these values in its step-size
  and convergence tests, so the differences grow into visibly different
  trajectories.

The default build (`SHUD_NVEC_HYBRID=1`, reported at startup as "Config E")
keeps the stock OpenMP element-wise operations and replaces every
accumulating operation with a serial loop owned by SHUD. The fastest build
(`SHUD_NVEC_DETRED=1`, "Config E2") replaces the summations with a parallel
sum over fixed blocks of `SHUD_NVEC_DETRED_B` entries (4096 by default)
combined in a fixed binary tree. The tree depends only on the vector length
and the block size, never on the thread count, so the result is again
independent of the thread count — but it is a different summation order,
hence not bit-identical to the serial-vector build.

Using the OpenMP vector backend without these replacements
(`SHUD_NVEC_HYBRID=0`) gives output that changes with the thread count. The
Makefile refuses that combination unless `SHUD_ALLOW_CONFIG_D=1` is given.

## Which operations are replaced

The audit below was made against the SUNDIALS 6.0.0 sources, not the
installed headers (the header `nvector_openmp.h` only declares the
functions; the loops are in the `.c` file):

- `cvode-6.0.0/src/nvector/openmp/nvector_openmp.c` — five operations use an
  OpenMP `reduction(` clause: `N_VDotProd_OpenMP` (line 718),
  `N_VWL2Norm_OpenMP` (837), `N_VL1Norm_OpenMP` (863),
  `N_VWSqrSumLocal_OpenMP` (1032), `N_VWSqrSumMaskLocal_OpenMP` (1060).
  Four more combine per-thread partial results inside `#pragma omp critical`
  and are just as order-dependent: `N_VMaxNorm_OpenMP` (751),
  `N_VMin_OpenMP` (808), `N_VMinQuotient_OpenMP` (1003),
  `N_VDotProdMulti_OpenMP` (1273). `N_VConstrMask_OpenMP` and
  `N_VInvTest_OpenMP` write a shared flag without synchronization; the
  returned boolean does not depend on the order, but they are replaced too,
  so that no stock parallel accumulation remains reachable.
- `cvode-6.0.0/src/sundials/sundials_nvector.c`, lines 103–134 — the fused
  and vector-array operations are `NULL` by default. SHUD never calls
  `N_VEnable*`, so they stay `NULL` and need no replacement.

**Aliased slots.** In the table created by `N_VNew_OpenMP`, a standard
operation and its `*local` sibling are the same function pointer
(`nvdotprod` and `nvdotprodlocal`, `nvmaxnorm` and `nvmaxnormlocal`,
`nvmin`, `nvl1norm`, `nvinvtest`, `nvconstrmask`, `nvminquotient` likewise).
Replacing only one of the two would leave the other pointing at the stock
parallel loop, so `MD_nvec_hybrid.cpp` writes the replacement into both.

**Order of installation.** The replacements are installed on the two
vectors SHUD creates (`udata`, `du`) immediately after `N_VNew_OpenMP` and
before `CVodeInit`. CVODE creates its internal work vectors with `N_VClone`,
which copies the operations table, so the replacements reach all of them.

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

In total, 9 element-wise operations keep the stock OpenMP implementation
and 20 populated accumulating slots (10 standard, 9 aliased `*local`, and
`nvdotprodmultilocal`) are replaced.

## The compiler-dependent part

Each replacement repeats the source of the corresponding function in
SUNDIALS' `nvector_serial.c`: the same left-to-right accumulation, the same
`min = x[0]` start value, the same per-entry square. That is necessary for a
bit-identical result, but on some toolchains it is not sufficient. The
serial-vector build calls the functions compiled into the **SUNDIALS
library**; the default build calls loops compiled as part of **SHUD**. The
two are compiled with different flags:

1. **Fused multiply-add (FMA).** SUNDIALS is built without an
   `-ffp-contract` flag, so on processors with FMA the compiler may turn
   `sum += x[i]*y[i]` into one fused operation with a single rounding. SHUD
   is built with `-ffp-contract=off`, which forbids this and gives two
   roundings.
2. **Vectorization.** At `-O2` the compiler may vectorize the SHUD loop and
   add partial sums across SIMD lanes, which is a different order from the
   scalar loop in the library.

Measured on 5000 random vectors, comparing the SHUD loop compiled with the
plain project flags against the library function:

| Toolchain | Library code | SHUD loop, plain `-O2 -ffp-contract=off` | Bit-identical? |
|---|---|---|---|
| Apple clang, ARM (macOS) | scalar, fused (`fmadd`) | not fused, vectorized | **No** — dot product differs in 4785 of 5000 cases, weighted square sum in 1819 of 5000 (about 1 ULP) |
| GCC 13, x86_64 (Linux) | scalar, not fused (`mulsd` + `addsd`) | scalar, not fused | **Yes** — 0 of 5000 |

On clang, `-ffp-contract=on` alone does not help (the loop is still
vectorized: 1823 of 5000 differ). The loop must be scalar **and** fused.

**The fix** is the `SHUD_NVEC_NOOPT` attribute carried by every replacement
that accumulates floating-point values:

- clang: `__attribute__((optnone))`. This switches off vectorization and
  also discards the function-level `-ffp-contract=off`, which restores the
  default contraction and thus the fused operation. Result: 0 of 5000
  differ.
- GCC: `__attribute__((optimize("O0","no-tree-vectorize")))`. Not needed for
  correctness on x86_64 today (the plain loop already matches), but it
  guards against a future GCC that vectorizes the loop.

**What this means for maintenance.** The bit-identity between the default
build and the serial-vector build holds because, on each tested compiler,
the SHUD loop and the library function happen to compile to the same
sequence of floating-point operations. It is not guaranteed by the language
or by the flags. It has been verified on Apple clang / ARM (macOS) and on
GCC 13 / x86_64 and GCC 13 / aarch64 (Linux). After any change of compiler, compiler version, compiler
flags, target architecture or SUNDIALS build options, verify it again:
build `make shud_omp` and `make shud_omp SHUD_USE_OPENMP_NVECTOR=0`, run the
same project with both at several thread counts, and compare the checksums
of the output files. They must all be equal.

Bit-identity across thread counts within one build does not depend on any
of this; it follows from the fixed order of the serial loops (default build)
or of the fixed tree (fastest build).

## SUNDIALS versions

The audit above was made on 6.0.0. Versions 6.1.1, 6.4.1 and 6.7.0 were
checked as well: `N_VNewEmpty_OpenMP` fills the same slots plus
`nvgetlocallength`, `nvprint` and `nvprintfile`, none of which accumulates,
and the number of `reduction(` and `omp critical` sites in
`nvector_openmp.c` is unchanged. With each of these versions the `ccw`
example (10 days, Apple clang / ARM) gives output files identical to those
of 6.0.0, and identical between the serial build and the OpenMP builds at
2 and 8 threads. A new major version needs the audit to be repeated.

## Records

The measurements quoted here, the disassembly, and the validation runs are
kept in the development repository:
<https://github.com/DankerMu/SHUD-OpenMP> (directory `docs/p12-nvec/` and
`docs/adr/0011-p12-nvec-tier1-verdict-and-tier2-gate.md`, tag
`cpu-accel-v1.1.1`).
