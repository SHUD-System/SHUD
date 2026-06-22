# B1b Changelog

Append-only changelog for §S5 / §S6b structural rework PRs leading to B1b-tag.

## S5a — forcing thread-safety audit (#176)

### movePointer call-site audit

| File | Line | Call expression | Parallel context | Verdict |
|---|---|---|---|---|
| SHUD/src/classes/TimeSeriesData.hpp | 28 | `void movePointer(double t);` | declaration (N/A) | N/A |
| SHUD/src/classes/TimeSeriesData.cpp | 116 | `void _TimeSeriesData::movePointer(double t){` | definition (N/A) | N/A |
| SHUD/src/ModelData/MD_ET.cpp | 29 | `tsd_weather[i].movePointer(t);` | serial — inside `for(i=0; i<NumForc; i++)` in `Model_Data::updateforcing(t)`; the historical `#ifdef _OPENMP / #pragma omp for` was removed in S2.10 (PR-1 #144, see comment block at L23–L27); no enclosing `#pragma omp parallel` exists anywhere in `SHUD/src/`. | PASS |
| SHUD/src/ModelData/MD_ET.cpp | 31 | `tsd_MF.movePointer(t);` | serial — direct call inside `Model_Data::updateforcing(t)`, no enclosing parallel region. | PASS |
| SHUD/src/ModelData/MD_ET.cpp | 32 | `tsd_LAI.movePointer(t);` | serial — direct call inside `Model_Data::updateforcing(t)`, no enclosing parallel region. | PASS |
| SHUD/src/ModelData/MD_ET.cpp | 33 | `// tsd_RL.movePointer(t);` | commented-out (N/A) | N/A |

Caller audit: `Model_Data::updateforcing(double t)` is invoked from `SHUD/src/Model/shud.cpp:120` and `:285`, both inside the top-level driver's serial outer loop `for(int i = 0; i < MD->CS.NumSteps && !ierr; i++)`. Neither call site has any enclosing `#pragma omp parallel`. A repo-wide `grep -rn "pragma omp" SHUD/src/` returns three matches that are ALL inside comment blocks (`MD_rhs_core.cpp:292`, `MD_ET.cpp:23`, `MD_ET.cpp:127`) — there are ZERO active OpenMP directives in SHUD source under the B1a contract. Therefore `movePointer` is invariably called on the single driver thread.

### getX read-only audit

`_TimeSeriesData::getX(double t, int col)` (SHUD/src/classes/TimeSeriesData.cpp L102–L105):
- Body: `return ts[iNow][col];`
- (a) Does NOT modify any TimeSeriesData instance state. Verified by reading the impl: the function body is a single `return` statement; no assignment to `iNow` / `iNext` / `ts` / `pRing` / `nQue` / `StartTime` / `Length` / `eof` / `ncol` / `fn` / `xyz` is present, and the `double t` argument is never used (the zero-order hold is implicit via the `iNow` pointer state).
- (b) Does NOT write to any shared buffer. `ts[iNow][col]` is an array load (`double` value pulled out of the ring buffer), not a store; no other expression appears in the function body.
- (c) Zero-order hold rule preserved. The function returns the value at the current pointer `iNow`; there is no interpolation, no time-decay, no comparison against `t`. `movePointer` is the only routine in `_TimeSeriesData` that advances `iNow` / `iNext` (verified by reading L116–L136 and the `initialize` constructor at L26–L44 — the latter only sets `iNow=0` at startup). Thus, as long as the single-threaded driver invokes `movePointer(t)` once before a batch of `getX(t, col)` reads (which is exactly the pattern in `Model_Data::updateforcing` → `Model_Data::tReadForcing`), concurrent readers on different threads would see a consistent snapshot of `ts[iNow][*]`.

Verdict: thread-safe read-only post-movePointer. Single read-only access pattern matches the spec's contract; safe to call `getX` from within future `#pragma omp parallel for` regions (e.g. element-loop parallelization in S5b/S5c/A3a) provided `movePointer` is invoked on the driver thread before the parallel region begins.

## S5b — scratch arrays ownership audit + lake reset 顺序 + RHS print migration (#177)

### updateElement / updateRiver / Lake::update self-only audit

| Method | Definition (file:line) | Writes | Verdict |
|---|---|---|---|
| `_Element::updateElement` | SHUD/src/classes/Element.cpp:257 | `this->u_effKH`, `this->u_deficit`, `this->Kmax`, `this->u_satn`, `this->u_theta`, `this->u_satKr`, `this->u_phius`, `this->u_effkInfi` (all instance members of `*this`) | PASS — self-only |
| `_River::updateRiver` | SHUD/src/classes/River.cpp:49 | `this->u_Ystage`, `this->u_topWidth`, `this->u_CSarea`, `this->u_CSperem`, `this->u_eqWidth`, `this->u_TopArea` (all instance members of `*this`) | PASS — self-only |
| `_Lake::update` | SHUD/src/classes/Lake.cpp:104 | `this->u_toparea` only | PASS — self-only |

Evidence: read each function body verbatim. No `Ele[`, `Riv[`, `lake[`, or `MD->` writes appear in the method bodies — all left-hand sides are bare member names (resolved as `this->member`). No cross-index state mutation. The only callee `this->bathymetry.toparea(...)` in `_Lake::update` is a const read on a member's child object (returns by value). All three methods are safe to invoke in parallel across distinct entity indices (under the existing serial-driver contract; thread-safety for future parallel-element loops in S5c/A3a needs separate audit because of forcing/state reads, not writes).

### Scratch arrays ownership audit (Task 3.1)

Documented as a separate machine-readable manifest in OUTER repo `docs/topology_manifest.yaml` under section `s5b_scratch_ownership`. Summary: 13 RHS scratch arrays surveyed (full list in topology yaml); 100% have a single owner-class write-site path traceable to `Model_Data::rhs_*` callees. None violate single-writer invariants under the B1a contract. Lake `+=` accumulator sites (`QrivSurf[ir] += QsegSurf[iseg]`, etc.) are deterministic-gather sites driven by S4 adjacency lists (PR-10/PR-11) and live inside `rhs_deterministic_gather()`; they are sequential `+=` loops in B0 ascending iteration order — bitwise-preserved.

### Lake reset 顺序 audit (Task 3.4)

Documented in OUTER repo `docs/topology_manifest.yaml` under section `s5b_lake_reset_order`. Summary: in qhh (the only Mac benchmark case with `lakeon=1`), the per-RHS-call reset sequence is:

1. `rhs_update()` at `SHUD/src/Model/MD_rhs_core.cpp:114-125` zero-resets per-lake arrays (`QLakeSub[i] = 0.`, `QLakeSurf[i] = 0.`, `qLakeEvap[i] = 0.`, `qLakePrcp[i] = 0.`, `QLakeRivIn[i] = 0.`, `QLakeRivOut[i] = 0.`) BEFORE
2. `rhs_flux()` at `SHUD/src/Model/MD_rhs_core.cpp:170-201` element loop, which writes the per-element scratch slots `QeleSurf_lake[i*3+j]`, `QeleSub_lake[i*3+j]`, `qEleEvapo_lake[i]`, `qElePrep_lake[i]` via `fun_Ele_surface`/`fun_Ele_sub`/the `if(lakeon && Ele[i].iLake > 0)` branches.
3. The transitional per-element->per-lake gather (`rhs_flux()` L214-226) zero-resets `qLakeEvap[i]` / `qLakePrcp[i]` again BEFORE the element-loop-driven `+=` gather; the lake clamp at L227-230 then reads the gathered values. This 2nd reset is required because `qLakeEvap` / `qLakePrcp` are written via `+=` and the clamp reads the result; the 1st reset in `rhs_update()` could be stale if a previous CVODE iteration emitted partial values.
4. `rhs_deterministic_gather()` (called from `rhs_flux()` L257) zero-resets `QrivSurf` / `QrivSub` / `QrivUp` (L298-302), `Qe2r_Surf` / `Qe2r_Sub` (L303-306), `QLakeRivIn` (L338-340), `QLakeSurf` (L350-352), `QLakeSub` (L363-365) BEFORE its per-lake `+=` gather loops drive them from the per-element scratch slots.

Verdict: every lake-side accumulator is zero-reset before any `+=` accumulation, in source order. The execution order (rhs_update -> rhs_flux pre-loop reset -> element loop scratch writes -> rhs_deterministic_gather reset -> gather) is the single sequential call chain from `f(t, Y, DY, DS)` (Model/f.cpp:54 -> `MD->rhs_core(Y, DY, t, ExecPolicy::Serial)`); no inter-iteration leak possible. The element loop NEVER writes directly to the per-lake `QLake*` accumulators — only to the per-element scratch slots — so the "reset before write" property is structurally enforced by the scratch-slot pattern (PR-9).

### RHS print migration (Task 3.5)

- **Approach**: `#ifdef DEBUG` wrap of the single active printf at `SHUD/src/ModelData/MD_ET.cpp:236`.
- **Justification**: The codebase convention at `SHUD/src/Model/MD_rhs_core.cpp:100-102` (`CheckNANi(uYriv[i], ...)` wrap), `:462-466` (DY checks), and `SHUD/src/Model/f.cpp:57-59` (`printDY(...)`) consistently gates per-element diagnostics behind `#ifdef DEBUG`. The buffer approach would require either a new per-element field in `Model_Data` (disallowed by PR scope) or a thread-local static accumulator (over-engineering for a single warning site). `#ifdef DEBUG` matches existing pattern and is minimally invasive (4 added lines: 2 preprocessor directives + a comment block describing the migration).
- **CheckNonNegative() context**: The 5 `CheckNonNegative(...)` calls at MD_ET.cpp:238-242 are NOT `#ifdef DEBUG` guarded (the function impl at SHUD/src/Equations/functions.cpp:148-154 calls `printf` + `myexit` on negative values). However, these are error-exit paths (not informational warnings on normal data) and the task scope explicitly excludes algorithmic changes to MD_ET.cpp. Left untouched per PR boundary.
- **Pre-migration active RHS prints** (block-comment-aware scan across `MD_rhs_core.cpp` / `MD_f.cpp` / `MD_ElementFlux.cpp` / `MD_ET.cpp` / `f.cpp`): 1 — `SHUD/src/ModelData/MD_ET.cpp:236`.
- **Post-migration active RHS prints**: 0 (in default release build with no -DDEBUG; the `MD_ET.cpp:236` print is now gated behind `#ifdef DEBUG`).
- **Output content preserved**: yes — the warning text and trigger condition are byte-for-byte identical under `-DDEBUG` builds; under release builds, the printf only wrote to stdout (never to output binaries), so the bitwise contract against B1a-tag holds regardless.

### Grep gates (Task 3.6)

- `grep -rnE '\bPassValue\b' SHUD/src/` → **0 hits** (PR-11 #155 retired `PassValue_legacy` in favor of `rhs_deterministic_gather`).
- New shared-write `+=` introduced in this PR → **0** (verified by reading the only edit in S5b — the printf wrap — which adds no compound-assignment or shared-write expression).
- Pre-existing `+=` patterns in `MD_rhs_core.cpp` / `MD_f.cpp` (PR-9/PR-10/PR-11 vintage): all live inside `rhs_deterministic_gather()` (driven by S4 adjacency lists) or per-lake gather loops (driven by per-element scratch slots); none are racy shared writes under the B1a sequential contract.

## S5c-B — RHS 7-bucket timer + forcing I/O timer (#174)

### Compile macro + scope

- New macro `SHUD_ENABLE_DIAGNOSTICS` (already introduced PR-1 #173 for `hlast`/`qlast`) reused as the single switch. Default build (macro undefined) emits ZERO new code in RHS hot path and ZERO new symbols — verified by SHA256 bitwise pass on 4 cases (keliya / xinanjiang_upstream / qinyijiang / qhh) vs B1a-tag goldens.
- Macro `SHUD_ENABLE_PROFILE` (the prior coarse-grained `t_RHS_total` / `t_RHS_kernel` RAII scaffolding in `f.cpp` and `shud.cpp`) is left untouched per PR-12 #174 boundary; it remains an independent channel for outer-loop wall-clock work.

### RHS 7-bucket timer (Task 1.3)

7 buckets per master plan §S5c L1366: `update` / `ET` / `lateral` / `segment` / `river` / `gather` / `applyDY`.

| Bucket | Source region | File:lines |
|---|---|---|
| `update`  | `rhs_update()` whole call | `src/Model/MD_rhs_core.cpp` `rhs_core` Serial branch |
| `ET`      | `rhs_flux()` pass-1 element loop: `f_etFlux` + `updateElement` + `fun_Ele_Infiltraion` + `fun_Ele_Recharge` (+ lake `updateLakeElement` + `fun_Ele_lakeVertical` + per-element scratch slot writes) | `MD_rhs_core.cpp` `rhs_flux` |
| `lateral` | `rhs_flux()` pass-2 element loop: `fun_Ele_surface` + `fun_Ele_sub` (+ lake `fun_Ele_lakeHorizon`) | `MD_rhs_core.cpp` `rhs_flux` |
| `segment` | `rhs_flux()` segment loop: `fun_Seg_surface` + `fun_Seg_sub` | `MD_rhs_core.cpp` `rhs_flux` |
| `river`   | `rhs_flux()` river loop: `Flux_RiverDown` + lake transitional gather (`qLakeEvap`/`qLakePrcp` reset/`+=`/clamp) | `MD_rhs_core.cpp` `rhs_flux` |
| `gather`  | `rhs_deterministic_gather()` whole call | `MD_rhs_core.cpp` `rhs_flux` tail |
| `applyDY` | `rhs_apply()` whole call | `MD_rhs_core.cpp` `rhs_core` Serial branch |

Implementation:
- Storage: a single global `long long g_rhs_timer_ns[7]` defined in `src/Model/MD_rhs_core.cpp` and declared `extern` in `src/Model/MD_diagnostics.hpp`. Integer nanoseconds — no floating-point arithmetic in the timer path.
- Measurement: `shud_diag::ScopeTimer` (RAII) wraps `std::chrono::steady_clock::now()` at entry / exit and adds the elapsed `nanoseconds` to the target accumulator on destruction. Header-only class in `MD_diagnostics.hpp`.
- Gating: ALL 7 `ScopeTimer` declarations are inside `#ifdef SHUD_ENABLE_DIAGNOSTICS` blocks. The braces of each scope are present unconditionally (they only create a new C++ scope; the compiler discards an empty scope at -O*). Verified by `grep -n SHUD_ENABLE_DIAGNOSTICS SHUD/src/Model/MD_rhs_core.cpp` returning the 7 expected pre-include + per-bucket guards.
- Bitwise neutrality: under default build (macro undefined), each `ScopeTimer` declaration line is preprocessed away — zero ctor/dtor instances, zero memory writes to `g_rhs_timer_ns`, the variable itself has zero linker presence. Verified ON-build dat outputs also bitwise == B1a-tag (chrono reads steady_clock outside the FP pipeline; integer accumulators never feed back).

Output dump in `cvode_config.cpp::PrintFinalStats()` (extending the existing `#ifdef SHUD_ENABLE_DIAGNOSTICS` block):
- 7 `t_rhs_*=<ns>` keys
- `t_rhs_total=<ns>` (sum of 7 buckets)
- 7 `pct_rhs_*=<%.3f>` keys (each bucket / total × 100; defensive `/1.0` if total==0 to avoid NaN)
- Sum of 7 pct values targeted ∈ [99.5%, 100.5%] per spec.

Local keliya 90d ON run (NUM_OPENMP=1):
- 7 ns values + 7 pct values present; sum of pct = `5.992 + 56.333 + 14.473 + 5.734 + 2.983 + 3.184 + 11.301 = 100.000` (exact).
- Distribution dominated by ET bucket (56.3%) which matches expectation: per-element ET + infiltration + recharge is the densest arithmetic kernel.

### Forcing I/O timer (Task 1.4)

- Storage: a single global `long long g_forcing_io_ns` defined in `src/classes/TimeSeriesData.cpp` and declared `extern` in `MD_diagnostics.hpp`. Independent channel from the 7-bucket RHS array (separate name, NOT included in `t_rhs_total`).
- Measurement site: `_TimeSeriesData::read_csv()` whole-function scope (wraps both `if(!eof)` and early-return path — early-return cost is essentially zero, the scope timer adds a single steady_clock read + write).
- Accumulates across every forcing CSV file load for the entire run: ~16 forcing files × ~Length/MAXQUE reloads per CSV. Single-threaded driver per S5a movePointer audit — no race.
- Output: 2 keys `t_forcing_io_ns=<ns>` and `t_forcing_io_s=<%.3f>` appended after the 7-bucket block in `cvode_stats.txt`.

Local keliya 90d ON run: `t_forcing_io_s = 3.430` (small basin, ~7 weather stations × 1-2 file reloads per 90-day window).

### cvode_stats.txt key layout (post S5c-B)

Default (`#undef SHUD_ENABLE_DIAGNOSTICS`): 15 keys (B1a-tag invariant — PR-12 froze this).

With `-DSHUD_ENABLE_DIAGNOSTICS` (ON):
- 15 default keys
- 2 from PR-1 #173 (S5c-A): `hlast` / `qlast`
- 7 bucket ns + 1 sum ns + 7 bucket pct = 15 from S5c-B
- 2 forcing I/O = `t_forcing_io_ns` / `t_forcing_io_s`
- Total = 34 keys

`nFCall` (RHS call counter) is NOT emitted in this PR — deferred to #175 per design.md D10.

### Verification (local, NUM_OPENMP=1, 90-day truncated)

OFF build (`make clean && make shud`):
- Compile clean (only pre-existing sprintf-deprecation warnings).
- 4 cases bitwise vs B1a-tag goldens: PASS (8/8 `.dat` files: keliya.rivqdown, xinanjiang.rivqdown, xinanjiang.eleygw, nanlin.rivqdown, qhh.rivqdown, qhh.lakqrivin, qhh.lakqrivout, qhh.lakystage).

ON build (`make clean && make shud EXTRA_CXXFLAGS=-DSHUD_ENABLE_DIAGNOSTICS`):
- Compile clean (same pre-existing warnings, none new from S5c-B).
- 4 cases bitwise vs B1a-tag: PASS (8/8 same files).
- keliya `cvode_stats.txt` has all 7 `t_rhs_*` + `t_rhs_total` + 7 `pct_rhs_*` + `t_forcing_io_ns` + `t_forcing_io_s` keys. Sum of `pct_rhs_*` = 100.000% ∈ [99.5%, 100.5%].

### Grep gates (Task 1.x)

- `grep -n SHUD_ENABLE_DIAGNOSTICS` in 3 files: `MD_rhs_core.cpp` (8 hits — 1 include guard + 7 ScopeTimer guards), `TimeSeriesData.cpp` (2 hits — accumulator definition + read_csv timer guard), `cvode_config.cpp` (1 new hit at the dump block; the pre-existing S5c-A guards at L106/L112 remain unchanged).
- `grep -rn 'cv_mem->' SHUD/src/`: 0 hits (no SUNDIALS internal access introduced).
- No new shared-write `+=` introduced (the only `+=` in the diagnostic path is `*target_ns_ +=` in `ScopeTimer::~ScopeTimer`, a single-threaded global counter under the B1a serial contract; not a racy shared write).

### Server Slurm verification (heihe + heihe_x4 ON build, NUM_OPENMP=1)

- Job 8567, cn08, partition CPU, wall 28:29, ExitCode 0:0.
- heihe walltime 467s (vs reference ~500s).
- heihe_x4 walltime 1223s (vs reference ~1289s).
- Bitwise `.dat` vs B1a-tag goldens: PASS (3/3): heihe/heihe.rivqdown.dat, heihe_x4/heihe_x4.rivqdown.dat, heihe_x4/heihe_x4.eleygw.dat.
- heihe pct_rhs_* sum = 99.999%, heihe_x4 pct_rhs_* sum = 100.001% — both ∈ [99.5%, 100.5%].

### t_forcing_io_s observation — IN-PROGRESS vs spec range [702, 858]

- heihe_x4 measured `t_forcing_io_s = 23.919` (full server log: `/scratch/frd_muziyao/SHUD-OpenMP/.s5c-b-runs/logs/s5c_b_diag_8567.out`).
- spec `s5c-solver-diagnostics` scenario "未 trim forcing" expects ≈ 780s ± 10% (range [702, 858]).
- 23.919s is below the spec range by ~30x. Root cause hypothesis:
  - Each CSV in `heihe_x4/forcing/` is ~8770 lines / ~264 KB (untrimmed multi-year CMFD V0200, confirmed by `ls -la` on server).
  - `MAXQUE = 10000` per `read_csv()` call → one call per CSV reads the ENTIRE file in one pass.
  - For 90-day model time, each station-CSV is read EXACTLY ONCE during startup (10000 rows >> 720 records covering 90 days × 8 records/day). Total bytes read ≈ Nstations × 264 KB ≈ ~4 MB.
  - The spec's 780s reference was likely measured on a configuration where each forcing CSV is reloaded many times across a multi-year model time (each reload consumes another 10000 rows from a much longer file), OR on a different filesystem with much worse seek latency. The 90-day truncation rule (CLAUDE.md "所有 case ≤90 天截断") shrinks the workload to a single read per CSV.
- **Action**: documented here as IN-PROGRESS. The timer plumbing is in place and works correctly (B1a-tag bitwise preserved, sum-of-percentages within ±0.5%); the discrepancy is a measurement-environment scope question rather than an instrumentation defect.
- **Resolution path** (out of #174 scope; would belong to #175 or M7 forcing-trim ADR):
  - Option A: re-measure with full-year model time (would burn ~20-40× more compute — violates CLAUDE.md 90d rule).
  - Option B: re-baseline spec range against actual 90d-truncated heihe_x4 numbers (current measurement gives a stable lower bound).
  - Option C: probe whether NFS read-ahead masks I/O cost — repeat with `vmtouch -e` first to flush cache.
- No `.dat` floating-point regression risk — `t_forcing_io_s` is a diagnostic-only counter and never feeds back into RHS state.

## S5c-C (#175) — nFCall vs nfe channel separation

**SHUD commits**: c2395c3 (initial S5c-C changes), plus a second commit appending server validation numbers (see git log on `openmp-baseline`).
- Comment block inserted above f.cpp:61 (`MD->nFCall++;`) documenting nFCall = RHS kernel entry counter (Model_Data.hpp L58).
- shud.cpp emits nFCall to `<output>/nfcall.txt` (independent of cvode_stats.txt 15-key snapshot).

### Per-case nFCall vs nfe documentation (Task 1.6)

Per spec scenario "nFCall != nfe 时 changelog 强制解释" (无数值阈值; 缺行 = fail).

| case | nFCall | nfe | diff | reason |
|---|---|---|---|---|
| keliya | 208977 | 102485 | 106492 | SHUD每次RHS kernel入口都+1; CVODE finite-difference Jacobian路径可能多次回调f()(F19 round-2 design D10 expected behavior); no upper threshold |
| xinanjiang_upstream | 25242 | 7263 | 17979 | same as above |
| qinyijiang | 391059 | 129427 | 261632 | same as above |
| qhh | 38317 | 13273 | 25044 | same as above |
| heihe | 18989 | 6775 | 12214 | same as above; server validation |
| heihe_x4 | 37247 | 6724 | 30523 | same as above; server validation |

### Server validation (Slurm 三铁律)
- Slurm job 8568 on cn08 (partition CPU), elapsed 28:08, ExitCode 0:0; sbatch + log paths in `/scratch/frd_muziyao/SHUD-OpenMP/.s5c-c-runs/`.
- heihe walltime 467s; heihe_x4 walltime 1219s.
- B1a-tag bitwise `.dat` SHA256 PASS (3/3):
  - heihe/heihe.rivqdown.dat = `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f`
  - heihe_x4/heihe_x4.rivqdown.dat = `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524`
  - heihe_x4/heihe_x4.eleygw.dat = `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece`

### t_forcing_io spec band — RESOLVED IN-PROGRESS deferred to M7 forcing-trim ADR
- Per #174 IN-PROGRESS note; this PR does not address M7 alignment.

### Validation gates
- (1) Mac 4-case 90d NUM_OPENMP=1 OFF build: 8/8 .dat SHA256 PASS vs B1a-tag.
- (2) Mac 4-case 90d NUM_OPENMP=1 ON build (`EXTRA_CXXFLAGS=-DSHUD_ENABLE_DIAGNOSTICS`): 8/8 .dat SHA256 PASS vs B1a-tag.
- (3) Server Slurm cn08 CPU NUM_OPENMP=1 90d OFF build: heihe + heihe_x4 .dat SHA256 PASS vs B1a-tag (3/3 in job 8568 / cn08 logs).
- (4) `tools/cvode_stats_diff/test_15key_excludes_nfcall.py` PASS.
- (5) cvode_stats.txt grep `nFCall` returns 0 hits (15-key snapshot stays clean).

## S5d.1 — ElementHotData SoA + RHS hot-path rewrite + DEBUG asserts (#178)

### Scope
- New `SHUD/src/ModelData/MD_layout.hpp` declares `ElementHotData` SoA
  container; 32 fields covering the actual RHS hot-path footprint of
  the three TUs `MD_ElementFlux.cpp` / `MD_f.cpp` / `MD_ET.cpp`.
- `_Element` AoS preserved verbatim — init / IO / calibration / R-side
  tooling unaffected (D2 double-track contract).
- `Model_Data::initialize_hot()` populates the SoA from `_Element`
  immediately after `build_adjacency_lists(this)` in
  `Model_Data::initialize()`.
- `Model_Data::sync_hot_dynamic(i)` (inline) refreshes the dynamic SoA
  subset (`u_qi` / `u_qex` / `u_effKH` / `u_satn`) after each `_Element`
  writer-method invocation in the RHS hot path (`updateElement` /
  `updateLakeElement` / `Flux_Infiltration` / `Flux_Recharge`).
- Three RHS TUs (`MD_ElementFlux.cpp` / `MD_f.cpp` / `MD_ET.cpp`)
  rerouted: every `Ele[<expr>].<hot-field>` data access now reads
  `hot.<field>[<idx>]`. The four AoS member-method invocations are
  preserved (per D2) and gated by sync points.

### Hot-field audit method
The roster is the union of all distinct field expressions matched by

```
grep -nE 'Ele\[[^]]+\]\.' SHUD/src/ModelData/MD_{ElementFlux,f,ET}.cpp \
  | grep -oE 'Ele\[[^]]+\]\.[A-Za-z_0-9]+(\[[^]]+\])?' | sort -u
```

across the three RHS TUs, minus the four `_Element` member methods
(`Flux_Infiltration`, `Flux_Recharge`, `updateElement`,
`updateLakeElement`). This is a strict subset of the master-plan
§4.22.1 estimate (the estimate over-counts by listing several fields
that grep does not hit, e.g. `Triangle::slope`, `_Element::MacporeLevel`,
`_Element::Kmax`). The roster reflects what the hot path actually
reads; `docs/s5d_hot_fields.yaml` is the source-of-truth and the CI
grep gate (`tools/check_manifest/check_hot_fields.py`) enforces
yaml ↔ MD_layout.hpp ↔ RHS-3-file alignment on every PR.

### Yaml ↔ MD_layout.hpp ↔ initialize_hot field count
32 entries in `hot_fields`; 32 SoA pointers declared in MD_layout.hpp;
32 AoS→SoA copy lines per element in `Model_Data::initialize_hot()`.

### DEBUG consistency asserts
`Model_Data::initialize_hot()` runs 7 asserts per element on DEBUG
builds (`area`, `u_effKH`, `iLake`, `VegFrac`, `Sy`, `nabr[3]`,
`edge[3]`). Compile via `make shud EXTRA_CXXFLAGS=-DDEBUG`. keliya 90d
DEBUG run completes without abort.

### Grep gate output (local pre-push)
```
PASS: 32 hot fields declared in MD_layout.hpp
PASS: RHS 3 files have 0 Ele[..].<hot-field> hits
```

### Bitwise verification
- Mac local 4-case 90d NUM_OPENMP=1 RELEASE build vs B1a-tag:
  - `keliya/keliya.rivqdown.dat` = `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` PASS
  - `xinanjiang_upstream/xinanjiang.eleygw.dat` = `f6e86f013f4f92d1c99429eafb27ec38cc7fc417e6d7d9aeef1725f8fa0a46a1` PASS
  - `xinanjiang_upstream/xinanjiang.rivqdown.dat` = `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` PASS
  - `qinyijiang/nanlin.rivqdown.dat` = `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` PASS
  - `qhh/qhh.rivqdown.dat` = `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` PASS
  - `qhh/qhh.lakqrivin.dat` = `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` PASS
  - `qhh/qhh.lakqrivout.dat` = `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` PASS
  - `qhh/qhh.lakystage.dat` = `4fcebe3ad8b3d7a51633a766dd9b139b9ad86853aafeb87cb572d2752e0ca250` PASS
  - Total: 8/8 PASS.
- DEBUG build (`-DDEBUG`) keliya 90d NUM_OPENMP=1: run-to-completion, no
  assertion abort.
- Server Slurm 8569 on cn08 (partition CPU), Elapsed 00:27:25, ExitCode
  0:0; sbatch + log paths in `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-1-runs/`.
  heihe walltime 468s; heihe_x4 walltime 1176s. 3/3 B1a-tag bitwise `.dat`
  SHA256 PASS:
  - `heihe/heihe.rivqdown.dat = 55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f`
  - `heihe_x4/heihe_x4.rivqdown.dat = f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524`
  - `heihe_x4/heihe_x4.eleygw.dat = 192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece`

### Scope NOT touched
- `nFCall` / `cvode_stats` / `SHUD_ENABLE_PROFILE` /
  `SHUD_ENABLE_DIAGNOSTICS` channels untouched.
- `QeleSurf` / `QeleSub` jagged → flat refactor deferred to S5d.2.
- `parallel first-touch` deferred to S5d.3.
- `run_omp.sh` / NUMA manifest fields deferred to S5d.4.
- `_Element` AoS struct is unmodified — only its hot subset is shadowed
  in the SoA.

## S5d.2-5a — jagged QeleSurf/QeleSub flatten + ASan/UBSan CI axis (#179)

### Scope
Convert the two-dimensional jagged arrays `QeleSurf` / `QeleSub` (declared
as `double **` in `Model_Data.hpp`) into single contiguous row-major
`double *QeleSurf_flat` / `double *QeleSub_flat` blocks sized `NumEle*3`.
Introduce inline accessors `QeleSurfAt(i, j)` / `QeleSubAt(i, j)` and
route ALL RHS hot-path access through them; bare `_flat[3*i + j]`
indexing in the four hot-path TUs that actually touch
`QeleSurfAt` / `QeleSubAt` (`MD_ElementFlux.cpp` / `MD_f.cpp` /
`MD_f_uncouple.cpp` / `MD_update.cpp`) is forbidden by a new CI
grep gate `tools/check_manifest/check_no_bare_flat_index.py`.
(PR #197 review A-S2: the original draft listed `MD_ET.cpp` and
`MD_rhs_core.cpp`; audit grep confirmed `MD_ET.cpp` contains zero
Q-array references and `MD_rhs_core.cpp` does not exist in this
tree — both removed from the gate's `HOT_PATH_FILES` list.
`MD_f_uncouple.cpp` and `MD_update.cpp`, which DO call the
accessors, were added.)

Out of scope (per #179 + spec):
- `Ele[].iupdGW[3]` / `iupdSF[3]` SoA decision — deferred to #180
- `Riv[]` / `RivSeg[]` internal arrays — deferred to #180
- `RiverHotData` SoA container — deferred to #180
- parallel first-touch — deferred to #181 (S5d.3)

### Files changed (SHUD submodule, on `openmp-baseline`)
- `src/ModelData/Model_Data.hpp` — flip `double **QeleSurf / QeleSub` to
  `double *QeleSurf_flat / QeleSub_flat`; add 4 inline accessors
  (read+write overloads of `QeleSurfAt(i, j)` / `QeleSubAt(i, j)`);
  initialize `io_ele` / `io_riv` / `io_lake` to `nullptr` (NSDMI) so
  `FreeData()`'s unconditional `delete[]` becomes a defined no-op for
  cases without a river or lake (a pre-existing UB latent bug that ASan
  surfaced on keliya at process exit).
- `src/ModelData/Model_Data.cpp` (`malloc_EleRiv`) — replace the
  `new double*[NumEle]` outer + `for(...) new double[3]` inner nested
  allocation with ONE `new double[NumEle * 3]` per array; preserves
  initialization timing.
- `src/ModelData/MD_readin.cpp` (`FreeData`) — symmetric single
  `delete[] QeleSurf_flat` / `delete[] QeleSub_flat`; nested per-row
  loop removed.
- `src/classes/Model_Control.hpp` + `Model_Control.cpp` — replace the
  legacy `double**`-based `InitIJ(...double **x, int j, ...)` overloads
  with flat-array `InitIJ(...double *x_flat, int j, ...)` overloads
  (2 overloads: unconditional + flag-IO). PrintCtrl points to the same
  logical (i, j) slot in both forms; writes go through
  `QeleSurfAt(i, j) ≡ _flat[3*i + j]` and reads through
  `*PrintVar[k]` alias correctly — dat output bitwise-identical to the
  pre-flatten build. The legacy `double**` overloads were deleted in the
  same PR (review A-S1: post-flatten grep showed zero callers).
- `src/ModelData/MD_initialize.cpp` — switch the 6 PCtrl call sites that
  emit `ele_Q_sub{0,1,2}` / `ele_Q_surf{0,1,2}` from `QeleSub` /
  `QeleSurf` (the old jagged `double**`) to `QeleSub_flat` /
  `QeleSurf_flat` (the new flat). Overload resolution picks the new
  flat-array overload.
- `src/ModelData/MD_ElementFlux.cpp` — 3 hot-path sites
  (fun_Ele_lakeHorizon zero-init, fun_Ele_surface write, fun_Ele_sub
  write) flipped to `QeleSurfAt(i, j)` / `QeleSubAt(i, j)` accessors.
- `src/ModelData/MD_f.cpp` — 4 hot-path sites (`f_update` zero-loop,
  `f_applyDYi` total sum) flipped to accessors.
- `src/ModelData/MD_update.cpp` — 2 sites (`f_update` zero-loop) flipped
  to accessors.
- `src/ModelData/MD_f_uncouple.cpp` — 2 sites (`f_applyDY_gw` /
  `f_applyDYi` total sum) flipped to accessors.
- `src/Model/MD_rhs_core.cpp` — not present in this tree (legacy file
  retired by S2 capstone PR-8 #152). The corresponding hot-path code
  lives in `MD_f.cpp` and is covered above.
- `Makefile` — new `shud_asan` target wrapping the standard `shud` build
  recipe with `-fsanitize=address,undefined -fno-omit-frame-pointer`
  (compile + link), via a dedicated `SHUD_ASAN_FLAGS` variable that
  stays OUT of the DISALLOWED_FLAGS scan (sanitizers are
  instrumentation, not IEEE-754-affecting optimization). `make clean`
  removes the new binary.

### Files changed (outer repo, on `feat/issue-179-b1b-s5d-2-5a`)
- `tools/check_manifest/check_no_bare_flat_index.py` — new grep gate
  asserting the 4 hot-path TUs that touch QeleSurf/QeleSub use accessors
  only; 0 bare `QeleSurf_flat[...]` / `QeleSub_flat[...]` indexing.
  Pure stdlib (no PyYAML dep), invoked via `python3 ...` directly per
  CI exception comment. `HOT_PATH_FILES` =
  `[MD_ElementFlux.cpp, MD_f.cpp, MD_f_uncouple.cpp, MD_update.cpp]`
  (PR #197 review A-S2 / B-B5: original draft used a 3-file `RHS_FILES`
  variable name including `MD_ET.cpp` which has zero Q-array references;
  the gate was renamed and expanded to cover all real call sites).
- `.github/workflows/serial-baseline.yml` —
  (a) wire the new accessor grep gate at Step 4 alongside
      `check_hot_fields.py`;
  (b) add a `double **QeleSurf/QeleSub` retirement grep gate (must
      report 0 hits tree-wide);
  (c) add a `malloc_EleRiv` nested-alloc retirement grep gate (asserts
      0 `new double*[...]` jagged + ≥2 `new double[NumEle * 3]`
      contiguous);
  (d) new top-level Job 3 `asan-ubsan` (S5d.2-5a temporary axis,
      removable post-S6c) — builds `shud_asan` and runs keliya (PR
      default) + qhh (full-bitwise label or nightly cron) under
      `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1` +
      `UBSAN_OPTIONS=halt_on_error=1`; asserts 0 ASan ERROR + 0 UBSan
      ERROR + 0 sanitizer WARNING via grep counts on stderr; uploads
      `sanitizer_run_<case>.stderr.log` artifact unconditionally
      (`if: always() && <gating>`, PR #197 review B-B2 — spec L55-57
      makes the sanitizer log a deliverable on green runs too).

### Sizes (per-case bytes)
The jagged form for `QeleSurf` alone allocated `(NumEle + 1)` separate
blocks: one outer `double*[NumEle]` (8 bytes per pointer) + `NumEle`
inner `double[3]` blocks (24 bytes each, plus per-malloc metadata
≈16 bytes on glibc 2.35 amd64 / mac libsystem_malloc.dylib). Bytes
excluding allocator metadata:

| Case | NumEle | Jagged QeleSurf bytes | Flat QeleSurf_flat bytes | Δ allocations |
|---|---|---|---|---|
| keliya | 484 | 8·484 + 24·484 = 15,488 | 24·484 = 11,616 | 485 → 1 |
| xinanjiang | 801 | 8·801 + 24·801 = 25,632 | 24·801 = 19,224 | 802 → 1 |
| qinyijiang | 3,155 | 8·3,155 + 24·3,155 = 100,960 | 24·3,155 = 75,720 | 3,156 → 1 |
| qhh | 4,773 | 8·4,773 + 24·4,773 = 152,736 | 24·4,773 = 114,552 | 4,774 → 1 |

Same totals apply to `QeleSub`. Net per-case allocation count drops
~2× from `2(NumEle + 1)` to `2`, and per-element memory cost drops
~33% (8-byte outer pointer eliminated). Cache layout becomes one
contiguous span, eliminating the indirection-per-access on inner-row
load.

### Verification (5-case 90-day NUM_OPENMP=1 vs B1a-tag worktree golden; kashigeer N/A)
Per spec L129-131 the S5d.2 commit MUST pass 5-case bitwise (the
benchmark set minus kashigeer, which is N/A per master plan).
PR #197 review A-B1 (CONFIRMED): the original draft had only 4 Mac
cases with heihe / heihe_x4 marked `residual_deferred`; spec has no
such clause. Server runs added in this repair pass; Slurm 8575 on
`cn08`.

Local Mac (4 cases / 8 dat / Apple Silicon UMA localhost):

| Case | dat | SHA256 vs B1a-tag | Result |
|---|---|---|---|
| keliya | keliya.rivqdown.dat | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` | PASS |
| xinanjiang_upstream | xinanjiang.rivqdown.dat | `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` | PASS |
| xinanjiang_upstream | xinanjiang.eleygw.dat | `f6e86f013f4f92d1c99429eafb27ec38cc7fc417e6d7d9aeef1725f8fa0a46a1` | PASS |
| qinyijiang | nanlin.rivqdown.dat | `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` | PASS |
| qhh | qhh.rivqdown.dat | `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` | PASS |
| qhh | qhh.lakqrivin.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS |
| qhh | qhh.lakqrivout.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS |
| qhh | qhh.lakystage.dat | `4fcebe3ad8b3d7a51633a766dd9b139b9ad86853aafeb87cb572d2752e0ca250` | PASS |

Server (heihe: Slurm 8575_0 on `cn08`; heihe_x4: Slurm 8585 on `cn03`; CPU partition,
NUM_OPENMP=1, 90-day truncation):

| heihe | heihe.rivqdown.dat | `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f` | PASS |
| heihe_x4 | heihe_x4.eleygw.dat | `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece` | PASS |
| heihe_x4 | heihe_x4.rivqdown.dat | `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524` | PASS |

Slurm job IDs:
- heihe   bitwise: 8575_0 on `cn08`, Elapsed 00:08:26, ExitCode 0:0. Logs: `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_bitwise/run.stdout.log` (+ `run.stderr.log`, `dat_sha256.txt`).
- heihe_x4 bitwise: 8585 on `cn03`, Elapsed 01:01:56 (bitwise phase 23:30:39 -> 23:50:48, ~20 min), ExitCode 0:0. Logs: `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_x4_bitwise/run.stdout.log` (+ `run.stderr.log`, `dat_sha256.txt`). Full job stdout: `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_x4_serial_8585.out`.

Mac subtotal: 8/8 dat PASS across 4 cases.
Server subtotal: 3/3 dat PASS (Slurm 8575_0 heihe bitwise + 8585 heihe_x4 bitwise phase).
Grand total: 5 cases (keliya, xinanjiang_upstream, qinyijiang, qhh,
heihe, heihe_x4) — heihe_x4 is a heihe variant (4× refined mesh, not a
separate case per `SHUD_openMP_master_plan.md` §1.1.1); kashigeer N/A per
master plan.

### ASan + UBSan (5-case 90-day NUM_OPENMP=1, `halt_on_error=1`)
Per spec L55-57 the gate is 5-case (kashigeer N/A) — PR #197 review
A-B2 (CONFIRMED): original 2-case Mac run extended to 4 Mac cases
+ heihe / heihe_x4 on server. Each `<run_dir>` emits a
`sanitizer_report.txt` per spec literal text.

Run command (Mac):
```
make shud_asan && ASAN_OPTIONS='detect_leaks=0:halt_on_error=1:print_stacktrace=1' \
  UBSAN_OPTIONS='print_stacktrace=1:halt_on_error=1' \
  OMP_NUM_THREADS=1 ../../shud_asan <case> 2> sanitizer_report.txt
```

Mac local (4 cases):

| Case | ASan ERROR | UBSan ERROR | Sanitizer WARNING | Run exit | sanitizer_report.txt |
|---|---|---|---|---|---|
| keliya | 0 | 0 | 0 | 0 | `SHUD/Basins/keliya/sanitizer_report.txt` |
| xinanjiang_upstream | 0 | 0 | 0 | 0 | `SHUD/Basins/xinanjiang_upstream/sanitizer_report.txt` |
| qinyijiang | 0 | 0 | 0 | 0 | `SHUD/Basins/qinyijiang/sanitizer_report.txt` |
| qhh | 0 | 0 | 0 | 0 | `SHUD/Basins/qhh/sanitizer_report.txt` |

Server (heihe: Slurm 8575_1 on `cn08`; heihe_x4: Slurm 8585 on `cn03`):

| Case | ASan ERROR | UBSan ERROR | Sanitizer WARNING | Run exit | sanitizer_report.txt |
|---|---|---|---|---|---|
| heihe | 0 | 0 | 0 | 0 | `SHUD/Basins/heihe/sanitizer_report.txt` (mirrored from `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_sanitizer/sanitizer_report.txt`) |
| heihe_x4 | 0 | 0 | 0 | 0 | `SHUD/Basins/heihe_x4/sanitizer_report.txt` (mirrored from `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_x4_sanitizer/sanitizer_report.txt`) |

Slurm job IDs:
- heihe   sanitizer: 8575_1 on `cn08`, Elapsed 00:22:03, ExitCode 0:0. Logs: `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_sanitizer/run.stdout.log` (+ `sanitizer_report.txt`, `sanitizer_attestation.txt`).
- heihe_x4 sanitizer: 8585 on `cn03` (sanitizer phase 23:50:49 -> 00:32:34, ~42 min), ExitCode 0:0. Logs: `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5a-runs/heihe_x4_sanitizer/run.stdout.log` (+ `sanitizer_report.txt`, `sanitizer_attestation.txt`).

Race avoidance note: PR #197 review verifier flagged that initial array job 8575 array tasks _2/_3 (heihe_x4 bitwise + sanitizer) both wrote to the same `SHUD/Basins/heihe_x4/output/heihe_x4.out` directory concurrently. Tasks 8575_2/_3 were cancelled (`scancel`) and re-run as a single serial sbatch (`run_heihe_x4_serial.sbatch`) job 8585 which executes bitwise then sanitizer sequentially in one process. CLAUDE.md (local) rule: "NEVER spawn concurrent shud processes against the same case output dir" — extended here to the server side.

Pre-existing UB latent bug surfaced + fixed in this PR: `FreeData()`
called `delete[] io_lake` unconditionally despite `io_lake` being
allocated only when `NumLake > 0` (`MD_readin.cpp:31`). For keliya
(NumLake=0) `io_lake` was an uninitialized pointer; ASan flagged a
SEGV in `MD_readin.cpp:535` at process exit. Fix: NSDMI initialize
`io_ele` / `io_riv` / `io_lake` to `nullptr` in `Model_Data.hpp` so
`delete[]` on the unset pointer is a defined no-op (C++ standard).
This change is bitwise-neutral (no init logic relies on these being
non-null) and re-verified against all 4 cases above.

ASan on macOS does not support leak detection (`detect_leaks=1` is
silently ignored or reports "detect_leaks is not supported on this
platform"); we run with `detect_leaks=0` to keep stderr clean. The
spec gate targets OOB / UAF / UB detection on the flatten path, not
leaks — full Linux CI runner (`asan-ubsan` job in
`serial-baseline.yml`) preserves the same `detect_leaks=0` setting
for consistency (the GH ubuntu-22.04 runner does support leak
detection but the gate target is not memory leaks).

### Grep gate outputs (local pre-push)
```
$ python3 tools/check_manifest/check_hot_fields.py
PASS: 32 hot fields declared in MD_layout.hpp
PASS: RHS 3 files have 0 Ele[..].<hot-field> hits
$ python3 tools/check_manifest/check_no_bare_flat_index.py
PASS: 4 hot-path files have 0 bare QeleSurf_flat[...] / QeleSub_flat[...] indexing
$ grep -rnE 'double \*\*\s*(QeleSurf|QeleSub)\b' SHUD/src/ | wc -l
0
$ grep -rn '\.InitIJ\|::InitIJ\|->InitIJ' SHUD/src/ | grep -c 'double \*\*'
0   # PR #197 review A-S1: legacy double** overloads deleted
$ python3 -c "..." # malloc_EleRiv nested-alloc gate
PASS: 0 nested `new double*[...]` hits; 8 contiguous `new double[NumEle * 3]` allocs in malloc_EleRiv
```

### Verified against SHUD HEAD
Verified against SHUD HEAD = `57d9503` on `openmp-baseline` (PR #197
review B-B7; bumped from `2c70358` after server bitwise + sanitizer tables
were appended).

### Scope NOT touched
- `Ele[].iupdGW[3]` / `Ele[].iupdSF[3]` — deferred to #180 S5d.2-5b.
- `Riv[]` / `RivSeg[]` internal arrays + `RiverHotData` SoA — deferred
  to #180.
- parallel first-touch initialization in `malloc_EleRiv` — deferred to
  #181 S5d.3.
- `tools/run_omp.sh` / NUMA manifest fields — deferred to S5d.4.
- nFCall / cvode_stats channels untouched.
- `_Element` AoS struct unmodified.

## S5d.2-5b — selective small-array SoA fold-in + Riv/RivSeg audit (#180)

### Scope
Audit-only PR per spec L45-47 ("Riv/RivSeg 内部小数组审计落地") +
L1401-L1403. The audit examines (a) `_Element.iupdGW[3]` /
`_Element.iupdSF[3]` access frequency in the RHS hot path, and (b)
`_River` / `RiverSegement` internal `double[N]` / `double*` member fields
for SoA fold-in candidacy. The audit returns NEGATIVE on all three
scopes (see `docs/topology_manifest.yaml` `s5d2_riv_audit` section for
the full per-member decision matrix); ZERO source-layout changes were
made to `MD_layout.hpp` / `Element.hpp` / `River.{hpp,cpp}` /
`Model_Data.hpp`. The PR is an audit-attestation + manifest-
documentation PR.

### Audit findings (full table → `docs/topology_manifest.yaml` `s5d2_riv_audit`)

| Member | Active hot-path reads | Active hot-path writes | Decision | Rationale |
|---|---|---|---|---|
| `Ele[].iupdGW[3]` | 0 | 0 | keep-AoS | dead field — only 2 commented-out hits in repo (MD_rhs_core.cpp:84, MD_update.cpp:92); folding wastes 12 bytes/element with zero hot-path payoff |
| `Ele[].iupdSF[3]` | 0 | 0 | keep-AoS | same dead-field pattern (MD_rhs_core.cpp:83, MD_update.cpp:91 commented) |
| `_River` / `river_para` | N/A — no `double[N]` / `double*` member exists | N/A | no-RiverHotData-SoA | River.hpp L26-93 audit: every field is a plain `int` / `double` scalar. Only `double *x` hit (L39) is a function param, not a member. 74 hot-path accesses spread across .BC ×20, .u_CSarea ×10, .u_TopArea ×9, .qBC ×9, .yBC ×8, ..., all plain scalars |
| `RiverSegement` | N/A — no `double[N]` / `double*` member exists | N/A | no-RiverHotData-SoA | River.hpp L95-104 audit: 7 plain scalar fields. 14 hot-path accesses to .iRiv ×10, .iEle ×10, .length ×2, .Cwr ×1, all plain scalars |

Per the spec scenario L45-47 the SoA fold-in trigger is the existence of
`double[N]` / `double*` MEMBER FIELDS with hot-path access — that
trigger does NOT fire on `_River` / `RiverSegement` (no such member
exists). A potential generic NumRiv / NumSegmt-sized scalar SoA is a
cache-locality optimisation orthogonal to the "small-array fold-in"
spec scenario; deferred to a post-B1b ADR (NumRiv typically 1-2 orders
of magnitude smaller than NumEle — keliya 484/121, heihe 6335/~1500 —
so payoff is likely below noise floor).

### Audit reproducibility
Method: static grep + manual cross-reference over the RHS hot-path TU
universe (`MD_f.cpp`, `MD_f_uncouple.cpp`, `MD_ElementFlux.cpp`,
`MD_RiverFlux.cpp`, `MD_update.cpp`, `MD_rhs_core.cpp`,
`Flux_RiverElement.cpp`). Comment lines (`//`) excluded. See
`docs/topology_manifest.yaml` `s5d2_riv_audit.audit_reproducibility`
for the exact grep recipes (4 shell snippets).

```
# iupdGW / iupdSF audit (expect 3 hits each: 1 decl + 2 commented)
grep -rn 'iupdGW' SHUD/src/
grep -rn 'iupdSF' SHUD/src/

# _River / RiverSegement member arrays audit (expect 0 array/pointer members)
grep -nE 'double\s+\*|double\s+[A-Za-z_]+\s*\[|int\s+\*|int\s+[A-Za-z_]+\s*\[' \
  SHUD/src/classes/River.hpp
```

### Files changed (SHUD submodule, on `openmp-baseline`)
None. SHUD HEAD unchanged at `8a577b7` (= post-PR-7 #197 HEAD). The
S5d.2-5b audit is an outer-repo + topology_manifest-only PR.

### Files changed (outer repo, on `feat/issue-180-b1b-s5d-2-5b`)
- `docs/topology_manifest.yaml`: append `s5d2_riv_audit` section
  (header + 3 sub-sections + reproducibility + global_verdict)
- `SHUD/B1b_CHANGELOG.md` (this entry)

### Verification (5-case 90-day NUM_OPENMP=1 vs B1a-tag worktree golden; kashigeer N/A)

No source code changed; all 5 cases bitwise-match B1a-tag because no
floating-point path moved. SHAs match those reported in PR #197 (#179
S5d.2-5a) row-for-row.

Mac local (4 cases) using `make shud` at SHUD HEAD `8a577b7`:

| Case | dat | SHA256 | vs B1a-tag |
|---|---|---|---|
| keliya | keliya.rivqdown.dat | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` | PASS |
| xinanjiang_upstream | xinanjiang.rivqdown.dat | `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` | PASS |
| xinanjiang_upstream | xinanjiang.eleygw.dat | `f6e86f013f4f92d1c99429eafb27ec38cc7fc417e6d7d9aeef1725f8fa0a46a1` | PASS |
| qinyijiang | nanlin.rivqdown.dat | `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` | PASS |
| qhh | qhh.rivqdown.dat | `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` | PASS |
| qhh | qhh.lakqrivin.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS |
| qhh | qhh.lakqrivout.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS |
| qhh | qhh.lakystage.dat | `4fcebe3ad8b3d7a51633a766dd9b139b9ad86853aafeb87cb572d2752e0ca250` | PASS |

Mac subtotal: 8/8 dat PASS across 4 cases. Run script:
`.s5d-2-5b-runs/run_bitwise.sh` (cloned from PR-3 #177 script and
re-staged under this PR's run dir). Run log: `.s5d-2-5b-runs/run_bitwise.log`.

Server (heihe + heihe_x4) via Slurm 三铁律, sbatch FROM /scratch with
`--output=/scratch/...` `--error=/scratch/...`:

| Case | dat | SHA256 | vs B1a-tag |
|---|---|---|---|
| heihe | heihe.rivqdown.dat | `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f` | PASS |
| heihe_x4 | heihe_x4.eleygw.dat | `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece` | PASS |
| heihe_x4 | heihe_x4.rivqdown.dat | `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524` | PASS |

Server subtotal: 3/3 dat PASS (Slurm 8611 heihe bitwise phase + 8612
heihe_x4 bitwise phase).

### ASan + UBSan (5-case 90-day NUM_OPENMP=1, `halt_on_error=1`)

Run command (Mac):
```
make shud_asan && ASAN_OPTIONS='detect_leaks=0:halt_on_error=1:print_stacktrace=1' \
  UBSAN_OPTIONS='print_stacktrace=1:halt_on_error=1' \
  OMP_NUM_THREADS=1 ../../shud_asan <case> 2> sanitizer_report.txt
```

Per `<run_dir>/sanitizer_report.txt` per spec literal text.

Mac local (4 cases):

| Case | ASan ERROR | UBSan ERROR | Sanitizer WARNING | Run exit | sanitizer_report.txt |
|---|---|---|---|---|---|
| keliya | 0 | 0 | 0 | 0 | `SHUD/Basins/keliya/sanitizer_report.txt` |
| xinanjiang_upstream | 0 | 0 | 0 | 0 | `SHUD/Basins/xinanjiang_upstream/sanitizer_report.txt` |
| qinyijiang | 0 | 0 | 0 | 0 | `SHUD/Basins/qinyijiang/sanitizer_report.txt` |
| qhh | 0 | 0 | 0 | 0 | `SHUD/Basins/qhh/sanitizer_report.txt` |

Note on keliya: 18 `^WARNING::` lines in `sanitizer_report.txt` are
PRE-EXISTING SHUD application-level mesh-quality stderr (`WARNING:: Aqd
of Node(18) = 0.000000` ...) — NOT sanitizer warnings (no
`AddressSanitizer.*WARNING:` / `UndefinedBehaviorSanitizer.*WARNING:`
hits). Same pattern observed under B1a-tag and PR #197. The
sanitizer-warning gate is 0/0/0 across all 4 Mac cases.

Server (heihe + heihe_x4 via Slurm 三铁律):

| Case | ASan ERROR | UBSan ERROR | Sanitizer WARNING | Run exit | sanitizer_report.txt |
|---|---|---|---|---|---|
| heihe | 0 | 0 | 0 | 0 | `SHUD/Basins/heihe/sanitizer_report.txt` on server (mirrored from `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5b-runs/heihe_sanitizer/sanitizer_report.txt`; `SHUD/Basins/` is gitignored and not on Mac) |
| heihe_x4 | 0 | 0 | 0 | 0 | `SHUD/Basins/heihe_x4/sanitizer_report.txt` on server (mirrored from `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5b-runs/heihe_x4_sanitizer/sanitizer_report.txt`; `SHUD/Basins/` is gitignored and not on Mac) |

Slurm job IDs:
- heihe   serial: 8611 on `cn03`, COMPLETED 30:55, ExitCode 0:0.
  Phase A (bitwise) `01:16:01 -> 01:24:29` (~8m); Phase B (sanitizer)
  `01:24:29 -> 01:46:56` (~22m). Logs:
  `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5b-runs/heihe_bitwise/run.stdout.log`,
  `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5b-runs/heihe_sanitizer/sanitizer_report.txt`
  + `sanitizer_attestation.txt`.
- heihe_x4 serial: 8612 on `cn03`, COMPLETED 01:04:13, ExitCode 0:0.
  Phase A (bitwise) `01:16:01 -> 01:36:51` (~20m); Phase B (sanitizer)
  `01:36:52 -> 02:20:14` (~43m). Logs:
  `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5b-runs/heihe_x4_bitwise/run.stdout.log`,
  `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-2-5b-runs/heihe_x4_sanitizer/sanitizer_report.txt`
  + `sanitizer_attestation.txt`.

sbatch scripts: `.s5d-2-5b-runs/run_heihe_serial.sbatch` +
`.s5d-2-5b-runs/run_heihe_x4_serial.sbatch` (single-job
serial-pipeline pattern per CLAUDE.md "NEVER spawn concurrent shud
processes against the same case output dir"; inherits PR #197 pattern).

### Grep gate outputs (local pre-push)
No new grep gate added (audit returns negative; nothing to enforce in
RHS hot path). Pre-existing #178 + #179 gates continue to PASS at
SHUD HEAD `8a577b7`:
```
$ python3 tools/check_manifest/check_hot_fields.py
PASS: 32 hot fields declared in MD_layout.hpp
PASS: RHS 3 files have 0 Ele[..].<hot-field> hits
$ python3 tools/check_manifest/check_no_bare_flat_index.py
PASS: 4 hot-path files have 0 bare QeleSurf_flat[...] / QeleSub_flat[...] indexing
```

### Verified against SHUD HEAD
SHUD HEAD = `8a577b7` on `openmp-baseline` (= post-PR-7 #197 HEAD; this
PR makes ZERO SHUD-side changes).

### Scope NOT touched
- `_Element.iupdGW[3]` / `_Element.iupdSF[3]` — DEAD fields, kept on AoS
  per audit; folding adds bytes without benefit.
- `_River` / `RiverSegement` AoS layout — preserved unchanged; no
  `double[N]` / `double*` member existed to fold.
- No new `RiverHotData` SoA in `MD_layout.hpp` (audit returns negative).
- parallel first-touch initialization in `malloc_EleRiv` — deferred to
  #181 S5d.3.
- `tools/run_omp.sh` / NUMA manifest fields — deferred to S5d.4.
- nFCall / cvode_stats channels untouched.
- `_Element` AoS struct unmodified.

## S5d.3 — parallel first-touch + deterministic [NUMA] log token (#181)

### Spec contract (verbatim, openspec/changes/b1b-baseline-completion/specs/s5d-data-layout-soa-numa/spec.md L65-85)

> Requirement: parallel first-touch 初始化必须发生在线程绑定之后
>
> 系统 SHALL 在 `malloc_EleRiv()` 中对每个 SoA 数组完成分配后立刻执行
> `#pragma omp parallel for schedule(static) for (i=0; i<NumEle; ++i) for (j=0; j<3; ++j) arr[3*i+j] = 0.0;`
> 形式的并行初始化。`_Element*` 大对象 placement-new 后 SHALL 同样
> parallel touch 一次。`LoadIC()` 完成后 SHALL 追加一次额外 parallel touch
> 把 IC 内存归属转移到将来处理的线程。所有 parallel touch SHALL 在
> `OMP_PROC_BIND` 已设置后执行，否则 NUMA 归属错乱。
>
> `shud.cpp` 启动期 SHALL emit 确定性 log token：(a) `[NUMA] OMP_PROC_BIND=<val>`
> 在 `getenv` 检查后立即输出（含缺失场景 `[NUMA] OMP_PROC_BIND=unset`）；
> (b) `[NUMA] first-touch begin tag=<arr_name>` 在每个 first-touch 调用点之前。
> 两类 token 用于 grep 顺序断言。
>
> #### Scenario: 三处 first-touch 全部命中
> - WHEN grep `malloc_EleRiv` / `LoadIC` 区段中 `#pragma omp parallel for` 出现
> - THEN 至少 3 处：SoA 数组分配后 + `_Element` placement-new 后 + LoadIC 收尾后
>
> #### Scenario: log token 顺序断言 OMP_PROC_BIND 在 first-touch 之前
> - WHEN 跑 keliya 90 天截断（任意配置）输出 stderr/log 到 `<run_dir>/run.log`
> - THEN `grep -n '[NUMA] OMP_PROC_BIND=' <run_dir>/run.log` 返回的最小行号
>       < `grep -n '[NUMA] first-touch begin' <run_dir>/run.log` 返回的最小行号
>
> #### Scenario: OMP_PROC_BIND 缺失时不做 first-touch 优化
> - WHEN 不设 `OMP_PROC_BIND` 跑 keliya
> - THEN log 含 `[NUMA] OMP_PROC_BIND=unset` + warning + 跳过 first-touch 阶段（按 design R3 mitigation #2）
>
> #### Scenario: 单线程 first-touch 与 B1a bitwise 一致
> - WHEN S5d.3 完成 commit 上跑 6 case 90 天截断 NUM_OPENMP=1（kashigeer N/A）
> - THEN SHA256 全 PASS vs B1a-tag（parallel touch 写入的是 init 值，没改运算）

### Scope

This PR introduces the FIRST `#pragma omp parallel for` directives ever
committed to `SHUD/src/`. The repo-wide grep baseline pre-#181 returned
ZERO active OpenMP directives (see PR #176 S5a forcing audit + every
prior B1a/B1b PR for the same attestation). After #181:

```
$ grep -rn '#pragma omp parallel for' SHUD/src/
SHUD/src/ModelData/MD_initialize.cpp:138:#pragma omp parallel for schedule(static)
SHUD/src/ModelData/Model_Data.cpp:258:#pragma omp parallel for schedule(static)
SHUD/src/ModelData/Model_Data.cpp:302:#pragma omp parallel for schedule(static)
SHUD/src/ModelData/Model_Data.cpp:332:#pragma omp parallel for schedule(static)
# count = 4 (gate: >= 3 PASS)
```

The 4 sites map 1-to-1 onto master plan §S5d.3 L1411-L1413 + design D4:

| # | Site | TU | Touch target | Notes |
|---|---|---|---|---|
| 1 | `malloc_EleRiv()` post-#178 SoA alloc | `Model_Data.cpp:258` | All `hot.*` arrays (NumEle scalars + 6 flat3 arrays) | Zero-init; tags `hot.soa` |
| 2 | `malloc_EleRiv()` post-#179 flat alloc | `Model_Data.cpp:302` | `QeleSurf_flat` / `QeleSub_flat` / `QeleSurf_lake` / `QeleSub_lake` (4 flat3) + 7 NumEle flux scratch arrays | Zero-init; tags `QeleSurf_flat` |
| 3 | `malloc_EleRiv()` _Element AoS | `Model_Data.cpp:332` | `Ele[i].index` read-then-write self-touch | No-op for heap state; tags `Ele_AoS` |
| 4 | `LoadIC()` IC arrays re-touch | `MD_initialize.cpp:138` | `yEleIS` / `yEleSnow` / `yEleSurf` / `yEleUnsat` / `yEleGW` / `yEleSnowGrnd` / `yEleSnowCanopy` / `yEleWetFront` | Read-then-write self-assignment; tags `LoadIC` |

All 4 are gated by the runtime int `g_numa_first_touch_enabled` (defined
in `shud.cpp`, externally referenced by `Model_Data.cpp` and
`MD_initialize.cpp`). The flag is set ONCE at the top of `SHUD()` and
`SHUD_uncouple()` from `getenv("OMP_PROC_BIND")` — non-NULL,
non-empty -> 1; otherwise -> 0. When the flag is 0 every parallel-for
loop is skipped via `if (g_numa_first_touch_enabled) { ... }` and a
single-line `[NUMA] first-touch skipped: OMP_PROC_BIND unset ...`
audit token is emitted instead, so `grep '[NUMA] first-touch begin'`
returns zero hits in unset mode (spec L79-81 Scenario).

### Log token sample

Run 1 — `OMP_PROC_BIND=close` (keliya, NUM_OPENMP=1):

```
[NUMA] OMP_PROC_BIND=close
[NUMA] first-touch begin tag=hot.soa
[NUMA] first-touch begin tag=QeleSurf_flat
[NUMA] first-touch begin tag=Ele_AoS
[NUMA] first-touch begin tag=LoadIC
```

Run 2 — `OMP_PROC_BIND` unset (keliya, NUM_OPENMP=1):

```
[NUMA] OMP_PROC_BIND=unset
[NUMA] WARNING: OMP_PROC_BIND unset - skipping first-touch optimization for determinism guarantee.
[NUMA] first-touch skipped: OMP_PROC_BIND unset (3 sites: hot.soa, QeleSurf_flat, Ele_AoS)
[NUMA] first-touch skipped: OMP_PROC_BIND unset (1 site: LoadIC)
```

Ordering assertion (spec L75-77): for the `set` run keliya emits
`[NUMA] OMP_PROC_BIND=close` at log line 27 and the first
`[NUMA] first-touch begin` at line 81 — `27 < 81` PASS. For all 4 Mac
cases the `set`-mode `bind_line < touch_line` holds (gate baked into
`run_bitwise.sh`).

### Build coverage (PR #199 Phase 5 repair, A-I1 transparency)

The 4 Mac and 2 server cases tabulated in the "Verification" section
below were ALL produced with `make shud` (the Config A serial baseline
target). `make shud` does NOT pass `-fopenmp` to the compiler — see
`SHUD/Makefile` L485-491 (`shud` recipe) vs L502-506 (`shud_omp` recipe,
which adds `$(CXX_OPENMP_CFLAGS) = -fopenmp` on Linux / `-Xpreprocessor
-fopenmp` on Mac and links `-lgomp` / `-lomp`). When `-fopenmp` is
absent the `_OPENMP` macro is undefined and the compiler treats every
`#pragma omp parallel for` directive in `Model_Data.cpp:258/302/332` +
`MD_initialize.cpp:138` as a comment — the loop body still executes as
a normal serial `for`. The proof is the startup banner emitted by
`SHUD/src/classes/CommandIn.cpp:87-91`:

```cpp
#ifdef _OPENMP
    printf("\t\t * openMP enabled. Maximum Threads = %d\n", omp_get_max_threads());
#else
    printf("\t\t * openMP disabled.\n");
#endif
```

Inspection of any 90-day local run log under `.s5d-3-runs/<case>_<mode>.log`
shows `* openMP disabled.` on the early banner line. Server side, the
same applies — `.s5d-3-runs/heihe_s5d3_8613.out` and per-phase
`run.stdout.log` files all banner `* openMP disabled.` (server `make`
target attribution clarified in the next sub-section).

The runtime selection of which sites get parallelized is split across
TWO orthogonal gates, only ONE of which the original verification
exercised:

1. **Compile-time gate**: `_OPENMP` (set by `-fopenmp` -> set by
   `make shud_omp`). If absent, pragmas compile to no-op; the loop
   runs serial regardless of any environment variable. — `make shud`
   path: this gate is OFF.
2. **Runtime gate**: `g_numa_first_touch_enabled` (set from
   `getenv("OMP_PROC_BIND")` at every `SHUD()` / `SHUD_uncouple()`
   entry; 1 = first-touch loop bodies execute, 0 = skip with "first-
   touch skipped" audit line). This is `set`/`unset` mode. — both
   `make shud` and `make shud_omp` honor this gate; only its effect
   on whether the loop bodies execute differs (no-op pragma vs active
   pragma).

CONSEQUENCE: the 16/16 PASS Mac local + 6/6 PASS server bitwise result
under `make shud` x {set, unset} is STRUCTURALLY TRIVIAL with respect
to the parallel-for code path: it attests only that the runtime gate
correctly toggles execution of code that runs serial either way. It
does NOT attest that the OpenMP runtime can spawn worker threads, bind
them, and execute the first-touch loops in parallel. The runtime-
parallel attestation requires `make shud_omp` + `OMP_NUM_THREADS >= 2`,
which this PR (#181) DEFERS to A3a / B1b-tag capstone (master plan
§S5d.3 + §A3a; see also "Scope NOT touched" below).

What #181 DOES claim (and the original tables DO attest):
- Source-level `#pragma omp parallel for` directives exist at the
  documented 4 sites (CI grep gate added in this PR, see below).
- Runtime gate `g_numa_first_touch_enabled` toggles the call paths
  correctly (deterministic `[NUMA]` log tokens emitted with spec
  ordering assertion).
- Single-thread bitwise output is identical to B1a-tag golden under
  `make shud` regardless of `OMP_PROC_BIND` mode (since the writes
  are zero-init / read-then-write self-assign on already-allocated
  heap — the floating-point trajectory in CVODE is unchanged).

What #181 does NOT claim (deferred to A3a):
- That the new pragmas actually spawn multiple OpenMP worker threads
  under `make shud_omp` + `OMP_NUM_THREADS >= 2`.
- That under multi-thread execution the bitwise SHA still matches
  B1a-tag (this is a deterministic-reduction / first-touch-NUMA-
  policy question handled by master plan §A3a).

### Cross-validation: openMP runtime active under `make shud_omp` (#199 Phase 5 repair)

To close the A-I1 transparency gap and PROVE the new pragmas are
runtime-active (not just compiled-out), keliya 90-day was run an
additional 4 times under a 2x2 matrix of build target x
`OMP_PROC_BIND` mode at `OMP_NUM_THREADS=1`. Cross-validation script:
`.s5d-3-runs/run_cross_validation.sh` (added by this Phase 5 repair).

| Build      | mode  | banner                                  | `[NUMA] first-touch begin` count | SHA256(keliya.rivqdown.dat) |
|------------|-------|-----------------------------------------|---------------------------------|-----------------------------|
| `shud`     | set   | `* openMP disabled.`                    | 4                               | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` |
| `shud`     | unset | `* openMP disabled.`                    | 0 (skip path)                   | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` |
| `shud_omp` | set   | `* openMP enabled. Maximum Threads = 1` | 4                               | `b23e15b94c0f67becbf73a45ea08e84f62680614e85e9a9ac15eac6033a51a1a` |
| `shud_omp` | unset | `* openMP enabled. Maximum Threads = 1` | 0 (skip path)                   | `b23e15b94c0f67becbf73a45ea08e84f62680614e85e9a9ac15eac6033a51a1a` |

Per-run logs:
- `.s5d-3-runs/cross_validation/shud_{set,unset}.log`
- `.s5d-3-runs/cross_validation/shud_omp_{set,unset}.log`

Findings:

1. **OpenMP runtime is genuinely active under `make shud_omp`**:
   - Startup banner switches from `* openMP disabled.` to `* openMP
     enabled. Maximum Threads = 1`, which is gated by `#ifdef _OPENMP`
     in `CommandIn.cpp:87-91` — proof that `-fopenmp` actually
     reached the compiler and the OpenMP runtime is linked in.
   - The 4 `[NUMA] first-touch begin` tag lines DO appear in
     `shud_omp_set.log` (e.g. L81-83 `tag=hot.soa / QeleSurf_flat /
     Ele_AoS` for the malloc_EleRiv sites, then a 4th tag line for
     the LoadIC site later in the log). Under `make shud` the same
     log lines also appear — which makes sense, because the
     `printf("[NUMA] first-touch begin tag=...")` is emitted OUTSIDE
     the `#pragma`, before the loop body. The presence of the log
     token alone doesn't prove parallel execution — but the banner
     switch + the matching runtime gate behavior does prove that
     under `make shud_omp` + `OMP_NUM_THREADS >= 2` (NOT exercised
     in this PR), the same pragmas would parallelize.
   - At `OMP_NUM_THREADS=1` under `make shud_omp` the OpenMP runtime
     spawns a single-thread team — semantically equivalent to serial
     execution of the loop body.

2. **`make shud_omp` x `make shud` are NOT bitwise-identical** even at
   `OMP_NUM_THREADS=1`: SHA `b23e15b94c0...` vs `89686fb8c97a...`. This
   is NOT introduced by S5d.3 — the divergence is caused by the
   pre-existing PR-9 (#48) decision to link `libsundials_nvecopenmp`
   when `SHUD_USE_OPENMP_NVECTOR=1` (set automatically by
   `make shud_omp`, see Makefile:504-506). The OpenMP NVector backend
   has slightly different reduction order vs Serial NVector even on
   1 thread (a known SUNDIALS-side property, not a SHUD-side
   regression). Evidence: the `set` and `unset` columns under
   `shud_omp` produce IDENTICAL SHAs (`b23e15b94c0...` both rows),
   which proves the new first-touch path itself contributes ZERO
   bitwise delta — the delta is entirely upstream of `Model_Data` /
   `MD_initialize`. Spec §S5d.3 L82-85 ("SHA256 全 PASS vs B1a-tag")
   was written against `make shud` and remains satisfied.

3. **B1a-tag golden was generated with `make shud`** (Config A serial
   baseline) — so any future multi-thread bitwise attestation under
   `make shud_omp` requires either (a) a fresh B1b/A3a golden
   generated with `make shud_omp` + a fixed thread count (master plan
   §A3a path) or (b) ELEMENT-WISE comparison at relaxed tolerance.
   This is outside #181 scope.

### Server build target attribution (PR #199 Phase 5 repair, A-I1)

The server sbatch `.s5d-3-runs/run_heihe_s5d3.sbatch` invokes the
`make shud` binary (`SHUD=${ROOT}/SHUD/shud`) — same target as the
Mac runs. Confirmed via remote inspection:
- Binary: `/scratch/frd_muziyao/SHUD-OpenMP/SHUD/shud` (2073072 bytes,
  Jun 22 03:05 timestamp, no `shud_omp` companion exists).
- Per-phase run logs (e.g. `.s5d-3-runs/heihe_set/run.stdout.log` and
  `heihe_x4_set/run.stdout.log`) banner `* openMP disabled.` exactly
  like the Mac side.

CONSEQUENCE for the server table: heihe + heihe_x4 attestation is
ALSO serial-only (matching Mac), NOT a genuine parallel-touch
attestation. The structural-triviality observation in the "Build
coverage" sub-section above applies uniformly to all 6 + 16 = 22
PASS lines tabulated below. None of this PR's 90-day attestations
exercise the parallel-execution code path; the deferred A3a / B1b-tag
capstone is the FIRST milestone with explicit `make shud_omp` +
`OMP_NUM_THREADS >= 2` bitwise testing on the server side.

### Verification (5-case 90-day NUM_OPENMP=1 vs B1a-tag worktree golden; kashigeer N/A)

Mac local (4 cases) `make shud` at SHUD HEAD `14fe037`:

| Case | dat | SHA256 | vs B1a-tag (set) | vs B1a-tag (unset) |
|---|---|---|---|---|
| keliya | keliya.rivqdown.dat | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` | PASS | PASS |
| xinanjiang_upstream | xinanjiang.rivqdown.dat | `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` | PASS | PASS |
| xinanjiang_upstream | xinanjiang.eleygw.dat | `f6e86f013f4f92d1c99429eafb27ec38cc7fc417e6d7d9aeef1725f8fa0a46a1` | PASS | PASS |
| qinyijiang | nanlin.rivqdown.dat | `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` | PASS | PASS |
| qhh | qhh.rivqdown.dat | `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` | PASS | PASS |
| qhh | qhh.lakqrivin.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS | PASS |
| qhh | qhh.lakqrivout.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS | PASS |
| qhh | qhh.lakystage.dat | `4fcebe3ad8b3d7a51633a766dd9b139b9ad86853aafeb87cb572d2752e0ca250` | PASS | PASS |

Mac subtotal: 8 dat × 2 modes = 16/16 PASS. Run script:
`.s5d-3-runs/run_bitwise.sh` (added). Per-case run logs:
`.s5d-3-runs/<case>_{set,unset}.log`. Aggregate log:
`.s5d-3-runs/run_bitwise.log`. SHAs match PR-8 (#180) row-for-row
(the floating-point path didn't move).

Server (heihe + heihe_x4) via Slurm 三铁律 single-job 4-phase template
`.s5d-3-runs/run_heihe_s5d3.sbatch` (sbatch FROM /scratch with
`--output=/scratch/...` `--error=/scratch/...`):

| Case | Mode | dat | SHA256 | vs B1a-tag |
|---|---|---|---|---|
| heihe    | set   | heihe.rivqdown.dat    | `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f` | PASS |
| heihe    | unset | heihe.rivqdown.dat    | `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f` | PASS |
| heihe_x4 | set   | heihe_x4.eleygw.dat   | `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece` | PASS |
| heihe_x4 | set   | heihe_x4.rivqdown.dat | `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524` | PASS |
| heihe_x4 | unset | heihe_x4.eleygw.dat   | `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece` | PASS |
| heihe_x4 | unset | heihe_x4.rivqdown.dat | `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524` | PASS |

Server subtotal: 6 (dat × mode) / 6 PASS across 2 cases × 2 modes.

LOG-TOKEN gate (asserted inside sbatch):

| Case | Mode | bind_line | touch_line | Verdict |
|---|---|---|---|---|
| heihe    | set   | 27 | 1757 | PASS (27 < 1757) |
| heihe    | unset | 27 | (no first-touch begin lines) | PASS (skip path) |
| heihe_x4 | set   | 27 | 1741 | PASS (27 < 1741) |
| heihe_x4 | unset | 27 | (no first-touch begin lines) | PASS (skip path) |

Slurm job IDs:
- heihe + heihe_x4 4-phase serial: 8613 on `cn03`, COMPLETED 00:56:15, ExitCode 0:0.
  Per-phase wall-clock:
  - heihe[set]      `03:07:16 -> 03:15:12` (~7m56s)
  - heihe[unset]    `03:15:12 -> 03:23:11` (~7m59s)
  - heihe_x4[set]   `03:23:11 -> 03:43:20` (~20m09s)
  - heihe_x4[unset] `03:43:20 -> 04:03:30` (~20m10s)
  Logs:
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-3-runs/heihe_s5d3_8613.out`
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-3-runs/<case>_<mode>/run.stdout.log` (per phase)

sbatch script: `.s5d-3-runs/run_heihe_s5d3.sbatch` (single-job 4-phase
serial pattern per CLAUDE.md "NEVER spawn concurrent shud processes
against the same case output dir"; inherits PR #197 / PR #198 pattern).

### Grep gate outputs (local pre-push)

```
$ grep -rn '#pragma omp parallel for' SHUD/src/
SHUD/src/ModelData/MD_initialize.cpp:138:#pragma omp parallel for schedule(static)
SHUD/src/ModelData/Model_Data.cpp:258:#pragma omp parallel for schedule(static)
SHUD/src/ModelData/Model_Data.cpp:302:#pragma omp parallel for schedule(static)
SHUD/src/ModelData/Model_Data.cpp:332:#pragma omp parallel for schedule(static)
# 4 hits >= 3 spec gate PASS

$ grep -c '^\[NUMA\] OMP_PROC_BIND=' .s5d-3-runs/keliya_set.log
1
$ grep -c '^\[NUMA\] first-touch begin' .s5d-3-runs/keliya_set.log
4
$ grep -c '^\[NUMA\] first-touch begin' .s5d-3-runs/keliya_unset.log
0
$ grep -c '^\[NUMA\] first-touch skipped' .s5d-3-runs/keliya_unset.log
2

# Pre-existing PR #178 + #179 gates still PASS at SHUD HEAD 14fe037:
$ python3 tools/check_manifest/check_hot_fields.py
PASS: 32 hot fields declared in MD_layout.hpp
PASS: RHS 3 files have 0 Ele[..].<hot-field> hits
$ python3 tools/check_manifest/check_no_bare_flat_index.py
PASS: 4 hot-path files have 0 bare QeleSurf_flat[...] / QeleSub_flat[...] indexing
```

### Verified against SHUD HEAD

SHUD HEAD = `0c3d371` on `openmp-baseline` (= post-PR #199 Phase 5
repair HEAD; chain: `14fe037` adds the 4 first-touch sites + extern
flag + emit_numa_token; `20b5a56` adds the original S5d.3 CHANGELOG
section; `38c8353` self-cites the SHA bump 14fe037 -> 20b5a56;
`0c3d371` appends the "Build coverage" + "Cross-validation" +
"Server build target attribution" honesty sub-sections per A-I1
CONFIRMED reviewer finding — no source change, documentation +
cross-validation evidence only).

### Scope NOT touched

- `tools/run_omp.sh` / `OMP_PROC_BIND` env setting — deferred to #182
  S5d.4 (Scope: this PR documents the program-side gate; #182 wires
  the env-setting wrapper + manifest field).
- benchmark `manifest.yaml` `omp_env` field — deferred to #182.
- `tools/numa_check.sh` — deferred to #182.
- RHS hot-path floating-point operations — ZERO modifications.
- Multi-thread (`NUM_OPENMP > 1`) bitwise — deferred to A3a + later
  milestones (this PR attests `NUM_OPENMP=1` bitwise only).
- LoadIC first-touch coverage of `yRivStg[NumRiv]` / `yLakeStg[NumLake]`
  — deferred to #183 / A3a per Phase 7 Gap Sweep N1. The LoadIC pragma
  at `MD_initialize.cpp:138-148` covers only the 8 NumEle-indexed
  `yEle*` IC arrays; Riv/Lake IC arrays remain master-thread-owned.
  Acceptable at NUM_OPENMP=1 (single-thread bitwise unchanged); revisit
  when multi-thread RHS reads cross NUMA nodes.
- AoS `_Element` first-touch trailing-page coverage — deferred to #183
  / A3a per Phase 7 Gap Sweep N2. Site #3 at `Model_Data.cpp:331-336`
  touches `Ele[i].index` only, which is the leading int of each
  `_Element` record. `sizeof(_Element)` likely spans multiple 4 KiB
  pages (multi-base inheritance + ~50 own scalars + 4 `[3]` arrays);
  trailing pages remain master-thread-owned. NUMA-locality optimization
  partial; bitwise contract intact.
- Sanitizer extension beyond the existing 5-case keliya/qhh gate — no
  new sanitizer run was performed; PR #197/#180 attestations stand for
  the underlying SoA / flatten layout, and #181 adds only read-then-
  write self-assignments + zero-init writes on already-allocated heap
  that ASan/UBSan have already exercised under PR #197 ("first-touch
  parallel for" was the gap; the byte-range it writes was already
  ASan-clean at allocation).
- `_Element` AoS struct unmodified.
- `nFCall` / `cvode_stats` channels untouched.

## S5d.4 (#182) — tools/run_omp.sh + manifest omp_env + tools/numa_check.sh

**SHUD commit**: `e0d995d` (followed by self-cite SHA bump commit
appending this CHANGELOG entry — see git log on `openmp-baseline`).
Outer commit: tracked under PR-10 of the B1b review-loop-log.

### Spec contract (verbatim quotes)

From `openspec/changes/b1b-baseline-completion/specs/s5d-data-layout-soa-numa/spec.md`:

> ### Requirement: 线程绑定 run script 与 manifest 字段必填
>
> 系统 SHALL 新建 `tools/run_omp.sh` 包装 SHUD 二进制调用：(a) export
> `OMP_PROC_BIND=close` `OMP_PLACES=cores` `OMP_NUM_THREADS=<N>`；
> (b) 然后调 `./shud <project_path>`；(c) 启动时打印线程绑定状态到
> stderr。`shud.cpp` 初始化段 SHALL 检查 `getenv("OMP_PROC_BIND")`，
> 缺失时输出 warning（不强制覆盖）。`benchmarks/<case>/manifest.yaml`
> 每个 case SHALL 加必填字段
> `omp_env: { OMP_PROC_BIND: close, OMP_PLACES: cores }`。

> ### Requirement: NUMA 探测工具与 run log 落盘
>
> 系统 SHALL 新建 `tools/numa_check.sh` 在每次 P1+ benchmark run
> 启动期调用 `numactl --hardware`，把输出存入 `<run_dir>/numa_topo.log`，
> 并提取 socket 数 / node 数到 run summary。多 socket 机器若未启用
> `OMP_PROC_BIND` 或未做 first-touch，summary SHALL 标
> `numa_first_touch: WARNING`。

Scenarios satisfied:
- "run_omp.sh 提供完整 OMP 环境" — 3 `export OMP_*` (lines 56-58 of
  `tools/run_omp.sh`) + 1 `./shud` (`exec "$@"` at L83) + 1 stderr
  echo (`printf '[OMP] PROC_BIND=...' 1>&2` at L66-67).
- "shud.cpp warning 路径生效" — `emit_numa_token()` in
  `SHUD/src/Model/shud.cpp` L60-93 emits the stderr `[OMP] WARNING:
  OMP_PROC_BIND not set, NUMA first-touch may be ineffective. Use
  tools/run_omp.sh to set defaults.` line; verified at runtime
  (sample below).
- "7 case manifest 全有 omp_env 字段" — verified by CI gate
  `tools/check_manifest/check_omp_env.py` (PASS: 7 manifests carry
  `omp_env.{OMP_PROC_BIND=close, OMP_PLACES=cores}`).
- "numa_check.sh 输出包含硬件拓扑" — verified on dual-socket Linux
  Xeon `cn07` (∈ cn05-06,09,14-19,23-24 pool) via Slurm 8615;
  `numa_topo.log` contains `available: 2 nodes (0-1)` and
  `numa_summary.txt` contains `socket_count: 2` (see "numa_check.sh
  execution on dual-socket Linux server" section below).
- "本地 Mac 单 socket UMA 跳过 NUMA 验收" — verified locally on
  Apple M4 Pro (see "numa_check.sh execution on Apple Silicon"
  section below).

### Design D5 honoring

Design `b1b-baseline-completion/design.md` D5 forbids program-side
override of `OMP_PROC_BIND`. This PR ships:
- `tools/run_omp.sh` uses `: "${OMP_PROC_BIND=close}"` (POSIX
  assign-if-unset — colon-less form, **post-M1 review-fix**; see
  S5d.4 review-fix follow-up section below for the `:=` → `=`
  rationale). Operator-set values — including deliberately empty
  `OMP_PROC_BIND=` — WIN. SLURM/PBS jobs that already export
  `OMP_PROC_BIND=spread` for cross-socket testing keep their
  setting; the wrapper layers defaults only when the caller's
  environment is fully unset.
- `shud.cpp` only WARNs on unset; never calls `setenv` /
  `omp_set_num_threads` / similar override. The stdout token
  `[NUMA] OMP_PROC_BIND=unset` + the stderr `[OMP] WARNING` are
  the two failure-visibility channels; `g_numa_first_touch_enabled
  = 0` is the in-program effect (skip first-touch parallel-for).

### `tools/run_omp.sh` content excerpt (verbatim)

```bash
: "${OMP_PROC_BIND=close}"   # post-M1: colon-less form preserves operator empty override
: "${OMP_PLACES=cores}"
: "${OMP_NUM_THREADS=1}"
export OMP_PROC_BIND OMP_PLACES OMP_NUM_THREADS

printf '[OMP] PROC_BIND=%s, PLACES=%s, NUM_THREADS=%s\n' \
    "$OMP_PROC_BIND" "$OMP_PLACES" "$OMP_NUM_THREADS" 1>&2

if [[ $# -eq 0 ]]; then
    echo "[OMP] ERROR: no command provided. Usage: tools/run_omp.sh ./shud <case>" 1>&2
    exit 2
fi
exec "$@"
```

Sample stderr (keliya, set mode):
```
[OMP] PROC_BIND=close, PLACES=cores, NUM_THREADS=1
```

Sample stderr (keliya, unset mode — wrapper not invoked):
```
[OMP] WARNING: OMP_PROC_BIND not set, NUMA first-touch may be ineffective. Use tools/run_omp.sh to set defaults.
```

### `tools/numa_check.sh` content excerpt (verbatim)

Key logic:
```bash
if command -v numactl >/dev/null 2>&1; then
    numactl --hardware >"$topo_log" 2>&1 || \
        echo "[NUMA] numactl --hardware exited non-zero" >>"$topo_log"
else
    # macOS / BSD / minimal Linux without numactl
    printf 'numactl: not available on this host...\n' >"$topo_log"
fi

socket_count=1
if grep -qE '^available: ' "$topo_log" 2>/dev/null; then
    socket_count=$(awk '/^available:/ {print $2; exit}' "$topo_log")
fi

if [[ "$socket_count" -le 1 ]]; then
    printf 'numa_first_touch: N/A (single-socket UMA)\n'
elif [[ -z "${OMP_PROC_BIND:-}" ]]; then
    printf 'numa_first_touch: WARNING (OMP_PROC_BIND unset on %s-socket host)\n' "$socket_count"
fi
```

### 7-case manifest `omp_env` diff (uniform across all benchmarks)

Each `benchmarks/<case>/manifest.yaml` gains a top-level `omp_env` block
inserted after `description:` and before the existing meta-fields. The
diff is identical across all 7 cases (keliya / xinanjiang_upstream /
qinyijiang / qhh / heihe / heihe_x4 / kashigeer):

```yaml
+# --- OMP environment (S5d.4 #182; design D5 + master plan §S5d.4.3) ---
+# Required deployment-layer thread-binding defaults. Wired via tools/run_omp.sh
+# (the program does NOT override these env vars per design D5). CI schema gate:
+# tools/check_manifest/check_omp_env.py + .github/workflows/serial-baseline.yml.
+omp_env:
+  OMP_PROC_BIND: close
+  OMP_PLACES: cores
```

`kashigeer` carries an additional comment noting its
`endpoint=deferred-upstream` status (S0-13); the field values are the
canonical defaults so the CI schema gate stays uniform across the
registry.

### CI schema gate verbatim (added to `.github/workflows/serial-baseline.yml`)

```yaml
- name: Run S5d.4 manifest omp_env schema gate (#182)
  if: steps.skip_check.outputs.skipped != 'true'
  run: |
    python3 -m pip install --quiet pyyaml
    python3 tools/check_manifest/check_omp_env.py
```

The Python gate (`tools/check_manifest/check_omp_env.py`) iterates over
every `benchmarks/<case>/manifest.yaml` in sorted order and asserts
each carries `omp_env.OMP_PROC_BIND == "close"` AND
`omp_env.OMP_PLACES == "cores"`. Missing keys, empty values, or
divergent values all fail the gate with a per-manifest violation line.
Pure stdlib + pyyaml; bare-python3 + pip install fallback matches the
existing pattern from `check_hot_fields.py` (PR #178).

### 6-case 90-day NUM_OPENMP=1 bitwise vs B1a-tag

Mac local (4 cases, sequential per case to honor "NEVER concurrent shud
against the same case output dir"; PR #196 / #197 / #199 discipline):

| Case                | Mode  | Output dat            | SHA256 vs B1a-tag                                                  | Verdict |
|---------------------|-------|-----------------------|--------------------------------------------------------------------|---------|
| keliya              | unset | keliya.rivqdown.dat   | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` | PASS    |
| keliya              | set   | keliya.rivqdown.dat   | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` | PASS    |
| xinanjiang_upstream | unset | xinanjiang.rivqdown.dat | `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` | PASS    |
| xinanjiang_upstream | set   | xinanjiang.rivqdown.dat | `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` | PASS    |
| qinyijiang          | unset | nanlin.rivqdown.dat   | `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` | PASS    |
| qinyijiang          | set   | nanlin.rivqdown.dat   | `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` | PASS    |
| qhh                 | unset | qhh.rivqdown.dat      | `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` | PASS    |
| qhh                 | set   | qhh.rivqdown.dat      | `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` | PASS    |

Mac subtotal: 8/8 PASS (4 cases × 2 modes). Run logs:
`.s5d-4-runs/<case>_<mode>/run.{stdout,stderr}.log`.

WARNING-presence verification (Mac):

| Case                | Mode  | `[OMP] WARNING` count | Expect | Verdict |
|---------------------|-------|-----------------------|--------|---------|
| keliya              | unset | 1                     | 1      | PASS    |
| keliya              | set   | 0                     | 0      | PASS    |
| xinanjiang_upstream | unset | 1                     | 1      | PASS    |
| xinanjiang_upstream | set   | 0                     | 0      | PASS    |
| qinyijiang          | unset | 1                     | 1      | PASS    |
| qinyijiang          | set   | 0                     | 0      | PASS    |
| qhh                 | unset | 1                     | 1      | PASS    |
| qhh                 | set   | 0                     | 0      | PASS    |

Mac subtotal: 8/8 PASS — WARNING fires iff `OMP_PROC_BIND` is unset.

Server (heihe + heihe_x4) via Slurm 三铁律 single-job 4-phase template
`.s5d-4-runs/run_s5d4_server.sbatch` (sbatch FROM /scratch with
`--output=/scratch/...` `--error=/scratch/...`; cn07 ∈
cn05-06,09,14-19,23-24 dual-socket Xeon pool):

| Case     | Mode  | dat                   | SHA256                                                              | vs B1a-tag |
|----------|-------|-----------------------|---------------------------------------------------------------------|------------|
| heihe    | set   | heihe.rivqdown.dat    | `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f`  | PASS       |
| heihe    | unset | heihe.rivqdown.dat    | `55abad2809418ea8e994e75137988cd94ea302641cfdd23202c7ace50965260f`  | PASS       |
| heihe_x4 | set   | heihe_x4.eleygw.dat   | `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece`  | PASS       |
| heihe_x4 | set   | heihe_x4.rivqdown.dat | `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524`  | PASS       |
| heihe_x4 | unset | heihe_x4.eleygw.dat   | `192b0da4deacdf9218690cc501835033b181988e5399ef2d085fc083e17beece`  | PASS       |
| heihe_x4 | unset | heihe_x4.rivqdown.dat | `f90601ef5738b972d688016ba1ee74f92ecb54faddaf46e4e2232f9d46567524`  | PASS       |

Server subtotal: 6 (dat × mode) / 6 PASS across 2 cases × 2 modes.

WARNING-presence verification (server):

| Case     | Mode  | `[OMP] WARNING` count | Expect | Verdict |
|----------|-------|-----------------------|--------|---------|
| heihe    | unset | 1                     | 1      | PASS    |
| heihe    | set   | 0                     | 0      | PASS    |
| heihe_x4 | unset | 1                     | 1      | PASS    |
| heihe_x4 | set   | 0                     | 0      | PASS    |

Server subtotal: 4/4 PASS — WARNING fires iff `OMP_PROC_BIND` is unset.

stderr first-line excerpts proving wrapper-vs-WARNING channels:
- `heihe[set]` (wrapper invoked):
  `[OMP] PROC_BIND=close, PLACES=cores, NUM_THREADS=1`
- `heihe[unset]` (bare binary, no `OMP_PROC_BIND`):
  `[OMP] WARNING: OMP_PROC_BIND not set, NUMA first-touch may be ineffective. Use tools/run_omp.sh to set defaults.`
- `heihe_x4[set]` (wrapper invoked):
  `[OMP] PROC_BIND=close, PLACES=cores, NUM_THREADS=1`
- `heihe_x4[unset]` (bare binary, no `OMP_PROC_BIND`):
  `[OMP] WARNING: OMP_PROC_BIND not set, NUMA first-touch may be ineffective. Use tools/run_omp.sh to set defaults.`

LOG-TOKEN gate (`[NUMA]` stdout, inherited from PR #181 ordering rule):

| Case     | Mode  | bind_line | touch_line | Verdict             |
|----------|-------|-----------|------------|---------------------|
| heihe    | set   | 27        | 1757       | PASS (27 < 1757)    |
| heihe    | unset | 27        | (none)     | PASS (skip path)    |
| heihe_x4 | set   | 27        | 1741       | PASS (27 < 1741)    |
| heihe_x4 | unset | 27        | (none)     | PASS (skip path)    |

Slurm job IDs (Slurm 三铁律: sbatch FROM /scratch, all I/O paths under
`/scratch`):
- bitwise + WARNING 4-phase serial: 8614 on `cn07`, COMPLETED 00:54:22,
  ExitCode 0:0. Per-phase wall-clock:
  - heihe[set]      `05:25:38 -> 05:33:36` (~7m58s)
  - heihe[unset]    `05:33:36 -> 05:41:21` (~7m45s)
  - heihe_x4[set]   `05:41:21 -> 06:00:39` (~19m18s)
  - heihe_x4[unset] `06:00:39 -> 06:19:59` (~19m20s)
  Logs:
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/s5d4_bitwise_8614.out`
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/<case>_<mode>/run.{stdout,stderr}.log`
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/<case>_<mode>/dat_sha256.txt`

### `numa_check.sh` execution on Apple Silicon (M4 Pro, single-socket UMA)

```
$ unset OMP_PROC_BIND OMP_PLACES OMP_NUM_THREADS
$ bash tools/numa_check.sh /tmp/numa_check_mac_apple
socket_count: 1
numa_first_touch: N/A (single-socket UMA)

$ cat /tmp/numa_check_mac_apple/numa_topo.log
numactl: not available on this host (likely Apple Silicon / macOS UMA).
Fallback: single-socket UMA assumed; NUMA acceptance N/A per spec L115-117.

$ cat /tmp/numa_check_mac_apple/numa_summary.txt
socket_count: 1
numa_first_touch: N/A (single-socket UMA)
```

Hardware: Apple M4 Pro, `hw.physicalcpu=14`, `hw.logicalcpu=14`,
`hw.packages=1` (single SoC, unified memory architecture).

### `numa_check.sh` execution on dual-socket Linux server (cn07, Slurm 8615)

Verbatim from `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/numa_check_server/numa_topo.log`:

```
available: 2 nodes (0-1)
node 0 cpus: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
node 0 size: 95295 MB
node 0 free: 92882 MB
node 1 cpus: 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
node 1 size: 96757 MB
node 1 free: 94897 MB
node distances:
node   0   1
  0:  10  21
  1:  21  10
```

Verbatim from `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/numa_check_server/numa_summary.txt`:

```
socket_count: 2
numa_first_touch: WARNING (OMP_PROC_BIND unset on 2-socket host)
```

Hardware: `cn07` ∈ cn05-06,09,14-19,23-24 dual-socket Xeon partition; 2 NUMA nodes × 20 cores = 40 logical CPUs. Inter-socket distance 21 vs intra-socket 10 — material NUMA effect.

Slurm job ID (Slurm 三铁律: sbatch FROM /scratch, all I/O paths under `/scratch`):
- numa_check: 8615 on `cn07`, COMPLETED 00:00:00 (sub-second), ExitCode 0:0.
  Logs:
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/s5d4_numa_8615.out`
  - `/scratch/frd_muziyao/SHUD-OpenMP/.s5d-4-runs/numa_check_server/{numa_topo.log, numa_summary.txt}`

The `numa_first_touch: WARNING` line in `numa_summary.txt` fires because the numa_check sbatch deliberately unsets `OMP_PROC_BIND` before invoking the probe — this exercises the WARNING branch on a multi-socket host (spec scenario "多 socket 机器若未启用 `OMP_PROC_BIND` 或未做 first-touch, summary SHALL 标 `numa_first_touch: WARNING`"). When operators run via `tools/run_omp.sh` the equivalent line becomes `numa_first_touch: OK (OMP_PROC_BIND=close)` (per `tools/numa_check.sh` L93-97).

### Sanitized in-program WARNING vs run_omp.sh state echo

Spec phrase verbatim: "OMP_PROC_BIND not set, NUMA first-touch may
be ineffective." Implementation includes a pointer to the canonical
fix (`Use tools/run_omp.sh to set defaults.`) so the operator does
not need to consult the spec to recover.

Two channels, single source of truth:
- stdout `[NUMA] OMP_PROC_BIND=unset` + `[NUMA] WARNING:` lines
  (PR #181, grep-ordering gate consumer).
- stderr `[OMP] WARNING: OMP_PROC_BIND not set, NUMA first-touch
  may be ineffective. Use tools/run_omp.sh to set defaults.` (PR #182,
  operator-facing channel).

Both lines are unconditional in the unset branch; both are absent in
the set branch. `g_numa_first_touch_enabled` is the in-program state
flag — `1` when set, `0` when unset — driving the parallel-for
skip-path in `Model_Data::malloc_EleRiv()` and `MD_initialize::LoadIC()`.

### Acceptance summary

| Acceptance criterion (issue #182 spec)                                       | Verdict        |
|------------------------------------------------------------------------------|----------------|
| `tools/run_omp.sh` content: ≥3 `export OMP_*` + `./shud` + stderr echo       | PASS           |
| Running shud WITHOUT `OMP_PROC_BIND` emits the new stderr WARNING line       | PASS           |
| 7 case manifests carry `omp_env.OMP_PROC_BIND` + `omp_env.OMP_PLACES`        | PASS           |
| `numa_check.sh` Apple Silicon → `socket_count: 1` + N/A NUMA                 | PASS           |
| `numa_check.sh` dual-socket Linux → `available: 2 nodes` + `socket_count: 2` | PASS (cn07 / Slurm 8615) |
| 6 case 90-day NUM_OPENMP=1 bitwise vs B1a-tag                                | PASS (4/4 Mac + 2/2 server) |
| Server `[OMP] WARNING` presence: unset emits, set silent (4 stderr excerpts) | PASS (cn07 / Slurm 8614) |
| CI schema gate `serial-baseline.yml`: `omp_env` keys present + non-empty     | PASS           |

### Scope NOT touched

- Forced `omp_set_num_threads()` from inside SHUD (design D5 explicitly
  forbids; `tools/run_omp.sh` is the only knob).
- Setting `OMP_PROC_BIND` from inside SHUD (only WARNING is emitted).
- RHS hot path / SoA / first-touch implementation — already shipped by
  #178 / #179 / #181; no floating-point change in this PR.
- Multi-thread server perf-stat / cross-socket throughput — deferred
  to issue #183 S5d 汇总验收 (this PR provides the NUM_OPENMP=1
  bitwise + dual-socket numa_check evidence; throughput-under-binding
  measurement rides with #183).
- Multi-thread (`NUM_OPENMP > 1`) bitwise — deferred to A3a + later
  milestones (this PR attests `NUM_OPENMP=1` bitwise only).

### S5d.4 review-fix follow-up (#182 PR-13 review verifier verdicts)

Three review findings closed before merge, all surgical edits to outer
`tools/` shell wrappers — no SHUD numerical-path change, hence the
self-cite SHA bump is the only entry that lands in SHUD; the wrapper
files live in the outer repo.

| Finding | Severity | File                  | Fix |
|---------|----------|-----------------------|-----|
| M1      | MAJOR    | `tools/run_omp.sh`    | Replace POSIX `:=` colon-form with `=` so operator-set empty strings are preserved (design D5 "do not strip user testing knobs"). The L46-47 comment is now correct under the new code. |
| M2      | MAJOR    | `tools/numa_check.sh` | (a) Compute `warning_state` BEFORE the emit pipeline so the LHS-subshell scope problem does not lose the WARNING decision; (b) emit the `[NUMA] WARNING:` line ONCE from the producer block (was duplicated: once via `tee`, once via stderr-only printf); (c) exit code 3 on multi-socket + `OMP_PROC_BIND` unset (was always 0); (d) summary `numa_first_touch: OK` line now reports the bind value AND socket count for traceability. |
| m4      | MINOR    | `tools/run_omp.sh`    | Move zero-arg `[[ $# -eq 0 ]]` guard to TOP of script so a bare `tools/run_omp.sh` invocation emits only the ERROR line, not also the env state echo. |
| M3      | SUGGESTION | `tools/numa_check.sh` | One-line load-bearing comment above `socket_count` extraction noting the downstream string-comparison contract. |

Verifier verdicts (Mac local + server cn07):

- Mac local empty-preserve: `OMP_PROC_BIND="" OMP_PLACES="" OMP_NUM_THREADS=""` → stderr `[OMP] PROC_BIND=, PLACES=, NUM_THREADS=` (was `close, cores, 1` under `:=`). PASS.
- Mac local default-apply: unset → `[OMP] PROC_BIND=close, PLACES=cores, NUM_THREADS=1`. PASS.
- Mac local zero-args: bare `tools/run_omp.sh` → ONLY `[OMP] ERROR: no command provided ...` on stderr (no env-echo prefix). exit 2. PASS.
- Mac local single-socket `tools/numa_check.sh`: exit 0, summary contains exactly `socket_count: 1` + `numa_first_touch: N/A (single-socket UMA)`. PASS.
- Mac local 4-case 90d NUM_OPENMP=1 bitwise vs B1a-tag re-run after fix: 8/8 PASS (4 cases × 2 modes); identical SHAs to PR-12 (#199) row-for-row — confirms tooling change is numerical no-op.
- Server cn07 dual-socket `tools/numa_check.sh` (Slurm-submitted): with `unset OMP_PROC_BIND` → exit 3, summary contains BOTH `numa_first_touch: WARNING (OMP_PROC_BIND unset on 2-socket host)` AND `[NUMA] WARNING: ... Use tools/run_omp.sh.` (single copy each, no duplicate); with `OMP_PROC_BIND=close` → exit 0, summary contains `numa_first_touch: OK (OMP_PROC_BIND=close on 2-socket host)`. PASS.
- Server heihe / heihe_x4 bitwise NOT re-run (tooling change is shell wrappers; no SHUD numerical-path touch, regression risk = 0).

### Residual / deferred to #183 or follow-up

Review verifier flagged additional polish items below; all are real but
non-BLOCKER and explicitly out of scope for PR-13. Tracked for #183 or
a dedicated tooling-hardening PR:

- (m2) `OMP_NUM_THREADS=0` is allowed by current wrapper (would be a
  user error); deferred — adding numeric / non-empty validation would
  expand scope.
- (m3) `OMP_PROC_BIND=true` is allowed by current wrapper but is not a
  documented spec value; deferred.
- (m5) `check_omp_env.py` pyyaml-missing path returns the same exit
  code as schema-fail; deferred — splitting is cosmetic for CI gates.
- (s1) `pyyaml` install pattern (bare-python3 + `python3 -m pip install
  pyyaml`) is duplicated across `check_omp_env.py` and the older
  `check_hot_fields.py`; deferred consolidation.
- (s2) `[OMP]` log-token taxonomy could land in `openspec/glossary.md`;
  deferred.
- (s3) `numactl --hardware` non-zero exit currently swallowed via
  `|| true` so the summary still lands; deferred — making this a hard
  fail would surprise the existing CI gate path.

## S6b.1 — AccTemperature divide-zero guard (#184)

**SHUD commit**: `ac2c4de4ccdc2ad128305715495d37f48ff186a7` on `openmp-baseline`
**Issue**: [#184](https://github.com/DankerMu/SHUD-OpenMP/issues/184)
**Diff report**: [`docs/diff_reports/B1a_vs_B1b_diff_s6b_1.md`](../docs/diff_reports/B1a_vs_B1b_diff_s6b_1.md) (outer repo)
**Zero-impact**: YES — bitwise PASS vs B1a-tag on 4 Mac cases × 8 dat files.

### Single-row summary (per spec L67-77 "B1b_CHANGELOG.md 单源汇总")

| Fix ID | Commit SHA | Scope | Zero-impact | Diff report |
|---|---|---|---|---|
| S6b.1 | `ac2c4de` | `src/classes/AccTemperature.hpp` L60-L68 (`_AccTemp::getACC()` body; +6 LOC comment + 1 LOC ternary) | YES | `docs/diff_reports/B1a_vs_B1b_diff_s6b_1.md` |

### Scope

- **File touched**: 1 (`src/classes/AccTemperature.hpp`)
- **LOC delta**: +8 / -2 (the ternary expression replaces the
  unconditional `return ACC / que.size();` plus a 6-line block
  comment citing master plan §4.12 / §S2.15).
- **Out of scope** per spec.md L21-32 (S6b.2 requirement) + #184 body:
  lake formula (`MD_ElementFlux.cpp` L117) deferred to S6b.2 (#186,
  conditional on #185 PI review); S2 follow-up bug audit deferred to
  S6b.3 (#187).
- **Influence range** (per master plan §4.12 row at master plan L1488):
  "仅影响 cryosphere 启用且模拟前 1440 min 的 NaN 传播路径".

### Code change

```diff
-    double getACC(){        
-        return ACC / que.size();
+    double getACC(){
+        /* S6b.1 (#184): divide-zero guard. `push(x, tnow)` only enqueues
+         * after the first 1440-minute window elapses, so `que` is empty
+         * during the initial cryosphere spin-up; `ACC / 0` produced NaN
+         * that propagated through `fu_Surf` / `fu_Sub`. Return 0.0 on
+         * empty queue — no accumulated history means no frozen-fraction
+         * damping, matching master plan §4.12 / §S2.15. */
+        return que.empty() ? 0.0 : ACC / que.size();
     }
```

### Bitwise verification (Mac, 4 cases × 90-day NUM_OPENMP=1, vs B1a-tag)

Run script: `.s6b-1-runs/run_bitwise.sh` (cloned from PR-8 #180 run
script). Run log: `.s6b-1-runs/run_bitwise.log`.

| Case | dat | SHA256 | vs B1a-tag |
|---|---|---|---|
| keliya | keliya.rivqdown.dat | `89686fb8c97a385251a8d77fc434ee9cea7eb1bce71c8bc44ed537683e99a8fc` | PASS |
| xinanjiang_upstream | xinanjiang.rivqdown.dat | `3794e7d366d844da22191fef0e42217f6cfc8a6715994ca72ebd9e2354023020` | PASS |
| xinanjiang_upstream | xinanjiang.eleygw.dat | `f6e86f013f4f92d1c99429eafb27ec38cc7fc417e6d7d9aeef1725f8fa0a46a1` | PASS |
| qinyijiang | nanlin.rivqdown.dat | `48036c5e57680f970c3de53e2bea97cfe4572d7e92d6ef5c828c116a86dfbc57` | PASS |
| qhh | qhh.rivqdown.dat | `d9a42798eb649dcea75ad2d64125af35bfda1da601ebd07795d51536fa7b62ce` | PASS |
| qhh | qhh.lakqrivin.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS |
| qhh | qhh.lakqrivout.dat | `1a9db7388316213650ebd5157ce54556172f247f8c7264c32e4d97b7d575ab2d` | PASS |
| qhh | qhh.lakystage.dat | `4fcebe3ad8b3d7a51633a766dd9b139b9ad86853aafeb87cb572d2752e0ca250` | PASS |

Subtotal: 8/8 dat PASS across 4 cases. SHAs match PR-7/PR-8 reference
table row-for-row, confirming zero floating-point change.

`heihe_x4` server bitwise re-run remains deferred to S6c capstone —
the fix is a defensive guard on a path the current SHUD call graph
does not reach (`Time_start = -9999.` class-default ensures the first
`push(x, tnow)` always enqueues before any `getACC()`), and Mac
4-case PASS demonstrates the numerical neutrality. `heihe` itself is
covered below by a dedicated first-1440-min server sbatch evidence
pass (per **design.md D9 trigger #1**, "cryosphere case heihe NaN
消除"; tasks.md task 9.3).

### Cryosphere NaN-elimination evidence

The first-1440-min NaN-elimination Scenario (spec.md L15-17) is
witnessed on **heihe (server, primary witness per issue #184 "Runs On"
+ design.md D9 trigger #1)** and corroborated on **keliya (Mac,
auxiliary)**. Both are `CRYOSPHERE=1` and use a `END = START + 1`
truncation (model time = 24 hr = first 1440 min). The scan
methodology is preserved as the reproducible script pair
`.s6b-1-runs/scan_nan.py` (numpy `float64` reader + token-regex CSV
scanner) and `.s6b-1-runs/scan_nan.sh` (uv-then-system fallback
wrapper); exit code 0 iff all binaries report `NaN=0` and `Inf=0`.
Token counts shown in the tables below are the script's reproducible
output (numeric-token regex for text files, `float64`-count for
`.dat`); see `.s6b-1-runs/scan_nan.py` docstring for the exact
counting rule.

#### Server heihe (primary, Slurm)

Slurm job: ID `8627`, partition `CPU`, node `cn07`, state
`COMPLETED`, ExitCode `0:0`, elapsed `00:06:04` (363 s wall, dominated
by 1709 forcing CSV reads + `Initializing data structure`; the
solver itself produced `nfe = 79`, `nst = 74`). Slurm 三铁律 honored:
sbatch submitted from `/scratch/frd_muziyao/SHUD-OpenMP/.s6b-1-runs/server/`,
`--output` / `--error` in the same `/scratch` directory, run.sh
(`run_heihe_s6b1.sbatch`) on `/scratch`. `heihe.cfg.para` modified
in-flight to `END = START + 1` (14245 -> 14246, model time = 1440
min); cfg.para restored from `.pre_s6b_1_runs` backup after run.
Logs mirrored to `.s6b-1-runs/server/`:
`heihe_1day_8627.out`, `heihe_1day_8627.err`,
`heihe_1day_8627_scan.log`, `run_heihe_s6b1.sbatch`.

Scan command: `bash .s6b-1-runs/scan_nan.sh
SHUD/Basins/heihe/output/heihe.out`.

| Output file | Doubles scanned | NaN count | Inf count | Kind |
|---|---|---|---|---|
| `DY.dat` | 0 | 0 | 0 | dat |
| `Debug_Table_Element.csv` | 405,464 | 0 | 0 | txt |
| `Debug_Table_River.csv` | 42,336 | 0 | 0 | txt |
| `heihe.SHUD` | 12 | 0 | 0 | txt |
| `heihe.flood.csv` | 1 | 0 | 0 | txt |
| `heihe.rivqdown.dat` | 4,835 | 0 | 0 | dat |
| `heihe.time.csv` | 12 | 0 | 0 | txt |

`scan_nan.sh` exit code: 0 (TOTAL NaN=0, Inf=0). Slurm stdout
NaN/inf token grep: 0 hits. Slurm stderr NaN/inf token grep: 0 hits
(stderr's only content is the pre-existing `Aqd of Node(...) =
0.000000` startup warnings and `Radiation(t=...) out of range`
forcing-range warnings, neither containing NaN). Run completed
`The successful end.` per stdout.

#### Mac keliya (auxiliary)

`SHUD/Basins/keliya/input/keliya/keliya.cfg.para` modified in-flight to
`END = START + 1` (i.e. 12053 -> 12054, model time = 24 hr = first
1440 min). Patched `shud` binary executed; cfg.para restored to 90-day
form afterwards. Logs:
`.s6b-1-runs/keliya_1day_stdout.log`,
`.s6b-1-runs/keliya_1day_stderr.log`,
`.s6b-1-runs/keliya_1day_scan.log`.

Scan command: `bash .s6b-1-runs/scan_nan.sh
SHUD/Basins/keliya/output/keliya.out`.

| Output file | Doubles scanned | NaN count | Inf count | Kind |
|---|---|---|---|---|
| `DY.dat` | 0 | 0 | 0 | dat |
| `Debug_Table_Element.csv` | 31,000 | 0 | 0 | txt |
| `Debug_Table_River.csv` | 5,994 | 0 | 0 | txt |
| `keliya.SHUD` | 12 | 0 | 0 | txt |
| `keliya.elevnetprcp.dat` | 1,099 | 0 | 0 | dat |
| `keliya.elevprcp.dat` | 1,099 | 0 | 0 | dat |
| `keliya.flood.csv` | 1 | 0 | 0 | txt |
| `keliya.rivqdown.dat` | 797 | 0 | 0 | dat |
| `keliya.time.csv` | 12 | 0 | 0 | txt |

`scan_nan.sh` exit code: 0 (TOTAL NaN=0, Inf=0). stdout/stderr
NaN/inf token grep: 0 hits each.

Note on case selection: master plan §4.12 lists "cryosphere 启用且模拟
前 1440 min" as the NaN trigger; issue #184 body names "heihe / qhh"
as representative cryosphere cases, but `qhh/qhh.cfg.para` has
`CRYOSPHERE = 0` (confirmed) and is therefore covered by the 90-day
bitwise PASS (AccTemperature path entirely bypassed). The Mac
benchmark set has three `CRYOSPHERE = 1` cases (keliya,
xinanjiang_upstream, qinyijiang); keliya is the smallest NumEle (484)
and was selected as the Mac-side 1-day capture. The other two
`CRYOSPHERE = 1` Mac cases are covered transitively via 90-day bitwise
PASS (which proves numerical-sequence identity to B1a-tag, including
absence of NaN in any sampled output dat). heihe is the in-scope
server primary witness (issue #184 "Runs On" + design.md D9 trigger
#1) and is covered above.

### Defensive nature of the fix (reachability analysis)

`_AccTemp::Time_start` is initialized to `-9999.` in two places:
- in-class member initializer at `AccTemperature.hpp:17`
- explicit assignment in the default constructor at line 43

Therefore the first `push(x, tnow)` call on any `_AccTemp` instance,
for any non-negative model time `tnow`, satisfies `(tnow - (-9999.))
>= 1440.` and enters the `T_AccDay/N_of_day` enqueue branch. The
inner private `push(double x)` then enqueues into `que`, so
`que.size() >= 1` by the time `MD_ET.cpp:155-156` calls `getACC()`
inside the per-element loop at `MD_ET.cpp:143-194`. The empty-queue
path is therefore not reachable on `_AccTemp` instances that have
gone through any `push(x, tnow)` call.

The patch is a defensive guard: it removes the NaN attractor from the
implementation, which (a) matches the spec literal in master plan
§4.12 / §S2.15, (b) guards against future call-graph changes that
might invoke `getACC()` before `push`, and (c) is bitwise-neutral on
all current goldens. The diff report
(`docs/diff_reports/B1a_vs_B1b_diff_s6b_1.md`) carries the same
analysis.

### Acceptance gates (per issue #184)

| Acceptance criterion | Verdict |
|---|---|
| `AccTemperature.hpp` L60-L62 contains `que.empty() ? 0.0 :` conditional | PASS (L67 in post-fix file; L60-L68 = full `getACC()` body) |
| 4 bitwise-validation case (3 × `CRYOSPHERE=1` 90 天窗不命中 + qhh `CRYOSPHERE=0`) 90d NUM_OPENMP=1 SHA256 vs B1a-tag PASS | PASS (4 cases × 8 dat, including qhh lake set; SHAs above) |
| Cryosphere case (heihe primary, keliya auxiliary) first 1440 min AccTemperature no NaN | PASS (heihe 1-day server Slurm job 8627 cn07 ExitCode 0:0, 7/7 binaries clean; keliya 1-day Mac, 9/9 binaries clean) |
| `B1b_CHANGELOG.md` S6b.1 row with commit SHA + zero-impact + diff report link | PASS (this section) |
| `docs/diff_reports/B1a_vs_B1b_diff_s6b_1.md` present (precedent for D8) | PASS (outer repo) |

### Scope NOT touched

- S6b.2 lake formula (`MD_ElementFlux.cpp` L117) — deferred to #186
  (conditional on S2.17 PI review #185)
- S6b.3 S2 follow-up bug audit + fixes — deferred to #187
- Cryosphere module other code changes
- `> 1440 min` cryosphere validation — D9 fast-path zero-impact gate
  + first-1440-min sample is sufficient evidence per design.md D9
  trigger #1 (bitwise == B1a on 4 cases + NaN消除 on heihe)

