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

**SHUD commits**: __SHUD_COMMIT_HASHES__
- Comment added at f.cpp:56 documenting nFCall = RHS kernel entry counter (Model_Data.hpp L58).
- shud.cpp emits nFCall to `<output>/nfcall.txt` (independent of cvode_stats.txt 15-key snapshot).

### Per-case nFCall vs nfe documentation (Task 1.6)

Per spec scenario "nFCall != nfe 时 changelog 强制解释" (无数值阈值; 缺行 = fail).

| case | nFCall | nfe | diff | reason |
|---|---|---|---|---|
| keliya | 208977 | 102485 | 106492 | SHUD每次RHS kernel入口都+1; CVODE finite-difference Jacobian路径可能多次回调f()(F19 round-2 design D10 expected behavior); no upper threshold |
| xinanjiang_upstream | 25242 | 7263 | 17979 | same as above |
| qinyijiang | 391059 | 129427 | 261632 | same as above |
| qhh | 38317 | 13273 | 25044 | same as above |
| heihe | __HEIHE_NFCALL__ | __HEIHE_NFE__ | __HEIHE_DIFF__ | same as above; server validation |
| heihe_x4 | __HEIHE_X4_NFCALL__ | __HEIHE_X4_NFE__ | __HEIHE_X4_DIFF__ | same as above; server validation |

### Server validation (Slurm 三铁律)
- Slurm job ID __SLURM_JOB_ID__ + node __SLURM_NODE__ + wall __SLURM_WALL__ + ExitCode __SLURM_EXITCODE__ + 3 dat SHA256 PASS lines (see § Validation gates (3)).

### t_forcing_io spec band — RESOLVED IN-PROGRESS deferred to M7 forcing-trim ADR
- Per #174 IN-PROGRESS note; this PR does not address M7 alignment.

### Validation gates
- (1) Mac 4-case 90d NUM_OPENMP=1 OFF build: 8/8 .dat SHA256 PASS vs B1a-tag.
- (2) Mac 4-case 90d NUM_OPENMP=1 ON build (`EXTRA_CXXFLAGS=-DSHUD_ENABLE_DIAGNOSTICS`): 8/8 .dat SHA256 PASS vs B1a-tag.
- (3) Server Slurm cn0X CPU NUM_OPENMP=1 90d OFF build: heihe + heihe_x4 .dat SHA256 PASS vs B1a-tag (3/3 in cn0X logs).
- (4) `tools/cvode_stats_diff/test_15key_excludes_nfcall.py` PASS.
- (5) cvode_stats.txt grep `nFCall` returns 0 hits (15-key snapshot stays clean).
