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
