## Summary

<!-- One paragraph: what changed and why. -->

## Type of change and evidence

Tick one. The evidence listed for it is required.

- [ ] **Bug fix** — a test that fails before the fix and passes after it.
- [ ] **Conservation / numerical / solver change** — an invariant test (identity, antisymmetry, pure right-hand side, unit invariance, ...).
- [ ] **Physical closure change** — invariant tests still pass; `tests/reference/` updated with a before/after comparison and the physical reason below.
- [ ] **Refactor** — `make regress` passes with `tests/reference/` unchanged.
- [ ] **Build, tooling, documentation**

## Runtime evidence

<!-- Commands actually run and their key output. State all four: commit, build
     variant (shud / shud_omp / options), which SUNDIALS build, configuration.
     Only accepted escape hatch, verbatim: "None — review-only change (reason: ...)". -->

## Results

- [ ] Model results are unchanged (`make regress` passes, reference untouched), **or**
- [ ] Results change: reason and size of the change are given here.

## Risk

<!-- What could break, who is affected, how to roll back. -->

## Checklist

- [ ] Touches the solver loop or right-hand side (`src/Model/`), fluxes and state update (`src/ModelData/`), `src/Equations/`, `src/classes/`, `Makefile`, `configure` or `tests/reference/`: yes / no. If yes, a human reviewed every diff line.
- [ ] Input formats, output file names/layout, command-line options and run-log lines used by scripts are unchanged, or the change is approved and documented.
- [ ] No `_v2`/`_new`/`_old` names, no commented-out code, no scratch directories.
- [ ] Over 400 changed lines only with the `diff-limit-exempt` label and a justification here.

## References

<!-- Issues. -->
