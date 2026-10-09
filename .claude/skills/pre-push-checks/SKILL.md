---
name: pre-push-checks
description: Use before pushing, opening a PR, or claiming checks pass on a SHUD branch — selects the smallest gates that cover the outgoing diff instead of reflexively running the full suite.
---

# Pre-Push Checks

Guidance, not a script. The hooks are intentionally narrow (file size, names, secrets, commit subject); CI owns exhaustive coverage. Run the narrowest evidence that would fail for the regression this change could cause, once.

Sources of truth: `AGENTS.md` (`## Verification Matrix`, `### Scoped verification`, `### Tests accompany changes`), `constraints.yaml`.

## Inspect the outgoing change

1. `git status --short --branch` — confirm checkout and branch.
2. `git diff --name-only <verified-base-ref>...HEAD` — explicit base, never guessed; re-run after a rebase.

## Select evidence by surface

| Change touches | Run |
|---|---|
| `src/` | `make lint`, `make test`, `make regress` |
| `tools/ci/`, `constraints.yaml` | `make test-guardrails` and the gate that changed |
| `Makefile`, `configure` | `make shud`, `make shud_omp`, `make regress` |
| Documentation | `make docs-check` |
| Any branch about to become a PR | `make pr-gates BASE=<verified-base-ref>` |

`make coverage` rebuilds the model and removes `./shud`; run it when tests or `src/` changed. Full `make ci` only on request, while diagnosing CI, or for a repo-wide change.

## What a green result does not prove

- `make regress` runs `ccw` for 3 days: no lake, no trapezoidal channel, no frozen soil. A change there is review-only; say so in the PR.
- Bit-identity of `shud` and `shud_omp` holds only with the SUNDIALS built by `./configure`. If `SUNDIALS_DIR` points elsewhere, a failure of layer 1 is expected and says nothing about the change.
- A passing reference comparison means "unchanged within tolerance", not "physically right".

## Failures and reporting

A relevant failure stops the push — never push and hope CI differs. If `make regress` fails because results are meant to change, do not run `UPDATE=1` on your own: that needs a physical reason, a before/after comparison and a record in `decisions/`. A suspected platform-only failure needs proof: exact command, exit code, the platform difference. Report only commands actually run, with commit, build variant, SUNDIALS build and configuration; report pending CI as pending.
