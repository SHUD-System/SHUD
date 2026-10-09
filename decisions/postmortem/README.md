# Postmortems

An incident write-up for a bug that reached where it should not have — a
published result, a merged PR, a release. The point is not the one-line fix; it
is *why our process let it through*. A postmortem is not a decision record:
records capture deliberate choices and beaten alternatives, postmortems capture
failures in hindsight.

## When to write one

Write a postmortem when all three hold:

1. **Subtle** — the mechanism is non-obvious; a careful engineer would
   re-derive it the hard way.
2. **Systemic** — it escaped because of a gap in tests, tooling, or
   conventions, not a one-off typo.
3. **Costly to rediscover** — real debugging time was spent and would be
   spent again.

Link every guardrail the postmortem produced (tests, rules, gates). Files are
`NNNN-short-slug.md`, numbered sequentially; register each in the table below.

| # | Title |
|---|---|

## Skeleton

```markdown
# Post-mortem NNNN: <one-line incident title>

Status: resolved (<fix reference>)

## Executive summary
One paragraph: what broke, the root cause in plain words, why it escaped, the durable lesson.

## Impact
Which results or users were affected, and what did NOT happen.

## Timeline
Evidence-driven: commits, commands, outputs.

## Root cause
The precise mechanism with file references, and why existing defenses missed it.

## Guardrails added
Each entry names a mechanism that exists: a test file, a gate in tools/ci/, an AGENTS.md rule.
A regression guard must fail when the bug is reintroduced.

## Lessons
```
