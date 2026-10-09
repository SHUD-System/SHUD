# Decision records

One file per non-trivial decision in `decisions/`, named `YYYY-MM-DD-topic.md`. Written in the same PR as the change it explains.

A record is required when a PR changes model results (it then also updates `tests/reference/` and `VersionUpdate.md`), or settles a choice between alternatives that a later reader would otherwise re-open.

Rules:
- `## Alternatives considered` is mandatory — alternatives are recorded, never invented after the fact.
- A decision is superseded by a NEW record cross-linked both ways; never edit an old record into its opposite.
- Update stale facts (paths, names, defaults) in place when code changes them; the decision and its rationale stay as written.
- No index file: the directory listing is the index.
- Incident write-ups go in `decisions/postmortem/`; they record failures in hindsight, not choices.
