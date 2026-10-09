"""Change-scoped gates: tests accompany source changes, PR size, commit subjects, documentation.

The pure functions take plain lists so test_guardrails.py can exercise them
without a git history.
"""
import os
import re

from common import run

# VersionUpdate.md is a changelog: it legitimately names files that were removed.
DOC_FILES = ("README.md", "OpenMP_Guide.md", "OpenMP_NVector_Determinism.md", "AGENTS.md", "CONTEXT.md")
BUILD_FILES = ("Makefile", "configure")
PATH_ROOTS = ("src/", "tests/", "tools/", "decisions/", "openspec/", ".github/", ".githooks/")


def changed_files(root, base):
    out = run(["git", "diff", "--name-only", "--diff-filter=ACMRD", "%s...HEAD" % base], cwd=root).stdout
    return [f for f in out.splitlines() if f]


def changed_line_count(root, base, excluded):
    """Added plus removed lines against `base`, ignoring whitespace and `excluded` path prefixes."""
    cmd = ["git", "diff", "--numstat", "--ignore-all-space", "%s...HEAD" % base]
    total = 0
    for line in run(cmd, cwd=root).stdout.splitlines():
        added, removed, path = line.split("\t", 2)
        if added == "-" or any(path.startswith(prefix) for prefix in excluded):
            continue
        total += int(added) + int(removed)
    return total


def check_tests_changed(files):
    """A change under src/ must come with a change under tests/."""
    if any(f.startswith("src/") for f in files) and not any(f.startswith("tests/") for f in files):
        return ["src/ changed without any change under tests/. Add the test, invariant check or "
                "reference update that covers it (AGENTS.md, Conventions -> Tests accompany changes)."]
    return []


def check_docs_changed(files):
    """A change to the build files must come with a documentation change."""
    if any(f in BUILD_FILES for f in files) and not any(f.endswith(".md") for f in files):
        return ["Makefile or configure changed without any .md change. Document the build change "
                "(README.md, OpenMP_Guide.md or VersionUpdate.md) in the same PR."]
    return []


def check_pr_diff(lines, limit, exempt):
    if exempt or lines <= limit:
        return []
    return ["PR changes %d lines, limit is %d. Split it, or add the `diff-limit-exempt` label "
            "and justify it in the PR description." % (lines, limit)]


def check_commit_subject(subject, allowed_types):
    """Conventional commit subject: `type(scope): text` or `type: text`."""
    if subject.startswith(("Merge ", "Revert ")):
        return []
    if re.match(r"^(%s)(\([a-z0-9_./-]+\))?!?: \S" % "|".join(allowed_types), subject):
        return []
    return ["commit subject `%s` is not a conventional commit. Use `<type>(<scope>): <text>` with type in: %s."
            % (subject, ", ".join(allowed_types))]


def commit_subjects(root, base):
    out = run(["git", "log", "--no-merges", "--format=%s", "%s..HEAD" % base], cwd=root).stdout
    return [s for s in out.splitlines() if s]


def make_targets(root):
    targets = set()
    with open(os.path.join(root, "Makefile"), encoding="utf-8") as handle:
        for line in handle:
            match = re.match(r"^([A-Za-z0-9_.][A-Za-z0-9_. -]*):(?!=)", line)
            if match:
                targets.update(t for t in match.group(1).split() if not t.startswith("."))
    return targets


def check_docs(root, doc_files):
    """`make <target>` and repository paths quoted in the documentation must exist."""
    targets = make_targets(root)
    failures = []
    for rel in doc_files:
        path = os.path.join(root, rel)
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
        code = re.findall(r"`([^`\n]+)`", text) + re.findall(r"(?m)^(?:\$ |    |\t)?(make [^\n`]+)$", text)
        for snippet in code:
            words = snippet.strip().split()
            if words and words[0] == "make":
                for word in words[1:]:
                    if "=" in word or word.startswith("-"):
                        continue
                    if not re.match(r"^[A-Za-z0-9_.]+$", word):
                        break
                    if word not in targets:
                        failures.append("%s: `%s` names make target `%s`, which does not exist." % (rel, snippet.strip(), word))
            token = snippet.strip().split(":")[0]
            if token.startswith(PATH_ROOTS) and not re.search(r"[*<>{}$ ]", token):
                if not os.path.exists(os.path.join(root, token)):
                    failures.append("%s: path `%s` does not exist." % (rel, token))
    return sorted(set(failures))
