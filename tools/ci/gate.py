#!/usr/bin/env python3
"""Entry point of the engineering gates. Called by the Makefile: `make <gate>`.

Thresholds and baselines come from constraints.yaml and tools/ci/baseline/.
AGENTS.md (Enforcement Index) lists what each gate enforces.
"""
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pr_checks  # noqa: E402
import runtime_checks  # noqa: E402
import static_checks  # noqa: E402
from common import ROOT, constraint, read_baseline, run, tool, tracked_files  # noqa: E402

SOURCE_DIRS = ["src"]


def ensure_tools():
    """Create .ci-venv with the pinned Python tools on first use."""
    requirements = os.path.join(ROOT, "tools", "ci", "requirements.txt")
    stamp = os.path.join(ROOT, ".ci-venv", "requirements.stamp")
    with open(requirements, encoding="utf-8") as handle:
        wanted = handle.read()
    if os.path.exists(stamp) and open(stamp, encoding="utf-8").read() == wanted:
        return
    run([sys.executable, "-m", "venv", os.path.join(ROOT, ".ci-venv")])
    run([tool("pip"), "install", "--quiet", "-r", requirements])
    with open(stamp, "w", encoding="utf-8") as handle:
        handle.write(wanted)


def require(program, hint):
    if shutil.which(program) is None:
        sys.exit("%s is not installed (%s). The gate cannot run without it." % (program, hint))


def naming_failures(files):
    return static_checks.check_naming(files, constraint("code_canonicality.forbidden_suffixes.pattern"),
                                      constraint("code_canonicality.scratchpad_directories.paths"))


def size_failures(files):
    return static_checks.check_size(ROOT, files, constraint("size_limits.max_file_lines.value"),
                                    read_baseline("size.txt"))


def baseline_count_failures():
    """constraints.yaml is the ledger: its counts must equal the baseline lists."""
    failures = []
    for key, name in (("size_violations", "size.txt"), ("complexity_violations", "complexity.txt"),
                      ("dead_code", "deadcode.txt")):
        listed, recorded = len(read_baseline(name)), constraint("baseline.counts." + key)
        if listed != recorded:
            failures.append("tools/ci/baseline/%s has %d entries but constraints.yaml baseline.counts.%s is %d. "
                            "Make them agree (counts may only go down)." % (name, listed, key, recorded))
    return failures


def gate_lint():
    ensure_tools()
    require("cppcheck", "apt-get install cppcheck / brew install cppcheck")
    files = tracked_files()
    return (baseline_count_failures() + size_failures(files) + naming_failures(files)
            + static_checks.check_complexity(ROOT, SOURCE_DIRS, constraint("size_limits.max_complexity.value"),
                                             read_baseline("complexity.txt"))
            + static_checks.check_duplicates(ROOT, SOURCE_DIRS,
                                             constraint("anti_drift.duplicate_code_threshold_percent.value"),
                                             constraint("baseline.counts.duplicate_code"))
            + static_checks.check_deadcode(ROOT, SOURCE_DIRS, read_baseline("deadcode.txt")))


def gate_secrets():
    require("gitleaks", "https://github.com/gitleaks/gitleaks")
    proc = run(["gitleaks", "git", "--no-banner", "--redact"], check=False)
    return [] if proc.returncode == 0 else ["gitleaks found secrets:\n" + proc.stdout[-2000:]]


def gate_precommit():
    staged = run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"]).stdout.split("\n")
    staged = [f for f in staged if f and os.path.isfile(os.path.join(ROOT, f))]
    baseline = {k: v for k, v in read_baseline("size.txt").items() if k in staged}
    failures = static_checks.check_size(ROOT, staged, constraint("size_limits.max_file_lines.value"), baseline)
    failures += naming_failures(staged)
    require("gitleaks", "https://github.com/gitleaks/gitleaks")
    proc = run(["gitleaks", "git", "--pre-commit", "--staged", "--no-banner", "--redact"], check=False)
    if proc.returncode != 0:
        failures.append("gitleaks found secrets in the staged changes:\n" + proc.stdout[-2000:])
    return failures


def gate_commit_msg():
    with open(sys.argv[2], encoding="utf-8") as handle:
        subject = handle.readline().strip()
    return pr_checks.check_commit_subject(subject, constraint("commits.conventional_commits_required.allowed_types"))


def gate_docs():
    return pr_checks.check_docs(ROOT, list(pr_checks.DOC_FILES) + ["decisions/README.md"])


def gate_pr():
    """Checks that need the PR base: `make pr-gates BASE=<ref>`."""
    base = os.environ.get("BASE") or "origin/master"
    files = pr_checks.changed_files(ROOT, base)
    failures = pr_checks.check_tests_changed(files) + pr_checks.check_docs_changed(files)
    lines = pr_checks.changed_line_count(ROOT, base, constraint("size_limits.max_pr_diff_lines.excluded_prefixes"))
    print("changed lines against %s: %d" % (base, lines))
    failures += pr_checks.check_pr_diff(lines, constraint("size_limits.max_pr_diff_lines.value"),
                                        os.environ.get("DIFF_LIMIT_EXEMPT") == "1")
    types = constraint("commits.conventional_commits_required.allowed_types")
    subjects = pr_checks.commit_subjects(ROOT, base)
    if os.environ.get("PR_TITLE"):
        subjects.append(os.environ["PR_TITLE"])
    for subject in subjects:
        failures += pr_checks.check_commit_subject(subject, types)
    return failures


def gate_guardrails():
    ensure_tools()
    require("cppcheck", "apt-get install cppcheck / brew install cppcheck")
    import test_guardrails
    return test_guardrails.main()


def gate_install_hooks():
    run(["git", "config", "core.hooksPath", ".githooks"])
    print("git hooks enabled: core.hooksPath = .githooks")
    return []


GATES = {
    "lint": gate_lint,
    "secrets": gate_secrets,
    "docs-check": gate_docs,
    "test-guardrails": gate_guardrails,
    "test": runtime_checks.gate_test,
    "regress": runtime_checks.gate_regress,
    "coverage": lambda: (ensure_tools(), runtime_checks.gate_coverage())[1],
    "pr-gates": gate_pr,
    "precommit": gate_precommit,
    "commit-msg": gate_commit_msg,
    "install-hooks": gate_install_hooks,
}
# `make ci`: everything that does not need a PR base. Coverage is last: it rebuilds ./shud.
CI_ORDER = ["lint", "docs-check", "test-guardrails", "test", "regress", "coverage"]


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else ""
    names = CI_ORDER if name == "ci" else [name]
    if any(n not in GATES for n in names):
        sys.exit("usage: gate.py {ci|%s}" % "|".join(GATES))
    failed = False
    for gate in names:
        print("== %s" % gate, flush=True)
        try:
            failures = GATES[gate]()
        except subprocess.SubprocessError as error:
            failures = [str(error)]
        for failure in failures:
            print("FAIL [%s] %s" % (gate, failure))
        print("%s: %s" % (gate, "FAIL (%d)" % len(failures) if failures else "PASS"), flush=True)
        failed = failed or bool(failures)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
