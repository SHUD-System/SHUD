"""Self-test of the gates: every guard must reject a violation and accept clean input.

A guard that passes on its own violation enforces nothing, so each case below
asserts both directions.
"""
import os
import shutil
import tempfile

import pr_checks
import runtime_checks
import static_checks
from common import constraint

CLEAN = """int add(int a, int b){
    return a + b;
}
int main(){
    return add(1, 2);
}
"""


def write(root, rel, text):
    path = os.path.join(root, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def branchy(name, branches):
    body = "".join("    if(x == %d){ y += %d; }\n" % (i, i) for i in range(branches))
    return "int %s(int x){\n    int y = 0;\n%s    return y;\n}\n" % (name, body)


def repeated_block(name):
    lines = "".join("    total += values[%d] * weights[%d] + offset[%d];\n" % (i, i, i) for i in range(40))
    return "double %s(double *values, double *weights, double *offset){\n    double total = 0.;\n%s    return total;\n}\n" % (name, lines)


def cases(tmp):
    limit = constraint("size_limits.max_file_lines.value")
    ccn = constraint("size_limits.max_complexity.value")
    suffix = constraint("code_canonicality.forbidden_suffixes.pattern")
    scratch = constraint("code_canonicality.scratchpad_directories.paths")
    types = constraint("commits.conventional_commits_required.allowed_types")
    pr_limit = constraint("size_limits.max_pr_diff_lines.value")

    write(tmp, "size/long.cpp", "int x;\n" * (limit + 1))
    write(tmp, "size/short.cpp", "int x;\n" * limit)
    yield ("file size", static_checks.check_size(tmp, ["size/long.cpp"], limit, {}),
           static_checks.check_size(tmp, ["size/short.cpp"], limit, {}))
    yield ("file size, baselined file grows",
           static_checks.check_size(tmp, ["size/long.cpp"], limit, {"size/long.cpp": str(limit)}),
           static_checks.check_size(tmp, ["size/long.cpp"], limit, {"size/long.cpp": str(limit + 1)}))

    yield ("name suffix", static_checks.check_naming(["src/flux_v2.cpp"], suffix, scratch),
           static_checks.check_naming(["src/flux.cpp", "src/MD_update.cpp"], suffix, scratch))
    yield ("scratch directory", static_checks.check_naming(["scratch/notes.cpp"], suffix, scratch),
           static_checks.check_naming(["src/notes.cpp"], suffix, scratch))

    write(tmp, "ccn_bad/a.cpp", branchy("tangled", ccn + 1))
    write(tmp, "ccn_ok/a.cpp", branchy("simple", ccn - 2))
    yield ("complexity", static_checks.check_complexity(tmp, ["ccn_bad"], ccn, {}),
           static_checks.check_complexity(tmp, ["ccn_ok"], ccn, {}))

    write(tmp, "dup_bad/a.cpp", repeated_block("first") + repeated_block("second"))
    write(tmp, "dup_ok/a.cpp", repeated_block("only") + branchy("other", 3))
    yield ("duplicate code", static_checks.check_duplicates(tmp, ["dup_bad"], 3.0, 0.0),
           static_checks.check_duplicates(tmp, ["dup_ok"], 3.0, 0.0))

    write(tmp, "dead_bad/a.cpp", CLEAN + "int orphan(int a){\n    return a;\n}\n")
    write(tmp, "dead_ok/a.cpp", CLEAN)
    yield ("dead code", static_checks.check_deadcode(tmp, ["dead_bad"], {}),
           static_checks.check_deadcode(tmp, ["dead_ok"], {}))

    yield ("tests accompany src changes", pr_checks.check_tests_changed(["src/Model/f.cpp"]),
           pr_checks.check_tests_changed(["src/Model/f.cpp", "tests/test_f.cpp"]))
    yield ("docs accompany build changes", pr_checks.check_docs_changed(["Makefile"]),
           pr_checks.check_docs_changed(["Makefile", "README.md"]))
    yield ("PR diff size", pr_checks.check_pr_diff(pr_limit + 1, pr_limit, False),
           pr_checks.check_pr_diff(pr_limit, pr_limit, False))
    yield ("commit subject", pr_checks.check_commit_subject("fixed stuff", types),
           pr_checks.check_commit_subject("fix(lake): correct bathymetry interpolation", types))

    write(tmp, "docs_bad/Makefile", "build:\n\ttrue\n")
    write(tmp, "docs_bad/README.md", "Run `make deploy` and read `src/missing.cpp`.\n")
    write(tmp, "docs_ok/Makefile", "build:\n\ttrue\n")
    write(tmp, "docs_ok/src/main.cpp", CLEAN)
    write(tmp, "docs_ok/README.md", "Run `make build OPT=1` and read `src/main.cpp:3`.\n")
    yield ("documentation references", pr_checks.check_docs(os.path.join(tmp, "docs_bad"), ["README.md"]),
           pr_checks.check_docs(os.path.join(tmp, "docs_ok"), ["README.md"]))

    reference = [1.0, 2.0, 3.0]
    yield ("regression tolerance", runtime_checks.compare_values("state", [1.0, 2.0, 3.1], reference, 1e-3, 1e-4),
           runtime_checks.compare_values("state", [1.0, 2.0, 3.0 + 1e-6], reference, 1e-3, 1e-4))
    yield ("regression length", runtime_checks.compare_values("state", [1.0, 2.0], reference, 1e-3, 1e-4),
           runtime_checks.compare_values("state", reference, reference, 1e-3, 1e-4))
    yield ("coverage floor", runtime_checks.check_coverage(49.0, 54.0, 80),
           runtime_checks.check_coverage(54.2, 54.0, 80))


def main():
    tmp = tempfile.mkdtemp(prefix="shud-guard-")
    failures = []
    try:
        for name, on_violation, on_clean in cases(tmp):
            rejected, accepted = bool(on_violation), not on_clean
            print("  %-34s rejects violation: %-4s accepts clean: %s"
                  % (name, "yes" if rejected else "NO", "yes" if accepted else "NO"))
            if not rejected:
                failures.append("guard `%s` did not reject its violation." % name)
            if not accepted:
                failures.append("guard `%s` rejected clean input: %s" % (name, on_clean[0]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return failures
