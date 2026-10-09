"""Static gates: file size, naming, scratch directories, complexity, duplication, dead code.

Every check returns a list of failure messages (empty = pass) and takes the
directory it inspects as an argument, so test_guardrails.py can run the same
code on a fixture that must fail and on one that must pass.
"""
import csv
import io
import os
import re

from common import SOURCE_EXT, run, tool


def count_lines(path):
    with open(path, "rb") as handle:
        return sum(1 for _ in handle)


def check_size(root, files, limit, baseline):
    """Source files may not exceed `limit` lines; baselined files may not grow."""
    failures, seen = [], set()
    for rel in files:
        if not rel.endswith(SOURCE_EXT):
            continue
        lines = count_lines(os.path.join(root, rel))
        if rel in baseline:
            seen.add(rel)
            ceiling = int(baseline[rel])
            if lines > ceiling:
                failures.append("%s: %d lines, baseline ceiling is %d (limit %d). "
                                "Move code out instead of growing this file." % (rel, lines, ceiling, limit))
            elif lines <= limit:
                failures.append("%s: now %d lines (<= %d). Remove it from tools/ci/baseline/size.txt "
                                "and lower baseline.counts.size_violations." % (rel, lines, limit))
        elif lines > limit:
            failures.append("%s: %d lines exceeds the %d-line limit. Split the file." % (rel, lines, limit))
    for rel in sorted(set(baseline) - seen):
        failures.append("%s: listed in tools/ci/baseline/size.txt but no longer exists. Remove the entry." % rel)
    return failures


def check_naming(files, suffix_regex, scratch_dirs):
    """Reject version-suffix file names and scratchpad directories."""
    failures = []
    pattern = re.compile(suffix_regex, re.IGNORECASE)
    for rel in files:
        parts = rel.split("/")
        stem = os.path.splitext(parts[-1])[0]
        if pattern.search(stem):
            failures.append("%s: forbidden name suffix. Edit the original file; git keeps the history." % rel)
        for part in parts[:-1]:
            if part in scratch_dirs:
                failures.append("%s: scratch directory `%s/` must not be committed." % (rel, part))
                break
            if pattern.search(part):
                failures.append("%s: directory `%s/` has a forbidden name suffix." % (rel, part))
                break
    return failures


def lizard_functions(root, paths):
    """{`file::function`: max cyclomatic complexity} for all functions under `paths`."""
    out = run([tool("lizard"), "-l", "cpp", "--csv"] + paths, cwd=root, check=False).stdout
    result = {}
    for row in csv.reader(io.StringIO(out)):
        if len(row) < 8 or not row[1].isdigit():
            continue
        key = "%s::%s" % (row[6], row[7])
        result[key] = max(result.get(key, 0), int(row[1]))
    return result


def check_complexity(root, paths, limit, baseline):
    """No function above `limit`; baselined functions may not get more complex."""
    functions = lizard_functions(root, paths)
    if not functions:
        return ["lizard reported no functions under %s; the complexity gate did not run." % " ".join(paths)]
    failures = []
    for key, ccn in sorted(functions.items()):
        if key in baseline:
            if ccn > int(baseline[key]):
                failures.append("%s: complexity %d, baseline ceiling is %s (limit %d)." % (key, ccn, baseline[key], limit))
            elif ccn <= limit:
                failures.append("%s: complexity now %d (<= %d). Remove it from tools/ci/baseline/complexity.txt "
                                "and lower baseline.counts.complexity_violations." % (key, ccn, limit))
        elif ccn > limit:
            failures.append("%s: complexity %d exceeds %d. Split the function." % (key, ccn, limit))
    for key in sorted(set(baseline) - set(functions)):
        failures.append("%s: listed in tools/ci/baseline/complexity.txt but no longer exists. Remove the entry." % key)
    return failures


def duplicate_rate(root, paths):
    out = run([tool("lizard"), "-l", "cpp", "-Eduplicate"] + paths, cwd=root, check=False).stdout
    match = re.search(r"Total duplicate rate:\s*([0-9.]+)%", out)
    return float(match.group(1)) if match else None


def check_duplicates(root, paths, limit, baseline_rate):
    """Duplicate rate may not exceed the limit, or the frozen baseline while that is higher."""
    rate = duplicate_rate(root, paths)
    if rate is None:
        return ["lizard printed no duplicate rate; the duplication gate did not run."]
    ceiling = max(limit, baseline_rate)
    if rate > ceiling + 1e-9:
        return ["duplicate rate %.2f%% exceeds %.2f%% (limit %.2f%%, frozen baseline %.2f%%). "
                "Reuse the existing code instead of copying it." % (rate, ceiling, limit, baseline_rate)]
    return []


def unused_functions(root, paths):
    """Set of `file::function` that cppcheck reports as never used."""
    proc = run(["cppcheck", "--enable=unusedFunction", "--quiet", "--language=c++", "--std=c++14",
                "--template={file}\t{id}\t{message}"] + paths, cwd=root, check=False)
    found = set()
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) == 3 and cols[1] == "unusedFunction":
            match = re.search(r"'([^']+)'", cols[2])
            if match:
                found.add("%s::%s" % (cols[0], match.group(1)))
    return found, proc.returncode


def check_deadcode(root, paths, baseline):
    """No unused function beyond the frozen baseline list."""
    found, code = unused_functions(root, paths)
    if code != 0:
        return ["cppcheck exited with %d; the dead-code gate did not run." % code]
    failures = ["%s: function is never used. Call it or delete it." % key for key in sorted(found - set(baseline))]
    return failures
