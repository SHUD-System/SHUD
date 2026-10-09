"""Runtime gates: self-tests, input/output contract, two-layer regression, line coverage."""
import array
import glob
import json
import math
import os
import re
import shutil
import sys
import tempfile

from common import ROOT, constraint, run, tool

REFERENCE_DIR = os.path.join(ROOT, "tests", "reference")
CONTRACT_DIR = os.path.join(ROOT, "tests", "io_contract")
STATE_FILES = ("eleysurf", "eleyunsat", "eleygw", "rivystage")


def updating():
    """Reference files are rewritten only on request, and never in CI."""
    if os.environ.get("UPDATE") != "1":
        return False
    if os.environ.get("CI"):
        sys.exit("UPDATE=1 is refused in CI: reference results are updated locally and reviewed in the PR.")
    return True


def make(*args):
    cmd = ["make"] + list(args)
    if os.environ.get("SUNDIALS_DIR"):
        cmd.append("SUNDIALS_DIR=" + os.environ["SUNDIALS_DIR"])
    return run(cmd)


def run_case(binary, project, days, threads=2):
    """Run `binary` on a private copy of input/<project> for `days`; return the output directory."""
    work = tempfile.mkdtemp(prefix="shud-gate-")
    shutil.copytree(os.path.join(ROOT, "input", project), os.path.join(work, "input", project))
    cfg = os.path.join(work, "input", project, project + ".cfg.para")
    with open(cfg, encoding="utf-8") as handle:
        text = handle.read()
    text = re.sub(r"(?m)^END\s.*$", "END\t%d" % days, text)
    text = re.sub(r"(?m)^NUM_OPENMP\s.*$", "NUM_OPENMP\t%d" % threads, text)
    with open(cfg, "w", encoding="utf-8") as handle:
        handle.write(text)
    env = {"OMP_NUM_THREADS": str(threads), "SHUD_RHS_THREADS": str(threads), "OMP_PROC_BIND": "close"}
    run([os.path.join(ROOT, binary), "-o", "out", project], cwd=work, env=env)
    return os.path.join(work, "out")


def read_dat(path):
    """SHUD binary output -> (number of columns, list of records without the time column)."""
    with open(path, "rb") as handle:
        handle.seek(1024)
        values = array.array("d")
        values.frombytes(handle.read())
    if sys.byteorder != "little":
        values.byteswap()
    if len(values) < 2:  # e.g. DY.dat, created empty unless a debug option is on
        return 0, []
    nvar = int(values[1])
    body = values[2 + nvar:]
    width = nvar + 1
    records = [body[i * width + 1:(i + 1) * width] for i in range(len(body) // width)]
    return nvar, records


def compare_values(name, values, reference, max_abs, rms):
    """Failure messages if `values` differ from `reference` by more than the tolerances."""
    if len(values) != len(reference):
        return ["%s: %d values, reference has %d." % (name, len(values), len(reference))]
    diffs = [abs(a - b) for a, b in zip(values, reference)]
    if any(math.isnan(d) for d in diffs):
        return ["%s: NaN in the result." % name]
    worst = max(diffs) if diffs else 0.0
    root_mean = math.sqrt(sum(d * d for d in diffs) / len(diffs)) if diffs else 0.0
    if worst > max_abs or root_mean > rms:
        return ["%s: differs from the reference (max %.3e, RMS %.3e; allowed %.1e, %.1e). If the change in "
                "results is intended, see AGENTS.md -> Reference results." % (name, worst, root_mean, max_abs, rms)]
    return []


def compare_snapshot(name, lines):
    """Compare `lines` with tests/io_contract/<name>; rewrite it when UPDATE=1."""
    path = os.path.join(CONTRACT_DIR, name)
    if updating():
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
        return []
    with open(path, encoding="utf-8") as handle:
        expected = handle.read().splitlines()
    if lines == expected:
        return []
    changed = sorted(set(lines) ^ set(expected))
    return ["%s changed (%s). Input/output formats are a public contract: restore them, or run "
            "`make test UPDATE=1` and justify the change in the PR." % (name, "; ".join(changed[:6]))]


def io_contract(out_dir, project):
    """Snapshot of the command-line options and of the output files with their column counts."""
    with open(os.path.join(ROOT, "src", "classes", "CommandIn.cpp"), encoding="utf-8", errors="replace") as handle:
        options = re.findall(r'getopt\s*\(\s*argc\s*,\s*argv\s*,\s*"([^"]+)"', handle.read())
    listing = []
    for path in sorted(glob.glob(os.path.join(out_dir, "*"))):
        name = os.path.basename(path)
        listing.append("%s\t%d" % (name, read_dat(path)[0]) if name.endswith(".dat") else name)
    return (compare_snapshot("cli_options.txt", options)
            + compare_snapshot("%s_outputs.txt" % project, listing))


def gate_test():
    make("test_adjacency_fallback")
    make("smoke_configd")
    make("shud")
    project = constraint("regression.project")
    return io_contract(run_case("shud", project, 1), project)


def gate_regress():
    """Layer 1: shud and shud_omp agree bit for bit. Layer 2: states match the committed reference."""
    project, days = constraint("regression.project"), constraint("regression.days")
    make("shud")
    make("shud_omp")
    serial = run_case("shud", project, days)
    parallel = run_case("shud_omp", project, days)
    failures = []
    names = sorted(os.path.basename(p) for p in glob.glob(os.path.join(serial, "*.dat")))
    if not names:
        return ["%s produced no .dat output; the regression gate did not run." % project]
    for name in names:
        other = os.path.join(parallel, name)
        with open(os.path.join(serial, name), "rb") as a:
            if not os.path.exists(other) or a.read() != open(other, "rb").read():
                failures.append("%s: shud and shud_omp outputs differ. They must be bit-identical "
                                "(OpenMP_NVector_Determinism.md)." % name)
    ref_dir = os.path.join(REFERENCE_DIR, project)
    for state in STATE_FILES:
        last = list(read_dat(os.path.join(serial, "%s.%s.dat" % (project, state)))[1][-1])
        ref_path = os.path.join(ref_dir, state + ".ref")
        if updating():
            with open(ref_path, "w", encoding="utf-8") as handle:
                handle.write("".join("%.12e\n" % v for v in last))
            continue
        with open(ref_path, encoding="utf-8") as handle:
            reference = [float(line) for line in handle if line.strip()]
        failures += compare_values(state, last, reference,
                                   constraint("regression.max_abs_m"), constraint("regression.rms_m"))
    return failures


def check_coverage(percent, baseline, target):
    """Coverage may not drop below the frozen baseline while that is under the target."""
    floor = min(baseline, target)
    if percent + 1e-9 < floor:
        return ["line coverage %.1f%% is below the floor %.1f%% (target %d%%). Add tests for the new code."
                % (percent, floor, target)]
    if percent >= baseline + 1.0 and baseline < target:
        return ["line coverage rose to %.1f%% (baseline %.1f%%). Raise baseline.coverage_percent in "
                "constraints.yaml so the gain is kept." % (percent, baseline)]
    return []


def gate_coverage():
    project, days = constraint("regression.project"), constraint("regression.days")
    for stale in glob.glob(os.path.join(ROOT, "*.gcno")) + glob.glob(os.path.join(ROOT, "*.gcda")):
        os.remove(stale)
    binary = os.path.join(ROOT, "shud")
    if os.path.exists(binary):
        os.remove(binary)
    try:
        make("shud", "EXTRA_CXXFLAGS=--coverage")
        # gcda files are written next to the gcno files, wherever the binary runs from.
        run_case("shud", project, days)
        summary = os.path.join(tempfile.mkdtemp(prefix="shud-cov-"), "summary.json")
        run([tool("gcovr"), "-r", ".", "--filter", "src/", "--json-summary", summary])
        with open(summary, encoding="utf-8") as handle:
            percent = float(json.load(handle)["line_percent"])
    finally:
        for leftover in glob.glob(os.path.join(ROOT, "*.gcno")) + glob.glob(os.path.join(ROOT, "*.gcda")):
            os.remove(leftover)
        if os.path.exists(binary):
            os.remove(binary)
    platform = "darwin" if sys.platform == "darwin" else "linux"
    print("line coverage (%s, %s %d days): %.1f%%" % (platform, project, days, percent))
    return check_coverage(percent, constraint("baseline.coverage_percent." + platform),
                          constraint("testing.min_line_coverage.value"))
