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
# The project targets macOS and Linux. Two years of integration amplify the rounding differences between
# them past any tolerance that still detects a change (PR #40), so each keeps its own reference and floor.
PLATFORM = "darwin" if sys.platform == "darwin" else "linux"
CONTRACT_DIR = os.path.join(ROOT, "tests", "io_contract")
STATE_FILES = ("eleysurf", "eleyunsat", "eleygw", "rivystage")
# diag_consistency.txt (src/Model/MD_consistency_diag.hpp): the volumes and counts compared with the reference,
# the counts that must be 0 (#16), and the interfaces whose residual must stay at rounding level of the gross volume.
DIAG_FILE = "diag_consistency.txt"
LEDGER_KEYS = ("prcp_m3", "netprcp_m3", "canopy_evap_m3", "d_canopy_m3", "d_snow_m3", "residual_m3")
COUNT_KEYS = ("land_elements", "landland_edges", "lakes", "cvode_overshoot_steps")
PURITY_KEYS = {
    "rhs_repeat_differing": "the right-hand side gave a different DY when called twice at the same (t, y)",
    "rhs_aba_differing": "the right-hand side gave a different DY at (t, yA) after a call at another state yB",
    "et_inputs_changed_steps": "the frozen step inputs qEleE_IC / qEleNetPrep of a land element were rewritten "
                               "during CVode",
}
INTERFACES = ("landland_sub", "landland_surf", "landlake_sub", "landlake_surf")


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


def run_case(binary, project, days, threads=2, extra_env=None):
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
    env.update(extra_env or {})
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
                "results is intended, run `make regress UPDATE=1` on macOS and on Linux and give the physical "
                "reason in the PR."
                % (name, worst, root_mean, max_abs, rms)]
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
    ref_dir = os.path.join(REFERENCE_DIR, project, PLATFORM)
    if updating():
        os.makedirs(ref_dir, exist_ok=True)
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
    return failures + check_diagnostics(serial, project, days)


def read_keys(path):
    """`key=value` lines -> {key: number}; a missing file gives {}."""
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as handle:
        return {key: float(value) for key, _, value in (line.strip().partition("=") for line in handle) if key}


def compare_ledger(values, reference, rel):
    """Failure messages if a canopy or snow volume differs from the reference by more than rel * precipitation."""
    allowed = rel * reference["prcp_m3"]
    return ["%s: %.12e m3, reference %.12e m3 (allowed difference %.1e m3). If the change in results is "
            "intended, run `make regress UPDATE=1` on macOS and on Linux and give the physical reason in the PR."
            % (key, values[key], reference[key], allowed)
            for key in LEDGER_KEYS if not abs(values[key] - reference[key]) <= allowed]


def check_closure(values, rel):
    """Failure message if the canopy and snow ledger leaves more than rel * precipitation unaccounted for."""
    if abs(values["residual_m3"]) <= rel * values["prcp_m3"]:
        return []
    return ["residual_m3 %.3e exceeds %.1e of prcp_m3 %.3e: the canopy and snow ledger does not close "
            "(prcp - netprcp - canopy_evap - d_canopy - d_snow)." % (values["residual_m3"], rel, values["prcp_m3"])]


def check_purity(values):
    """Failure messages for every count of PURITY_KEYS that is not 0: f must be a function of (t, y) and the
    inputs frozen for the macro step, and must not write those inputs."""
    return ["%s=%d: %s (#16). f(t, y) must depend only on (t, y) and the inputs frozen for the macro step."
            % (key, values[key], message) for key, message in PURITY_KEYS.items() if values[key] != 0]


def compare_outputs(default_out, other_out, skip, gate, why):
    """Failure messages unless `other_out` holds the files of `default_out` with the same bytes, apart from
    *.time.csv (wall-clock times), the diagnostics file and the names in `skip`."""
    names = [{name for name in os.listdir(out)
              if name != DIAG_FILE and name not in skip and not name.endswith(".time.csv")}
             for out in (default_out, other_out)]
    failures = ["%s: written with or without %s only; %s." % (name, gate, why) for name in sorted(names[0] ^ names[1])]
    for name in sorted(names[0] & names[1]):
        with open(os.path.join(default_out, name), "rb") as a, open(os.path.join(other_out, name), "rb") as b:
            if a.read() != b.read():
                failures.append("%s: changed by %s; %s." % (name, gate, why))
    return failures


def check_diagnostics(default_out, project, days):
    """Consistency diagnostics: read-only run and repeat-check run against the invariants and the reference."""
    ledger_out = run_case("shud", project, days, extra_env={"SHUD_DIAG_CONSISTENCY": "1"})
    repeat_out = run_case("shud", project, days, extra_env={"SHUD_DIAG_RHS_REPEAT": "1"})
    values = read_keys(os.path.join(ledger_out, DIAG_FILE))
    repeat = read_keys(os.path.join(repeat_out, DIAG_FILE))
    if not values or not repeat:
        return ["%s was not written with SHUD_DIAG_CONSISTENCY=1 / SHUD_DIAG_RHS_REPEAT=1." % DIAG_FILE]
    values["rhs_repeat_differing"] = repeat["rhs_repeat_differing"]
    values["rhs_aba_differing"] = repeat["rhs_aba_differing"]
    print("diagnostics (%s %d days): residual_frac_prcp=%.6f cvode_max_overshoot_min=%.3e rhs_repeat_max_abs=%.3e"
          % (project, days, values["residual_frac_prcp"], values["cvode_max_overshoot_min"],
             repeat["rhs_repeat_max_abs"]))
    print("diagnostics (%s %d days): %s"
          % (project, days, " ".join("%s=%d" % (key, values[key]) for key in PURITY_KEYS)))
    failures = compare_outputs(default_out, ledger_out, (), "SHUD_DIAG_CONSISTENCY=1", "the gate must be read-only")
    # nfcall.txt counts the extra f() calls of the repeat check.
    failures += compare_outputs(default_out, repeat_out, ("nfcall.txt",), "SHUD_DIAG_RHS_REPEAT=1",
                                "the right-hand side depends on its call history (#16)")
    failures += check_purity(values)
    with open(os.path.join(ROOT, "input", project, project + ".cfg.para"), encoding="utf-8") as handle:
        step_min = float(re.search(r"(?m)^MAX_SOLVER_STEP\s+(\S+)", handle.read()).group(1))
    rel = constraint("regression.ledger_rel")
    for ok, message in (
            (values["macro_steps"] == days * 1440 / step_min, "macro_steps is not days * 1440 / MAX_SOLVER_STEP"),
            (repeat["rhs_repeat_calls"] == repeat["macro_steps"], "rhs_repeat_calls differs from macro_steps"),
            (values["prcp_m3"] > 0, "prcp_m3 is not positive; the ledger saw no precipitation"),
            (values["landland_edges"] > 0, "landland_edges is not positive"),
            (values["landland_sub_gross_m3"] > 0, "landland_sub_gross_m3 is not positive"),
            (values["landland_surf_gross_m3"] > 0, "landland_surf_gross_m3 is not positive")):
        if not ok:
            failures.append("%s: %s." % (DIAG_FILE, message))
    for name in INTERFACES:
        residual, gross = values[name + "_residual_m3"], values[name + "_gross_m3"]
        if not residual <= rel * gross:
            failures.append("%s_residual_m3 %.3e exceeds %.1e of the gross volume %.3e m3: the two sides of "
                            "the interface book different fluxes." % (name, residual, rel, gross))
    failures += check_closure(values, rel)
    ref_path = os.path.join(REFERENCE_DIR, project, PLATFORM, "diag_consistency.ref")
    if updating():
        if not failures:  # a run that breaks an invariant is not recorded
            with open(ref_path, "w", encoding="utf-8") as handle:
                handle.write("".join("%s=%.12e\n" % (key, values[key]) for key in LEDGER_KEYS)
                             + "".join("%s=%d\n" % (key, values[key]) for key in COUNT_KEYS))
        return failures
    reference = read_keys(ref_path)
    failures += compare_ledger(values, reference, rel)
    failures += ["%s: %d, reference %d." % (key, values[key], reference[key])
                 for key in COUNT_KEYS if values[key] != reference[key]]
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
        # Temporary tool, off by default: left out so the floor measures model code and stays put when it is removed.
        # suspicious_hits: gcovr rejects a line count above 2**32 as a gcov bug; the two-year run really gets there.
        run([tool("gcovr"), "-r", ".", "--filter", "src/", "--exclude", "src/Model/MD_consistency_diag.hpp",
             "--gcov-ignore-parse-errors", "suspicious_hits.warn_once_per_file", "--json-summary", summary])
        with open(summary, encoding="utf-8") as handle:
            percent = float(json.load(handle)["line_percent"])
    finally:
        for leftover in glob.glob(os.path.join(ROOT, "*.gcno")) + glob.glob(os.path.join(ROOT, "*.gcda")):
            os.remove(leftover)
        if os.path.exists(binary):
            os.remove(binary)
    print("line coverage (%s, %s %d days): %.1f%%" % (PLATFORM, project, days, percent))
    return check_coverage(percent, constraint("baseline.coverage_percent." + PLATFORM),
                          constraint("testing.min_line_coverage.value"))
