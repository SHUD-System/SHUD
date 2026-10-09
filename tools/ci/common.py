"""Shared helpers for the engineering gates: constraints.yaml reader, paths, subprocess."""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASELINE_DIR = os.path.join(ROOT, "tools", "ci", "baseline")
SOURCE_EXT = (".cpp", ".hpp", ".h", ".py")


def _scalar(text):
    text = text.strip()
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        return [_scalar(p) for p in inner.split(",")] if inner else []
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    if text in ("true", "false"):
        return text == "true"
    try:
        return float(text) if "." in text else int(text)
    except ValueError:
        return text


def _strip_comment(line):
    out, quote = [], None
    for ch in line:
        if quote:
            quote = None if ch == quote else quote
        elif ch in "\"'":
            quote = ch
        elif ch == "#":
            break
        out.append(ch)
    return "".join(out).rstrip()


def _parse(lines, pos, indent):
    """Parse the block of `lines` starting at `pos` whose indentation is `indent`."""
    if lines[pos][1].startswith("- "):
        items = []
        while pos < len(lines) and lines[pos][0] == indent and lines[pos][1].startswith("- "):
            body = lines[pos][1][2:]
            if re.match(r"^[A-Za-z_][\w.-]*:(\s|$)", body):
                lines[pos] = (indent + 2, body)
                value, pos = _parse(lines, pos, indent + 2)
            else:
                value, pos = _scalar(body), pos + 1
            items.append(value)
        return items, pos
    mapping = {}
    while pos < len(lines) and lines[pos][0] == indent:
        key, _, rest = lines[pos][1].partition(":")
        rest = rest.strip()
        pos += 1
        if rest:
            mapping[key.strip()] = _scalar(rest)
        elif pos < len(lines) and lines[pos][0] > indent:
            mapping[key.strip()], pos = _parse(lines, pos, lines[pos][0])
        else:
            mapping[key.strip()] = None
    return mapping, pos


def load_yaml(path):
    """Read the YAML subset used by constraints.yaml (maps, lists, scalars, inline lists)."""
    lines = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            text = _strip_comment(raw.rstrip("\n"))
            if text.strip():
                lines.append((len(text) - len(text.lstrip()), text.strip()))
    return _parse(lines, 0, 0)[0] if lines else {}


_CONSTRAINTS = None


def constraint(dotted):
    """Return constraints.yaml value at `a.b.c`; a missing key is a hard error."""
    global _CONSTRAINTS
    if _CONSTRAINTS is None:
        _CONSTRAINTS = load_yaml(os.path.join(ROOT, "constraints.yaml"))
    node = _CONSTRAINTS
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            sys.exit("constraints.yaml: missing key %s" % dotted)
        node = node[part]
    return node


def run(cmd, cwd=ROOT, check=True, env=None, capture=True):
    merged = dict(os.environ)
    merged.update(env or {})
    proc = subprocess.run(cmd, cwd=cwd, env=merged, text=True,
                          stdout=subprocess.PIPE if capture else None,
                          stderr=subprocess.STDOUT if capture else None)
    if check and proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").splitlines()[-25:])
        sys.exit("FAILED (%d): %s\n%s" % (proc.returncode, " ".join(cmd), tail))
    return proc


def tracked_files(root=ROOT):
    """Tracked plus untracked-not-ignored files, so new files are checked before `git add`."""
    out = run(["git", "ls-files", "-co", "--exclude-standard"], cwd=root).stdout
    return [f for f in out.splitlines() if os.path.isfile(os.path.join(root, f))]


def read_baseline(name):
    """Baseline file -> {key: value}. Lines are `key<TAB>value`; `#` starts a comment."""
    path = os.path.join(BASELINE_DIR, name)
    entries = {}
    if not os.path.exists(path):
        return entries
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line and not line.startswith("#"):
                key, _, value = line.partition("\t")
                entries[key] = value
    return entries


def tool(name):
    """Path of a Python tool installed in the gate virtualenv (.ci-venv)."""
    return os.path.join(ROOT, ".ci-venv", "bin", name)
