#!/usr/bin/env python3
"""Map which script reads and writes which file, by parsing the code.

    python tools/dataflow.py                  # write the manifest, print the findings
    python tools/dataflow.py --reads X        # who reads a file
    python tools/dataflow.py --writes X       # who writes it

Writes tools/dataflow.tsv: one row per (script, direction, file).

These scripts are coupled through TSV files on disk, not through imports, so an import or
call graph does not show the pipeline. This does. It needs no execution and cannot drift
from the code, unlike a hand-written diagram.

Paths are resolved through the `_layout.py` constants so `READS_MANIFEST` and
`SILVER / "reads_manifest.tsv"` are recognised as the same file. Each script's docstring is
read as well, and a declared "Reads X / Writes Y" that the code contradicts is reported: the
docstring is what a person believes, the AST is what runs.
"""

import argparse
import ast
import collections
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "tools" / "dataflow.tsv"
SKIP = {"superseded", "__pycache__", ".git", "graphify-out", "site", "studies"}
DATA_EXT = (".tsv", ".csv", ".json", ".txt", ".jsonl", ".db", ".sqlite", ".fasta",
            ".fna", ".report", ".xml", ".k2d", ".gz")
WRITE_MODES = ("w", "a", "x")
# Calls that take a path but do no I/O. Without this, `print(f"-> {OUT.name}")` alone
# produced 50 spurious edges.
NOT_IO = {"print", "len", "str", "int", "float", "sorted", "set", "list", "dict", "tuple",
          "exit", "format", "repr", "sum", "max", "min", "any", "all", "zip", "enumerate",
          "range", "isinstance", "getattr", "setattr", "relative_to", "name", "exists",
          "is_file", "is_dir", "mkdir", "stem", "suffix", "parent", "resolve", "joinpath",
          "SystemExit", "ValueError", "RuntimeError", "sys", "f", "A", "B"}
COLS = ["script", "direction", "file", "expression", "line", "found_by"]


def constants_in(tree, src):
    """ALL_CAPS constant -> filename, for one module's own top-level assignments.

    Most scripts define their own OUT, FUNNEL, SILVER and so on rather than importing them,
    so resolving only the shared _layout.py constants missed the majority of the pipeline:
    98 edges instead of the real figure, and build_s3.py reading depth_funnel.tsv through a
    local FUNNEL constant looked like a dead output.
    """
    out = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        tgt = node.targets[0]
        if not isinstance(tgt, ast.Name) or not tgt.id.isupper():
            continue
        expr = ast.get_source_segment(src, node.value) or ""
        names = re.findall(r'"([^"]+)"', expr) + re.findall(r"'([^']+)'", expr)
        hit = next((n for n in reversed(names) if n.endswith(DATA_EXT)), None)
        if hit:
            out[tgt.id] = hit
    return out


def layout_constants():
    """ALL_CAPS constant -> filename, from every _layout.py / _paths.py in the tree.

    Read with the AST rather than imported: importing executes module-level code, and some
    of these modules touch the filesystem when loaded.
    """
    out = {}
    for p in ROOT.rglob("_*.py"):
        if p.name not in ("_layout.py", "_paths.py") or any(s in p.parts for s in SKIP):
            continue
        try:
            tree = ast.parse(p.read_text())
        except SyntaxError:
            continue
        src = p.read_text()
        for node in tree.body:
            if not isinstance(node, ast.Assign) or len(node.targets) != 1:
                continue
            t = node.targets[0]
            if not isinstance(t, ast.Name) or not t.id.isupper():
                continue
            expr = ast.get_source_segment(src, node.value) or ""
            names = re.findall(r'"([^"]+)"', expr) + re.findall(r"'([^']+)'", expr)
            hit = next((n for n in reversed(names) if n.endswith(DATA_EXT)), None)
            if hit:
                out[t.id] = hit
    return out


def target_of(expr, consts):
    """The filename an expression refers to, for matching one script's path against another's."""
    names = re.findall(r'"([^"]+)"', expr) + re.findall(r"'([^']+)'", expr)
    hit = next((n for n in reversed(names) if n.endswith(DATA_EXT)), None)
    if hit:
        return Path(hit).name
    for c in re.findall(r"\b([A-Z][A-Z0-9_]+)\b", expr):
        if c in consts:
            return Path(consts[c]).name
    return ""


def mode_of(call):
    """The mode string of an open() call, defaulting to read."""
    if len(call.args) > 1 and isinstance(call.args[1], ast.Constant):
        return str(call.args[1].value)
    for kw in call.keywords:
        if kw.arg == "mode" and isinstance(kw.value, ast.Constant):
            return str(kw.value.value)
    return "r"


def scan(path, consts):
    src = path.read_text(errors="ignore")
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return []
    rows, rel = [], str(path.relative_to(ROOT))
    # A module's own constants win over the shared ones where both define a name.
    consts = {**consts, **constants_in(tree, src)}
    # A parameter whose default is a resolving constant carries that file inside the
    # function body: `def manifest(manifest_path=READS_MANIFEST)` then writes through
    # `manifest_path`. Without this, every helper-mediated path is invisible.
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            a = node.args
            for arg, default in zip(a.args[len(a.args) - len(a.defaults):], a.defaults):
                expr = ast.get_source_segment(src, default) or ""
                tgt = target_of(expr, consts)
                if tgt:
                    consts.setdefault(arg.arg, tgt)

    def add(direction, node, expr, found_by):
        tgt = target_of(expr, consts)
        if tgt:
            rows.append({"script": rel, "direction": direction, "file": tgt,
                         "expression": " ".join(expr.split())[:70],
                         "line": getattr(node, "lineno", ""), "found_by": found_by})

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        name = f.id if isinstance(f, ast.Name) else (f.attr if isinstance(f, ast.Attribute) else "")
        if name == "open" and node.args:
            expr = ast.get_source_segment(src, node.args[0]) or ""
            d = "write" if mode_of(node)[:1] in WRITE_MODES else "read"
            add(d, node, expr, "open")
        elif name in ("atomic_write",) and node.args:
            add("write", node, ast.get_source_segment(src, node.args[0]) or "", name)
        elif name in ("write_text", "write_bytes", "to_csv", "to_parquet"):
            expr = (ast.get_source_segment(src, f.value) if isinstance(f, ast.Attribute) else "")
            if name.startswith("to_") and node.args:
                expr = ast.get_source_segment(src, node.args[0]) or ""
            add("write", node, expr or "", name)
        elif name in ("read_csv", "read_table", "read_json") and node.args:
            add("read", node, ast.get_source_segment(src, node.args[0]) or "", name)
        elif name and name not in NOT_IO and node.args:
            # Any other call handed something that resolves to a data file. Direction is
            # UNKNOWN, so it is labelled `touches` rather than guessed: "no readers" is a
            # conclusion, "could not determine" is not, and reporting the second as the
            # first is how this manifest would mislead. `write_manifest(tree, path)` in
            # classify.py is the case that forced this.
            for a in node.args:
                expr = ast.get_source_segment(src, a) or ""
                if target_of(expr, consts):
                    add("touches", node, expr, name)

    # What the docstring claims, for comparison with what the code does.
    doc = ast.get_docstring(tree) or ""
    for m in re.finditer(r"^\s*(Reads|Writes)\s+(\S+)", doc, re.M | re.I):
        tgt = target_of(m.group(2), consts) or Path(m.group(2)).name
        if tgt.endswith(DATA_EXT):
            rows.append({"script": rel, "direction": m.group(1).lower().rstrip("s"),
                         "file": tgt, "expression": m.group(2), "line": "",
                         "found_by": "docstring"})
    return rows


# R has no AST parser here, and the figures are all R, so these are matched textually.
# Without them every figure's data table looked like a dead output: prep_depth_funnel.py
# writes depth_funnel.tsv and only depth_funnel.R reads it.
# R has no AST parser here, and every figure is R, so these are matched textually. Without
# them each figure's data table looked like a dead output: prep_depth_funnel.py writes
# depth_funnel.tsv and only depth_funnel.R reads it.
#   The whole call is captured, not its first argument. `read.delim(file.path(HERE,
# "depth_funnel.tsv"))` puts the filename in file.path's SECOND argument, so a pattern that
# stops at the first comma captures `file.path(HERE` and finds nothing.
R_CALL = re.compile(r"\b(read\.delim|read\.csv|read\.table|read_tsv|read_csv|fread|"
                    r"write\.csv|write\.table|write_tsv|write_csv|ggsave)\s*\(([^\n]*)",
                    re.I)


def scan_r(path, consts):
    src = path.read_text(errors="ignore")
    rel, rows = str(path.relative_to(ROOT)), []
    local = dict(consts)
    for m in re.finditer(r"(\w+)\s*<-\s*([^\n]+)", src):
        names = re.findall(r'"([^"]+)"', m.group(2)) + re.findall(r"'([^']+)'", m.group(2))
        hit = next((n for n in reversed(names) if n.endswith(DATA_EXT)), None)
        if hit:
            local[m.group(1)] = hit
    seen = set()
    for m in R_CALL.finditer(src):
        fn, rest = m.group(1).lower(), m.group(2)
        direction = "read" if fn.startswith(("read", "fread")) else "write"
        lits = [n for n in re.findall(r'"([^"]+)"', rest) + re.findall(r"'([^']+)'", rest)
                if n.endswith(DATA_EXT) or n.endswith((".png", ".pdf", ".svg"))]
        for v in re.findall(r"\b(\w+)\b", rest):
            if v in local and local[v] not in lits:
                lits.append(local[v])
        for name in lits:
            tgt = Path(name).name
            if (direction, tgt) in seen:
                continue
            seen.add((direction, tgt))
            rows.append({"script": rel, "direction": direction, "file": tgt,
                         "expression": " ".join(rest.split())[:70],
                         "line": src[:m.start()].count("\n") + 1, "found_by": "R"})
    return rows


def build():
    consts = layout_constants()
    rows = []
    for p in sorted(ROOT.rglob("*.py")):
        if any(s in p.parts for s in SKIP):
            continue
        rows.extend(scan(p, consts))
    for p in sorted(list(ROOT.rglob("*.R")) + list(ROOT.rglob("*.r"))):
        if any(s in p.parts for s in SKIP):
            continue
        rows.extend(scan_r(p, consts))
    seen, uniq = set(), []
    for r in rows:
        k = (r["script"], r["direction"], r["file"])
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(uniq)
    print(f"-> {OUT.relative_to(ROOT)}  {len(uniq)} edges, "
          f"{len({r['script'] for r in uniq})} scripts, {len({r['file'] for r in uniq})} files")
    print(f"   {len(consts)} path constants resolved from _layout.py")
    return uniq


def report(rows):
    readers = collections.defaultdict(set)
    writers = collections.defaultdict(set)
    touchers = collections.defaultdict(set)
    for r in rows:
        if r["direction"] == "write":
            writers[r["file"]].add(r["script"])
        elif r["direction"] == "read":
            readers[r["file"]].add(r["script"])
        else:
            touchers[r["file"]].add(r["script"])

    print("\n  files with MORE THAN ONE writer (who wins?):")
    multi = {f: w for f, w in writers.items() if len(w) > 1}
    for f, w in sorted(multi.items(), key=lambda x: -len(x[1]))[:8]:
        print(f"    {f:<34} {len(w)} writers: {', '.join(sorted(w))[:88]}")
    if not multi:
        print("    none")

    print("\n  files read but NEVER written here (external inputs):")
    ext = sorted(set(readers) - set(writers))
    print(f"    {len(ext)}: {', '.join(ext[:12])}")

    # Not "dead": a manuscript table, a status log and a figure all legitimately end here.
    # It means no SCRIPT consumes it, which is a question to ask, not a verdict.
    print("\n  files no script consumes (terminal artefacts, or genuinely dead):")
    dead = sorted(set(writers) - set(readers) - set(touchers))
    print(f"    {len(dead)}: {', '.join(dead[:12])}")
    unclear = sorted((set(writers) & set(touchers)) - set(readers))
    print(f"  written, and only TOUCHED by something else (direction undetermined): "
          f"{len(unclear)}")

    print("\n  scripts whose every output is read by nobody (candidate dead ends):")
    out_by = collections.defaultdict(set)
    for r in rows:
        if r["direction"] == "write":
            out_by[r["script"]].add(r["file"])
    orph = [s for s, fs in out_by.items()
            if fs and not any(readers.get(f) or touchers.get(f) for f in fs)]
    for s in sorted(orph)[:14]:
        print(f"    {s}")
    print(f"    ({len(orph)} of {len(out_by)} writing scripts)")

    print("\n  declared in a docstring but not found in the code:")
    code = {(r["script"], r["direction"], r["file"]) for r in rows if r["found_by"] != "docstring"}
    gap = [r for r in rows if r["found_by"] == "docstring"
           and (r["script"], r["direction"], r["file"]) not in code]
    for r in gap[:10]:
        print(f"    {r['script']}  {r['direction']}s {r['file']}")
    print(f"    ({len(gap)} such claims)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", help="which scripts read this file")
    ap.add_argument("--writes", help="which scripts write it")
    args = ap.parse_args()

    rows = build()
    if args.reads or args.writes:
        want, d = (args.reads, "read") if args.reads else (args.writes, "write")
        key = Path(want).name
        hits = [r for r in rows if r["file"] == key and r["direction"] == d]
        maybe = [r for r in rows if r["file"] == key and r["direction"] == "touches"]
        print(f"\n  {len(hits)} script(s) {d} {key}:")
        for r in sorted(hits, key=lambda x: x["script"]):
            print(f"    {r['script']}:{r['line']}  {r['expression']}")
        if maybe:
            # Shown because "nobody" and "could not tell" must not look the same.
            print(f"  plus {len(maybe)} that reference it through a helper, "
                  "direction undetermined:")
            for r in sorted(maybe, key=lambda x: x["script"]):
                print(f"    {r['script']}:{r['line']}  via {r['found_by']}()")
        if not hits and not maybe:
            print("    nothing found. A path built inline rather than from a constant "
                  "is invisible here.")
        return
    report(rows)


if __name__ == "__main__":
    main()
