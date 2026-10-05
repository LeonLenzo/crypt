#!/usr/bin/env python3
"""kallisto_select.py — pick the stratified set of detections to test, and so the module's scope.

Do not quant the whole cohort. This chooses a set that spans the ways a detection can be real
or artefactual, which doubles as the calibration set for the 03_kraken/05_filter thresholds
(`score`, `coverage`) that nothing has ever been able to set objectively.

Six strata, each a class the three read-level facts should separate:

  strong      high score and coverage, no measured competitor in the same run -> confirm
  borderline  coverage within 2x of the floor -> what the floor gets calibrated on
  congeneric  >=2 species of one genus flagged in the same run -> collapse if cross-mapping
  off_host    the species' flagged runs sit overwhelmingly on other hosts -> reject
  discordant  coverage high but score below the strong bar: the genome looks largely present
              while the accumulation shape looks marginal. The two filter criteria disagree
              here, which makes it the most informative class for calibrating them. 1,630 of
              the 3,980 flagged detections, dominated by P. striiformis.
  mid         coverage between the borderline band and the strong bar, at neither end.

The first four span the ends and the known artefact mechanisms, but on their own they label
only 45% of the flagged set; `discordant` and `mid` close that gap so the calibration set can
speak for the whole flagged population rather than its extremes.

A detection can carry several labels; all are kept, because a congeneric borderline case is
exactly the kind the two criteria disagree about.

Quant is per run, so a run carries every detection it holds. Selection is deterministic (no
RNG) and reproducible from the committed inputs.

**Every flagged host is represented.** Each stratum first takes one detection from each host
that qualifies, then fills the remainder round-robin over host x species. Without the first
pass a flat quota silently drops whole hosts: at 40 per stratum it lost Zea mays, the second
largest with 421 flagged detections.

    python 04_kallisto/01_select/kallisto_select.py
    python 04_kallisto/01_select/kallisto_select.py --n-per 60
"""
import argparse
import collections
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import READS_MAN, STRATA, flagged, num, step_log

BORDER_MAX    = 0.02   # within 2x of the coverage floor
STRONG_COV    = 0.10
STRONG_SCORE  = 0.90
OFFHOST_SHARE = 0.10   # this host holds <10% of the species' flagged runs
OFFHOST_MIN_N = 20     # ...and it has enough runs for a dominant host to be real
STRATA_ORDER  = ("strong", "borderline", "congeneric", "off_host", "discordant", "mid")
EXPECT = {"strong": "confirm", "congeneric": "collapse-if-cross-map",
          "off_host": "reject", "borderline": "unknown-calibrates-floor",
          "discordant": "unknown-criteria-conflict", "mid": "unknown-mid-range"}


def read_paths():
    """run -> [(mate, path, mtime)] from the committed 03_select manifest.

    The manifest is tracked, so read layout and purge exposure are answerable locally without
    touching Setonix.
    """
    pat = re.compile(r"reads/([A-Z]RR\d+)(?:_([12]))?\.fastq(?:\.gz)?$")
    out = collections.defaultdict(list)
    if not READS_MAN.exists():
        return out
    for r in csv.DictReader(open(READS_MAN), delimiter="\t"):
        m = pat.match(r["path"])
        if m:
            out[m.group(1)].append((m.group(2) or "0", r["path"], int(r["mtime_epoch"])))
    for k in out:
        out[k].sort()
    return out


def label(dets):
    """Attach the strata each detection belongs to."""
    per_run = collections.defaultdict(list)
    for d in dets:
        per_run[(d["host"], d["run"])].append(d)

    sp_host = collections.defaultdict(collections.Counter)   # species -> host -> n runs
    for (h, _), rs in per_run.items():
        for sp in {d["species"] for d in rs}:
            sp_host[sp][h] += 1

    for (h, _), rs in per_run.items():
        genus_n = collections.Counter(d["species"].split()[0] for d in rs)
        for d in rs:
            cov, sc = num(d, "coverage"), num(d, "score")
            competing = genus_n[d["species"].split()[0]] > 1
            tags = ["congeneric"] if competing else []
            total = sum(sp_host[d["species"]].values())
            if total >= OFFHOST_MIN_N and sp_host[d["species"]][h] / total < OFFHOST_SHARE:
                tags.append("off_host")
            if cov < BORDER_MAX:
                tags.append("borderline")
            if cov >= STRONG_COV and sc >= STRONG_SCORE and not competing:
                tags.append("strong")
            if cov >= STRONG_COV and sc < STRONG_SCORE:
                tags.append("discordant")
            if BORDER_MAX <= cov < STRONG_COV:
                tags.append("mid")
            d["strata"] = tags
    return dets


def _spread(items, n):
    """n items evenly spaced along coverage, so a stratum covers its span, not its head."""
    items = sorted(items, key=lambda d: num(d, "coverage"))
    if len(items) <= n:
        return items
    step = (len(items) - 1) / (n - 1) if n > 1 else 0
    return [items[round(i * step)] for i in range(n)]


def pick(dets, stratum, n_per, have_reads):
    """n_per detections for one stratum, guaranteeing every qualifying host appears."""
    pool = [d for d in dets if stratum in d["strata"] and d["run"] in have_reads]
    if not pool:
        return []
    by_host = collections.defaultdict(list)
    for d in pool:
        by_host[d["host"]].append(d)

    # pass 1: one per host, mid-coverage, so no host can be squeezed out by the quota
    chosen, seen = [], set()
    for h in sorted(by_host):
        d = _spread(by_host[h], 3)[len(_spread(by_host[h], 3)) // 2]
        chosen.append(d)
        seen.add((d["run"], d["species"]))

    # pass 2: fill the remainder round-robin over host x species, spread along coverage
    groups = collections.defaultdict(list)
    for d in pool:
        if (d["run"], d["species"]) not in seen:
            groups[(d["host"], d["species"])].append(d)
    order = sorted(groups)
    while len(chosen) < n_per and order:
        progressed = False
        for g in list(order):
            if len(chosen) >= n_per:
                break
            cand = [d for d in groups[g] if (d["run"], d["species"]) not in seen]
            if not cand:
                order.remove(g)
                continue
            d = _spread(cand, 1)[0]
            chosen.append(d)
            seen.add((d["run"], d["species"]))
            progressed = True
        if not progressed:
            break
    return chosen


def select(n_per=40, out=STRATA):
    dets = label(flagged())
    have_reads = read_paths()
    picked = {}
    for st in STRATA_ORDER:
        for d in pick(dets, st, n_per, have_reads):
            picked.setdefault((d["run"], d["species"]),
                              dict(d, picked_for=[]))["picked_for"].append(st)

    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["stratum", "expectation", "host", "run", "biosample", "species", "reads",
            "coverage", "score", "all_strata", "layout", "read_files"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for (run, sp), d in sorted(picked.items(),
                                   key=lambda kv: (kv[1]["picked_for"][0], kv[0])):
            rf = have_reads.get(run, [])
            layout = "paired" if {m for m, _, _ in rf} >= {"1", "2"} else "single"
            fh.write("\t".join([",".join(d["picked_for"]), EXPECT[d["picked_for"][0]],
                                d["host"], run, d["biosample"], sp, d["reads"],
                                d["coverage"], d["score"], ",".join(d["strata"]), layout,
                                ",".join(p for _, p, _ in rf)]) + "\n")

    hosts = {d["host"] for d in picked.values()}
    all_hosts = {d["host"] for d in dets}
    print(f"{out}: {len(picked)} detections, {len({r for r, _ in picked})} runs, "
          f"{len(hosts)} hosts")
    for st in STRATA_ORDER:
        sel = [d for d in picked.values() if st in d["picked_for"]]
        print(f"  {st:11s} {len(sel):4d} detections  {len({d['run'] for d in sel}):4d} runs  "
              f"{len({d['species'] for d in sel}):3d} species  "
              f"{len({d['host'] for d in sel}):2d} hosts  -> expect {EXPECT[st]}")
    missed = sorted(all_hosts - hosts)
    print(f"\nhost coverage: {len(hosts)}/{len(all_hosts)} flagged hosts"
          + (f"\n  MISSING: {', '.join(missed)}" if missed else " (all represented)"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-per", type=int, default=40, help="detections per stratum (default 40)")
    ap.add_argument("--out", type=Path, default=STRATA)
    args = ap.parse_args()
    log = step_log(__file__, "kallisto_select")
    try:
        select(args.n_per, args.out)
    finally:
        log.close()


if __name__ == "__main__":
    main()
