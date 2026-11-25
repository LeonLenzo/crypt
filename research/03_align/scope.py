#!/usr/bin/env python3
"""What job 1 is aimed at and how big it is, computed rather than quoted.

    python research/03_align/scope.py

Reads  release/03_classify/data/outputs/supp_data_s5_detections.tsv   the db2 calls
       release/03_classify/data/outputs/supp_data_s6_runs.tsv         the cluster flag
       release/01_cohort/data/outputs/supp_data_s2_cohort.tsv         base_count per run
       release/02_reference/data/outputs/supp_data_s3_reference.tsv   what the index can hold
Writes nothing. Prints the sizing the README quotes.

A cluster here is >= 3 congeners of one genus called in one run, which is the shape the
classifier cannot adjudicate: the genus' k-mer space is shared by construction, so uniformity
across congeners is a diagnostic rather than a co-infection. The CV printed beside it is the
within-run spread of coverage over the genus, and Puccinia is the contrast, three species that
are separable by exact k-mers across a 16 point ANI gap.

THE READ VOLUME IS THE POINT OF THE LAST BLOCK. bowtie2 -k over this many bases is the
module's first real cost and nothing in the project has measured it, so this prints the volume
and not a runtime. Size the SLURM request from a timed dry run on three runs of different
depth.
"""

import collections
import csv
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLS = ROOT / "release/03_classify/data/outputs"
ARM = "db2"
MIN_CLUSTER = 3


def rd(path):
    with open(path, newline="", encoding="utf-8") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main():
    by_run = collections.defaultdict(lambda: collections.defaultdict(list))
    by_genus = collections.defaultdict(list)
    for r in rd(CLS / "supp_data_s5_detections.tsv"):
        if r["arm"] != ARM or r["taxon_class"] == "host":
            continue
        by_run[r["run"]][r["genus"]].append(r)
        by_genus[r["genus"]].append(r)

    flagged = [r for r in rd(CLS / "supp_data_s6_runs.tsv")
               if r["congener_clusters"] not in ("", "0")]
    print(f"runs with a congener cluster flagged: {len(flagged)}")
    print("  by genus:", collections.Counter(
        r["congener_clusters"] for r in flagged).most_common(6))

    clusters = collections.Counter()
    species = collections.defaultdict(set)
    for run, genera in by_run.items():
        for genus, rows in genera.items():
            if len(rows) >= MIN_CLUSTER:
                clusters[genus] += 1
                species[genus].update(r["taxon"] for r in rows)
    print(f"\ngenera with >= {MIN_CLUSTER} congeners called in one run:")
    for genus, n in clusters.most_common(6):
        cv = [float(r["genus_cov_cv"]) for r in by_genus[genus]
              if r["genus_cov_cv"] not in ("", "nan")]
        med = f"{statistics.median(cv):.3f}" if cv else "n/a"
        print(f"  {genus:<14} {n:>4} runs  {len(species[genus]):>3} species  "
              f"median within-run coverage CV {med}")

    target = {run for run, genera in by_run.items()
              if len(genera.get("Alternaria", [])) >= MIN_CLUSTER}
    coh = {r["run"]: r for r in rd(ROOT / "release/01_cohort/data/outputs/"
                                         "supp_data_s2_cohort.tsv")}
    bases = [int(coh[r]["base_count"]) for r in target
             if coh.get(r, {}).get("base_count")]
    print(f"\nJOB 1 READS: {len(target)} runs with an Alternaria cluster, "
          f"base_count known for {len(bases)}")
    print(f"  {sum(bases) / 1e9:.1f} Gbases total, mean {statistics.mean(bases) / 1e9:.2f}")
    print("  by project:", collections.Counter(
        coh[r]["BioProject"] for r in target if r in coh).most_common(6))

    s3 = list(rd(ROOT / "release/02_reference/data/outputs/supp_data_s3_reference.tsv"))
    for label, rows in (("Pleosporaceae", [r for r in s3 if r["family"] == "Pleosporaceae"]),
                        ("Alternaria", [r for r in s3 if r["genus"] == "Alternaria"])):
        annotated = sum(1 for r in rows if r["n_annotated"] not in ("", "0"))
        predicted = sum(1 for r in rows if r["n_predicted"] not in ("", "0"))
        print(f"JOB 1 INDEX: {len(rows):>4} {label} taxa in the reference, "
              f"{annotated} NCBI-annotated, {predicted} predicted here")

    for taxon in ("Alternaria triticina", "Stemphylium lycopersici"):
        n = sum(len([x for x in rows if x["taxon"] == taxon])
                for rows in (by_genus[taxon.split()[0]],))
        print(f"TIER A: {taxon:<24} {n} {ARM} calls")
    return 0


if __name__ == "__main__":
    sys.exit(main())
