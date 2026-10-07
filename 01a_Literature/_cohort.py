#!/usr/bin/env python3
"""THE definition of the cohort. One place, so two answers are impossible.

Added 2026-10-06. The three predicates below were being retyped into every ad-hoc query, and
on that day alone they were written out four times. One copy differed, which produced a
locality count of 425 in one script and 424 in another and cost three tool calls to chase a
discrepancy that meant nothing. A definition used in four places belongs in one.

A run is in the cohort when all three hold:

    setting starts with `field`      the TISSUE was field-grown. A `field*` setting covers
                                    excised-and-forced and transferred-outdoors variants,
                                    which are field samples because the organisms in them
                                    were acquired in the field (leon 2026-10-06).
    library_representative != "no"   one row per biological sample. Only assessed where
                                    duplication was actually found, so most rows are blank
                                    and blank means "representative".
    tissue is not root-like          aerial tissue only. Soil community makes a fungal hit
                                    in root material as likely a saprotroph or endophyte as
                                    a pathogen. Also drops the 3 insect runs.

What is deliberately NOT a criterion: `sampling_selection`. Every stratum is in the cohort;
the flag decides what counts as a finding, not whether the sample counts at all. See
triage.py's sampling_selection comment and crypt-two-evidence-standards.

The counts and `gold/cohort.tsv` come from `report/cohort.py`, which is the CLI over this.
"""

from __future__ import annotations

import csv, re
from pathlib import Path

from _layout import RUNS

NON_AERIAL = re.compile(r"root|rhizo|tuber|nodule|whole.?(plant|seedling)|insect", re.I)


def is_cohort(r: dict) -> bool:
    return (r.get("setting", "").startswith("field")
            and r.get("library_representative", "") != "no"
            and not NON_AERIAL.search(r.get("tissue", "") or ""))


def all_runs(path: Path = RUNS) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cohort(runs: list[dict] | None = None) -> list[dict]:
    return [r for r in (runs if runs is not None else all_runs()) if is_cohort(r)]
