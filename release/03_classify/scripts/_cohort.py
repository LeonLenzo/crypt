"""Which runs a figure in this module analyses.

The pipeline classifies EVERY run on disk, identically and with no selection, because a
detector that sees a filtered cohort cannot be audited. The analysed set is narrowed here
instead, at figure time, where the narrowing is visible and can be stated in a caption.

**The analysed set is FIELD-COLLECTED *Triticum aestivum*.** Two groups are deliberately
retained in the cohort and excluded from the figures, and neither is noise:

- **Lab libraries** (71 runs of *T. aestivum*, 15 with no host recorded) are the axenic and
  inoculated controls that anchor the detection floors. They have no field geography and a
  pathogen in them was put there, so they cannot be read as an unreported co-infection.
- **Related cereal hosts** deposited inside cohort BioProjects (*Hordeum vulgare* 57,
  *Triticosecale* 10, *Triticum turgidum* 8, bare *Triticum* 8, *Secale cereale* 5,
  *Hordeum jubatum* 4, and three singletons) are the HOST CONTRAST. They share a deposit and
  a protocol with the wheat runs, which is exactly what makes them a contrast, but a
  *Ramularia* on barley is not an observation about wheat.

Bare `Triticum` is excluded with the rest: it is 8 runs whose host was never resolved below
genus, so calling them *T. aestivum* would be an assumption rather than a reading.

This removes 181 of 1,560 runs and 7 BioProjects, 6 of which carry no curated declaration
anyway. See crypt-cohort-is-wheat and crypt-clean-pipeline-runs.
"""

HOST = "Triticum aestivum"
SETTING = "field"


def is_analysed(row):
    """True if this `supp_data_s2_cohort.tsv` row is in the analysed set."""
    return row["setting"] == SETTING and row["host"] == HOST


def describe():
    return f"field-collected {HOST}"
