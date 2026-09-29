# Complexity by host — per-species detail

One image per pathogen species (61), each showing reads vs distinct genome fragments
(k-mers) with one panel per host. This is the per-species detail behind the
tightness model in `../true_detections/`.

**How to read each plot:** each panel is one host. Points that rise **tightly** along
the fit (R² ≥ 0.5, blue, ✓ real) are a genuine infection — coverage grows cleanly with
sequencing depth. A **flat or scattered** fit (red, ✗ artefact) is a pile-up on one
fragment or cross-mapping from a relative. Tightness, not the amount of signal, so
genuine low-abundance infections are kept. Hosts come from the read-based call
(`kraken_host_call.py`), which resolves every run — no "unresolved".

The same organism can be real on one host and an artefact on another — see
`Bipolaris_maydis.png` (real on maize and barley, artefact on wheat).

*Regenerate: `kraken/assign/figures/prep_by_host_individual.py` then
`complexity_by_host_individual.R`.*
