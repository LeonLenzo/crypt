"""
Four-panel comparison: masked vs unmasked Kraken2 DB, split by host and pathogen.

Output:
  masked_vs_unmasked.png   – 2×2 panel figure
  masked_vs_unmasked_caption.txt

Inputs:
  /tmp/setonix_kraken_metrics.json   – extracted from Setonix (re-run extract if stale)
  output/00_build/data/phibase_db.json

Run from crypt/:
  python output/figures/db_comparison/compare_host_pathogen.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── paths ─────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).parents[3]
KRAKEN_JSON = Path("/tmp/setonix_kraken_metrics.json")
PHIBASE     = ROOT / "output/00_build/data/phibase_db.json"
OUT_DIR     = Path(__file__).parent

# ── load reference DB for species classification ──────────────────────────────
db = json.loads(PHIBASE.read_text())
pathogen_taxids = (
    set(db["fungal_to_seed"])
    | set(db["oomycete_to_seed"])
    | set(db["nematode_to_seed"])
)
host_taxids   = set(db["host_to_seed"])
name_to_taxid = db["name_to_taxid"]
viridiplantae = set(db["viridiplantae_names"])

def classify(name: str) -> str:
    low = name.lower()
    tid = str(name_to_taxid.get(low, ""))
    if tid in pathogen_taxids:
        return "pathogen"
    if tid in host_taxids or low in viridiplantae:
        return "host"
    return "unknown"

# ── load Kraken data ──────────────────────────────────────────────────────────
kraken = json.loads(KRAKEN_JSON.read_text())
runs   = sorted(kraken)

def summarise(species_dict, kind):
    pct, n = 0.0, 0
    for name, p in species_dict.items():
        if classify(name) == kind:
            pct += p
            n   += 1
    return pct, n

rows = {}
for run in runs:
    m, u = kraken[run]["masked_species"], kraken[run]["unmasked_species"]
    rows[run] = {
        "m_host_pct":  summarise(m, "host")[0],
        "u_host_pct":  summarise(u, "host")[0],
        "m_host_sp":   summarise(m, "host")[1],
        "u_host_sp":   summarise(u, "host")[1],
        "m_path_pct":  summarise(m, "pathogen")[0],
        "u_path_pct":  summarise(u, "pathogen")[0],
        "m_path_sp":   summarise(m, "pathogen")[1],
        "u_path_sp":   summarise(u, "pathogen")[1],
        "top_host":    next((n for n in m if classify(n) == "host"),     "Unknown"),
        "top_path":    next((n for n in m if classify(n) == "pathogen"), "Unknown"),
    }

# ── colour maps ───────────────────────────────────────────────────────────────
HOST_COLOURS = {
    "Triticum aestivum":       "#c9a227",
    "Triticum turgidum":       "#e9c46a",
    "Cicer arietinum":         "#e76f51",
    "Zea mays":                "#52b788",
    "Hordeum vulgare":         "#a7c957",
    "Glycine max":             "#6a994e",
    "Medicago truncatula":     "#457b9d",
    "Brachypodium distachyon": "#8ecae6",
    "Arachis hypogaea":        "#f4a261",
}
PATHOGEN_COLOURS = {
    "Puccinia striiformis":     "#e6700a",
    "Puccinia graminis":        "#9b2226",
    "Puccinia triticina":       "#d62828",
    "Zymoseptoria tritici":     "#0077b6",
    "Sclerotinia sclerotiorum": "#2d6a4f",
    "Alternaria alternata":     "#7b2d8b",
    "Exserohilum turcicum":     "#3a86ff",
    "Monilinia fructicola":     "#f4a261",
    "Fusarium graminearum":     "#c77dff",
}
DEFAULT = "#aaaaaa"

def colours_patches(key, cmap):
    clrs = [cmap.get(rows[r][key], DEFAULT) for r in runs]
    seen = {}
    for r in runs:
        seen.setdefault(cmap.get(rows[r][key], DEFAULT), rows[r][key])
    ptchs = [mpatches.Patch(color=c, label=lbl) for c, lbl in seen.items()]
    return clrs, ptchs

hc, hp = colours_patches("top_host", HOST_COLOURS)
pc, pp = colours_patches("top_path", PATHOGEN_COLOURS)

# ── data arrays ───────────────────────────────────────────────────────────────
mhp = np.array([rows[r]["m_host_pct"]  for r in runs])
uhp = np.array([rows[r]["u_host_pct"]  for r in runs])
mhs = np.array([rows[r]["m_host_sp"]   for r in runs], dtype=float)
uhs = np.array([rows[r]["u_host_sp"]   for r in runs], dtype=float)
mpp = np.array([rows[r]["m_path_pct"]  for r in runs])
upp = np.array([rows[r]["u_path_pct"]  for r in runs])
mps = np.array([rows[r]["m_path_sp"]   for r in runs], dtype=float)
ups = np.array([rows[r]["u_path_sp"]   for r in runs], dtype=float)

print(f"n = {len(runs)} runs")
print(f"Host    reads:  masked {mhp.mean():.1f}%  unmasked {uhp.mean():.1f}%")
print(f"Host    species: masked {mhs.mean():.1f}   unmasked {uhs.mean():.1f}")
print(f"Pathogen reads: masked {mpp.mean():.1f}%  unmasked {upp.mean():.1f}%")
print(f"Pathogen sp:    masked {mps.mean():.1f}   unmasked {ups.mean():.1f}")

# ── scatter helper ─────────────────────────────────────────────────────────────
def panel(ax, x, y, clrs, xlabel, ylabel, lim=None):
    lim = lim or (max(x.max(), y.max()) * 1.08 or 1)
    ax.plot([0, lim], [0, lim], "--", color="#999999", lw=0.9, zorder=1)
    ax.scatter(x, y, c=clrs, s=42, alpha=0.88,
               edgecolors="white", linewidths=0.4, zorder=2)
    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.tick_params(labelsize=7.5)
    if np.std(x) > 0 and np.std(y) > 0:
        r = np.corrcoef(x, y)[0, 1]
        ax.text(0.97, 0.04, f"r = {r:.2f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="#666666")

# ── 2×2 figure ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.5))
fig.subplots_adjust(hspace=0.45, wspace=0.30, left=0.09, right=0.72,
                    top=0.93, bottom=0.08)

panel_labels = [["A", "B"], ["C", "D"]]
row_titles   = ["Host species", "Pathogen species"]
col_titles   = ["Reads assigned (%)", "Species detected"]

clrs = [[hc, hc], [pc, pc]]
xs   = [[mhp, mhs], [mpp, mps]]
ys   = [[uhp, uhs], [upp, ups]]

for ri in range(2):
    for ci in range(2):
        ax = axes[ri, ci]
        panel(ax, xs[ri][ci], ys[ri][ci], clrs[ri][ci],
              "Masked (%)" if ci == 0 else "Masked (# spp.)",
              "Unmasked (%)" if ci == 0 else "Unmasked (# spp.)",
              lim=100)
        # panel label (A–D) top-left
        ax.text(0.03, 0.97, panel_labels[ri][ci], transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")
        # column title on top row only
        if ri == 0:
            ax.set_title(col_titles[ci], fontsize=10, fontweight="bold", pad=5)
        # y-axis label on left column only
        if ci == 0:
            ax.set_ylabel("Unmasked (%)", fontsize=8.5)

# bold row labels to the left of each row (mirrors column title style)
for ri, label in enumerate(row_titles):
    # get the centre y of the row in figure coords
    y0 = axes[ri, 0].get_position().y0
    y1 = axes[ri, 0].get_position().y1
    fig.text(0.01, (y0 + y1) / 2, label, va="center", ha="left",
             rotation="vertical", fontsize=10, fontweight="bold", color="#222222")

# legends to the right of the figure
fig.legend(handles=hp, loc="upper left", fontsize=8, frameon=False,
           title="Host species", title_fontsize=8.5,
           bbox_to_anchor=(0.74, 0.95))
fig.legend(handles=pp, loc="upper left", fontsize=8, frameon=False,
           title="Pathogen species", title_fontsize=8.5,
           bbox_to_anchor=(0.74, 0.52))

out = OUT_DIR / "masked_vs_unmasked.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out}")

# ── caption ───────────────────────────────────────────────────────────────────
caption = """\
## Figure X — Masked vs unmasked Kraken2 database: host and pathogen read assignment

![Masked vs unmasked Kraken2 DB comparison](masked_vs_unmasked.png)

### Background and motivation

NCBI STAT pre-computed k-mer taxonomy was used in the upstream pipeline to screen
~593k SRA runs for co-infection signal. Screening microbe-as-library (MAL) runs —
runs where the sequenced organism is itself a PHI-base plant pathogen, and which
should by definition contain pathogen reads — revealed that a large proportion
returned zero detected eukaryotic pathogen reads in STAT. Investigation of the most
affected group, wheat stripe rust (*Puccinia striiformis* f. sp. *tritici*, PST),
confirmed the blind spot: euk_pct = 0% across all PST pilot runs, while Kraken2
and kallisto independently detected *P. striiformis* at 10–68% of reads. STAT's
reference k-mer database inadequately covers Basidiomycota plant pathogens, making
it an unreliable screen for some of the most agronomically important fungal diseases.

To replace STAT, a custom Kraken2 database was built from the coding sequences (CDS)
of all PHI-base eukaryotic pathogens (fungi, oomycetes, nematodes; 205 seed taxa,
~15,000 taxids after taxonomic expansion) plus Viridiplantae host CDS for the 180
PHI-base host species. Two versions were evaluated:

- **Unmasked** — reference sequences used as-is.
- **Masked** — bidirectional BBDuk k-mer masking (`k=31`, `kmask=N`): pathogen
  sequences masked against all host CDS and against k-mers shared with any other
  pathogen in the database; host sequences masked against all pathogen CDS. Only
  species-diagnostic k-mers are retained.

### What this figure shows

Four-panel scatter plot comparing the two databases across **41 high-confidence
field RNA-seq runs** (biosample-representative, same-genus secondary pathogens
excluded, streamed from ENA FTP at 500,000 reads per run). Each point is one run;
the dashed line is *y* = *x*; points above the diagonal indicate greater assignment
by the unmasked database. Points are coloured by the dominant host (panels A–B) or
pathogen (panels C–D) species detected in the masked run.

**Panels A–B (Host species).**
Reads classified to host plant taxa and count of distinct host species detected.
Even the masked database detects a mean of 4.4 host species per run (unmasked: 11.0)
— biologically impossible given that each run comes from a single host species.
This noise arises from k-mer similarity between related plant genomes (e.g.
*Triticum aestivum*, *T. turgidum*, *Hordeum vulgare*, and *Brachypodium distachyon*
share large blocks of conserved sequence). Masking reduces but cannot eliminate this
problem because the shared sequence is intrinsic to plant genome evolution.

**Panels C–D (Pathogen species).**
Reads classified to PHI-base eukaryotic pathogens and count of pathogen species
detected. The masked database is more conservative (mean 18.9% vs 31.9% reads
assigned), but more specific: spurious cross-genus hits visible in the unmasked
results (e.g. *Melampsora laricis-populina* in wheat rust runs, *Sorghum bicolor*
in maize runs) are eliminated by masking. Species counts remain well-correlated
between databases (*r* = 0.95), confirming that masking does not suppress genuine
detections. Runs showing extreme displacement above the diagonal in panel C have
high k-mer ambiguity between closely related pathogen species (notably *Puccinia*
spp.) and produce unreliable estimates from either database.

### Key finding

Including host sequences in the database is fundamentally counterproductive.
Host genomes are too similar to each other (and to pathogen sequences) for k-mer
classification to reliably identify a single host species, and the noise they
introduce inflates species counts and competes with pathogen k-mer assignment.
Host identity is more accurately and efficiently obtained from SRA run metadata
(submitter-declared organism), which is available for all runs.

### Path forward

1. **Rebuild the database with pathogen CDS only** — no host sequences.
   Masking becomes pathogen-vs-pathogen: each species is masked against k-mers
   shared with any other pathogen in the database, leaving only species-diagnostic
   signal. The resulting database will be smaller (faster load, no Lustre mmap
   issues on Setonix), and pct_classified will represent pathogen burden directly.

2. **Infer host species from metadata** — use the `host` / `named_host` columns
   already present in `output/02_filter_runs/data/runs.tsv`.

3. **Filter unreliable runs** — exclude runs where unmasked pathogen % greatly
   exceeds masked pathogen % (threshold TBD after pathogen-only DB results are in).
   These runs have insufficient species-diagnostic signal to support confident
   pathogen identification.
"""

# remove old caption files if they exist
for old_name in ("masked_vs_unmasked_caption.txt", "masked_vs_unmasked_caption.md"):
    old = OUT_DIR / old_name
    if old.exists():
        old.unlink()

cap_path = OUT_DIR / "README.md"
cap_path.write_text(caption)
print(f"Saved {cap_path}")

# ── TSV export ────────────────────────────────────────────────────────────────
tsv_cols = [
    "run",
    "top_host", "top_pathogen",
    "masked_host_pct",    "unmasked_host_pct",
    "masked_host_sp",     "unmasked_host_sp",
    "masked_path_pct",    "unmasked_path_pct",
    "masked_path_sp",     "unmasked_path_sp",
    "masked_total_pct",   "unmasked_total_pct",
]
tsv_path = OUT_DIR / "masked_vs_unmasked.tsv"
with open(tsv_path, "w") as f:
    f.write("\t".join(tsv_cols) + "\n")
    for run in runs:
        r = rows[run]
        d = kraken[run]
        f.write("\t".join(str(x) for x in [
            run,
            r["top_host"],          r["top_path"],
            round(r["m_host_pct"],  3), round(r["u_host_pct"],  3),
            r["m_host_sp"],             r["u_host_sp"],
            round(r["m_path_pct"],  3), round(r["u_path_pct"],  3),
            r["m_path_sp"],             r["u_path_sp"],
            round(d["masked_pct"],  3), round(d["unmasked_pct"], 3),
        ]) + "\n")
print(f"Saved {tsv_path}")
