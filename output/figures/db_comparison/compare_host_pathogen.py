"""
DB comparison figures: masked+hosts vs unmasked, and masked+hosts vs pathogens-only.

Outputs:
  masked_vs_unmasked.png   – 2×2 panel (existing; motivation for redesign)
  masked_vs_pathogens.png  – 2×2 panel (new; result of redesign)
  masked_vs_unmasked.tsv   – full data table

Inputs:
  /tmp/setonix_kraken_metrics.json   – merged masked + unmasked + pathogens-only results
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

# All runs with masked+unmasked results
runs_mu = sorted(kraken)
# Runs with all three DB results
runs_all = sorted(r for r in kraken if kraken[r].get("pathogens_pct") is not None)

print(f"masked vs unmasked:      n = {len(runs_mu)}")
print(f"masked vs pathogens-only: n = {len(runs_all)}")

def summarise(species_dict, kind):
    pct, n = 0.0, 0
    for name, p in species_dict.items():
        if classify(name) == kind:
            pct += p
            n   += 1
    return pct, n

def build_rows(run_list):
    rows = {}
    for run in run_list:
        d = kraken[run]
        m  = d["masked_species"]
        u  = d["unmasked_species"]
        p  = d.get("pathogens_species", {})
        rows[run] = {
            "m_host_pct":  summarise(m, "host")[0],
            "u_host_pct":  summarise(u, "host")[0],
            "p_host_pct":  summarise(p, "host")[0],
            "m_host_sp":   summarise(m, "host")[1],
            "u_host_sp":   summarise(u, "host")[1],
            "p_host_sp":   summarise(p, "host")[1],
            "m_path_pct":  summarise(m, "pathogen")[0],
            "u_path_pct":  summarise(u, "pathogen")[0],
            "p_path_pct":  summarise(p, "pathogen")[0],
            "m_path_sp":   summarise(m, "pathogen")[1],
            "u_path_sp":   summarise(u, "pathogen")[1],
            "p_path_sp":   summarise(p, "pathogen")[1],
            "top_host":    next((n for n in m if classify(n) == "host"),     "Unknown"),
            "top_path":    next((n for n in m if classify(n) == "pathogen"), "Unknown"),
            "top_path_po": next((n for n in p if classify(n) == "pathogen"), "Unknown"),
        }
    return rows

rows_mu  = build_rows(runs_mu)
rows_all = build_rows(runs_all)

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

def colours_patches(run_list, rows, key, cmap):
    clrs = [cmap.get(rows[r][key], DEFAULT) for r in run_list]
    seen = {}
    for r in run_list:
        seen.setdefault(cmap.get(rows[r][key], DEFAULT), rows[r][key])
    ptchs = [mpatches.Patch(color=c, label=lbl) for c, lbl in seen.items()]
    return clrs, ptchs

# ── scatter helper ─────────────────────────────────────────────────────────────
def panel(ax, x, y, clrs, xlabel, ylabel, lim=100):
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

def make_2x2(runs, rows, y_label_short, y_label_long, hc, hp, pc, pp, outpath):
    mhp = np.array([rows[r]["m_host_pct"]  for r in runs])
    yhp = np.array([rows[r][f"{y_label_short}_host_pct"]  for r in runs])
    mhs = np.array([rows[r]["m_host_sp"]   for r in runs], dtype=float)
    yhs = np.array([rows[r][f"{y_label_short}_host_sp"]   for r in runs], dtype=float)
    mpp = np.array([rows[r]["m_path_pct"]  for r in runs])
    ypp = np.array([rows[r][f"{y_label_short}_path_pct"]  for r in runs])
    mps = np.array([rows[r]["m_path_sp"]   for r in runs], dtype=float)
    yps = np.array([rows[r][f"{y_label_short}_path_sp"]   for r in runs], dtype=float)

    print(f"\n{outpath.name}  (n={len(runs)})")
    print(f"  Host reads:    masked {mhp.mean():.1f}%  {y_label_long} {yhp.mean():.1f}%")
    print(f"  Host species:  masked {mhs.mean():.1f}   {y_label_long} {yhs.mean():.1f}")
    print(f"  Path reads:    masked {mpp.mean():.1f}%  {y_label_long} {ypp.mean():.1f}%")
    print(f"  Path species:  masked {mps.mean():.1f}   {y_label_long} {yps.mean():.1f}")

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.5))
    fig.subplots_adjust(hspace=0.45, wspace=0.30, left=0.09, right=0.72,
                        top=0.93, bottom=0.08)

    clrs = [[hc, hc], [pc, pc]]
    xs   = [[mhp, mhs], [mpp, mps]]
    ys   = [[yhp, yhs], [ypp, yps]]
    ylabels_pct = f"{y_label_long} (%)"
    ylabels_sp  = f"{y_label_long} (# spp.)"

    for ri in range(2):
        for ci in range(2):
            ax = axes[ri, ci]
            panel(ax, xs[ri][ci], ys[ri][ci], clrs[ri][ci],
                  "Masked (%)" if ci == 0 else "Masked (# spp.)",
                  ylabels_pct  if ci == 0 else ylabels_sp)
            ax.text(0.03, 0.97, "ABCD"[ri*2+ci], transform=ax.transAxes,
                    fontsize=10, fontweight="bold", va="top")
            if ri == 0:
                ax.set_title(["Reads assigned (%)", "Species detected"][ci],
                             fontsize=10, fontweight="bold", pad=5)
            if ci == 0:
                ax.set_ylabel(ylabels_pct, fontsize=8.5)

    row_titles = ["Host species", "Pathogen species"]
    for ri, label in enumerate(row_titles):
        y0 = axes[ri, 0].get_position().y0
        y1 = axes[ri, 0].get_position().y1
        fig.text(0.01, (y0 + y1) / 2, label, va="center", ha="left",
                 rotation="vertical", fontsize=10, fontweight="bold", color="#222222")

    fig.legend(handles=hp, loc="upper left", fontsize=8, frameon=False,
               title="Host species", title_fontsize=8.5,
               bbox_to_anchor=(0.74, 0.95))
    fig.legend(handles=pp, loc="upper left", fontsize=8, frameon=False,
               title="Pathogen species", title_fontsize=8.5,
               bbox_to_anchor=(0.74, 0.52))

    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")

# ── Figure 1: masked vs unmasked ──────────────────────────────────────────────
hc_mu, hp_mu = colours_patches(runs_mu,  rows_mu,  "top_host", HOST_COLOURS)
pc_mu, pp_mu = colours_patches(runs_mu,  rows_mu,  "top_path", PATHOGEN_COLOURS)
make_2x2(runs_mu, rows_mu, "u", "Unmasked",
         hc_mu, hp_mu, pc_mu, pp_mu,
         OUT_DIR / "masked_vs_unmasked.png")

# ── Figure 2: masked vs pathogens-only ────────────────────────────────────────
hc_po, hp_po = colours_patches(runs_all, rows_all, "top_host",  HOST_COLOURS)
pc_po, pp_po = colours_patches(runs_all, rows_all, "top_path_po", PATHOGEN_COLOURS)
make_2x2(runs_all, rows_all, "p", "Pathogens-only",
         hc_po, hp_po, pc_po, pp_po,
         OUT_DIR / "masked_vs_pathogens.png")

# ── TSV export ────────────────────────────────────────────────────────────────
tsv_cols = [
    "run",
    "top_host", "top_pathogen",
    "masked_host_pct",    "unmasked_host_pct",    "pathogens_host_pct",
    "masked_host_sp",     "unmasked_host_sp",     "pathogens_host_sp",
    "masked_path_pct",    "unmasked_path_pct",    "pathogens_path_pct",
    "masked_path_sp",     "unmasked_path_sp",     "pathogens_path_sp",
    "masked_total_pct",   "unmasked_total_pct",   "pathogens_total_pct",
]
tsv_path = OUT_DIR / "masked_vs_unmasked.tsv"
with open(tsv_path, "w") as f:
    f.write("\t".join(tsv_cols) + "\n")
    for run in runs_mu:
        r = rows_mu[run]
        d = kraken[run]
        f.write("\t".join(str(x) for x in [
            run,
            r["top_host"], r["top_path"],
            round(r["m_host_pct"], 3), round(r["u_host_pct"], 3), round(r["p_host_pct"], 3),
            r["m_host_sp"],            r["u_host_sp"],            r["p_host_sp"],
            round(r["m_path_pct"], 3), round(r["u_path_pct"], 3), round(r["p_path_pct"], 3),
            r["m_path_sp"],            r["u_path_sp"],            r["p_path_sp"],
            round(d["masked_pct"],  3), round(d["unmasked_pct"], 3),
            round(d.get("pathogens_pct") or 0, 3),
        ]) + "\n")
print(f"\nSaved {tsv_path}")

# ── README ────────────────────────────────────────────────────────────────────
readme = """\
## Figure — Kraken2 database comparison: masked+hosts vs unmasked vs pathogens-only

### Figure 1: Masked+hosts vs unmasked (`masked_vs_unmasked.png`)

![Masked vs unmasked Kraken2 DB comparison](masked_vs_unmasked.png)

### Figure 2: Masked+hosts vs pathogens-only (`masked_vs_pathogens.png`)

![Masked+hosts vs pathogens-only](masked_vs_pathogens.png)

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

To replace STAT, a custom Kraken2 database was built from the CDS of all PHI-base
eukaryotic pathogens (fungi + oomycetes). Three versions were evaluated across
high-confidence field RNA-seq runs (biosample-representative, ≥1 HC pathogen):

- **Unmasked** — pathogen + host CDS used as-is.
- **Masked** — bidirectional BBDuk k-mer masking (`k=35`): pathogens masked against
  shared pathogen k-mers + all host CDS; hosts masked against all pathogen CDS.
- **Pathogens-only** — pathogen CDS only, masked against k-mers shared between
  any two pathogens (pathogen-vs-pathogen masking). No host sequences.

### Figure 1: Why host sequences were removed

**Panels A–B (Host species).**
Even the masked database detects a mean of 4.4 host species per run — biologically
impossible (one host per run). This noise arises from k-mer similarity between
related plant genomes. Masking reduces but cannot eliminate the problem.

**Panels C–D (Pathogen species).**
The masked database is more conservative (mean 18.9% vs 31.9% reads assigned) and
more specific: spurious cross-genus hits in the unmasked results (e.g.
*Melampsora laricis-populina* in wheat rust runs) are eliminated. Species counts
remain well-correlated (*r* = 0.95), confirming masking does not suppress genuine
detections.

### Figure 2: Effect of removing host sequences

**Panels A–B (Host species).**
Pathogens-only assigns effectively zero reads to host species (mean ~0% vs 4.4%
for masked). The host noise problem is completely resolved by excluding host CDS.

**Panels C–D (Pathogen species).**
Pathogen read assignment is broadly consistent between masked+hosts and
pathogens-only (the same dominant species are detected in each run). Pathogens-only
tends to assign slightly more reads to pathogens, as k-mers previously competing
with host sequences are now retained in the pathogen index.

### Key finding

Including host sequences in the database is fundamentally counterproductive.
Host identity is more accurately and efficiently obtained from SRA run metadata
(submitter-declared organism). The pathogens-only database is smaller (~1.9 GB vs
3.4 GB), loads faster on Lustre, and assigns 100% of classified reads to the
pathogen signal of interest.
"""

for old in (OUT_DIR / "masked_vs_unmasked_caption.txt",
            OUT_DIR / "masked_vs_unmasked_caption.md"):
    if old.exists():
        old.unlink()

(OUT_DIR / "README.md").write_text(readme)
print(f"Saved {OUT_DIR / 'README.md'}")
