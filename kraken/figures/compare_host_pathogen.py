"""
DB comparison figures: masked+hosts vs unmasked, masked+hosts vs pathogens-only,
and STAT euk% vs pathogens-only.

Outputs:
  masked_vs_unmasked.png   – 2×2 panel (motivation for redesign)
  masked_vs_pathogens.png  – 2×2 panel (result of redesign)
  stat_vs_pathogens.png    – 1-panel scatter (STAT vs pathogens-only DB)
  masked_vs_unmasked.tsv   – full data table

Inputs:
  /tmp/setonix_kraken_metrics.json        – merged masked + unmasked + pathogens-only
  output/00_build/data/phibase_db.json
  output/02_filter_runs/data/runs.tsv     – STAT euk% per run

Run from crypt/:
  python output/figures/db_comparison/compare_host_pathogen.py
"""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── paths ─────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
KRAKEN_JSON = Path("/tmp/setonix_kraken_metrics.json")
PHIBASE     = ROOT / "stat/output/data/phibase_db.json"
OUT_DIR     = Path("kraken/output/figures/db_comparison")

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

# ── Figure 3: STAT euk% vs pathogens-only ────────────────────────────────────
RUNS_TSV = ROOT / "stat/output/data/runs.tsv"
stat_data = {}
with open(RUNS_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["Run"] in runs_all:
            sp_str = row.get("stat_pathogens", "")
            # parse "Name:pct%; Name2:pct2%" into {name: pct}
            stat_sp = {}
            for entry in sp_str.split(";"):
                entry = entry.strip()
                if ":" in entry:
                    name, pct = entry.rsplit(":", 1)
                    try:
                        stat_sp[name.strip()] = float(pct.strip().rstrip("%"))
                    except ValueError:
                        pass
            stat_data[row["Run"]] = {
                "euk_pct": float(row.get("fungi_pct", 0) or 0) + float(row.get("oomycete_pct", 0) or 0),
                "stat_sp": stat_sp,
                "stat_top": max(stat_sp, key=stat_sp.get) if stat_sp else "",
            }

stat_runs = [r for r in runs_all if r in stat_data]

# Does STAT's top species match the DB's top pathogen (genus level)?
def genus(name):
    return name.split()[0] if name else ""

def stat_matches_db(run):
    db_top   = rows_all[run]["top_path_po"]
    stat_top = stat_data[run]["stat_top"]
    return genus(db_top) == genus(stat_top)

sx      = np.array([stat_data[r]["euk_pct"]           for r in stat_runs])
sy      = np.array([kraken[r]["pathogens_pct"]         for r in stat_runs])
s_match = [stat_matches_db(r) for r in stat_runs]
spc     = [PATHOGEN_COLOURS.get(rows_all[r]["top_path_po"], DEFAULT) for r in stat_runs]

n_match    = sum(s_match)
n_mismatch = len(stat_runs) - n_match
r_val = float(np.corrcoef(sx, sy)[0, 1]) if np.std(sx) > 0 and np.std(sy) > 0 else 0.0

# legend patches for top DB pathogens
seen_sp = {}
for r in stat_runs:
    name = rows_all[r]["top_path_po"]
    seen_sp.setdefault(PATHOGEN_COLOURS.get(name, DEFAULT), name)
sp_patches = [mpatches.Patch(color=c, label=lbl) for c, lbl in seen_sp.items()]

# raw STAT species counts (full specific_hits, not post-filter column)
LEAF_FRAC = 0.9
def _node_count(table, node):
    for e in table:
        if e.get("org") == node:
            return e.get("total_count", 0)
    return 0

def _specific_hits(table, node, analyzed, min_pct=0.5):
    top = _node_count(table, node)
    if not top:
        return []
    min_count = max(1, int(analyzed * min_pct / 100))
    under = [(e["org"], e["total_count"]) for e in table
             if 0 < e.get("total_count", 0) <= top
             and e.get("org") != node
             and e.get("total_count", 0) >= min_count]
    if not under:
        return [(node, top / analyzed * 100)]
    count_vals  = sorted({c for _, c in under}, reverse=True)
    leaf_counts = {c for c in count_vals
                   if not any(c2 < c and c2 >= LEAF_FRAC * c for c2 in count_vals)}
    best = {}
    for name, count in under:
        if count not in leaf_counts:
            continue
        pct = count / analyzed * 100
        if count not in best or (len(name.split()) == 2 and len(best[count][0].split()) != 2):
            best[count] = (name, pct)
    return sorted(best.values(), key=lambda x: -x[1])

stat_raw = {}
with open(ROOT / "stat/output/data/stat_cache.jsonl") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        run, _, rest = line.partition("\t")
        if run in runs_all:
            stat_raw[run] = json.loads(rest)

stat_sp_raw = {}
for run, entries in stat_raw.items():
    entry    = entries[0]
    table    = entry["tax_table"]
    totals   = entry["tax_totals"]
    analyzed = totals.get("analysed") or totals.get("analyzed") or totals.get("total", 1)
    hits = []
    for kingdom, thresh in [("Fungi", 0.5), ("Oomycota", 0.5)]:
        hits.extend(_specific_hits(table, kingdom, analyzed, min_pct=thresh))
    stat_sp_raw[run] = {name: pct for name, pct in hits}

# species count arrays (only runs with stat_raw)
sp_stat_runs = [r for r in stat_runs if r in stat_sp_raw]
sn_stat = np.array([len(stat_sp_raw[r]) for r in sp_stat_runs], dtype=float)
sn_db   = np.array([len(kraken[r].get("pathogens_species", {})) for r in sp_stat_runs], dtype=float)
sn_clrs = [PATHOGEN_COLOURS.get(rows_all[r]["top_path_po"], DEFAULT) for r in sp_stat_runs]

lim3 = 70
lim3b = max(sn_stat.max(), sn_db.max()) * 1.1

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(13.0, 5.5))
fig3.subplots_adjust(left=0.08, right=0.73, top=0.91, bottom=0.12, wspace=0.32)

# Panel A: reads %
ax3a.plot([0, lim3], [0, lim3], "--", color="#999999", lw=0.9, zorder=1)
match_runs    = [r for r, m in zip(stat_runs, s_match) if m]
mismatch_runs = [r for r, m in zip(stat_runs, s_match) if not m]
ax3a.scatter(
    [stat_data[r]["euk_pct"] for r in match_runs],
    [kraken[r]["pathogens_pct"] for r in match_runs],
    c=[PATHOGEN_COLOURS.get(rows_all[r]["top_path_po"], DEFAULT) for r in match_runs],
    s=50, alpha=0.88, edgecolors="white", linewidths=0.4, zorder=2,
    label=f"Top genus matches (n={n_match})"
)
ax3a.scatter(
    [stat_data[r]["euk_pct"] for r in mismatch_runs],
    [kraken[r]["pathogens_pct"] for r in mismatch_runs],
    c=[PATHOGEN_COLOURS.get(rows_all[r]["top_path_po"], DEFAULT) for r in mismatch_runs],
    s=70, alpha=0.88, edgecolors="#333333", linewidths=1.2, marker="D", zorder=3,
    label=f"Top genus mismatch (n={n_mismatch})"
)
ax3a.text(0.97, 0.04, f"r = {r_val:.2f}", transform=ax3a.transAxes,
          ha="right", va="bottom", fontsize=8, color="#666666")
ax3a.set_xlim(0, lim3); ax3a.set_ylim(0, lim3)
ax3a.set_xlabel("STAT eukaryotic pathogen (%)", fontsize=9)
ax3a.set_ylabel("Pathogens-only DB (%)", fontsize=9)
ax3a.set_title("A   Reads assigned (%)", fontsize=10, fontweight="bold", pad=5, loc="left")
ax3a.tick_params(labelsize=8)
ax3a.legend(loc="lower right", fontsize=7.5, frameon=False, handletextpad=0.4)

# Panel B: species count
r_sp = float(np.corrcoef(sn_stat, sn_db)[0, 1]) if np.std(sn_stat) > 0 else 0.0
ax3b.plot([0, lim3b], [0, lim3b], "--", color="#999999", lw=0.9, zorder=1)
ax3b.scatter(sn_stat, sn_db, c=sn_clrs, s=50, alpha=0.88,
             edgecolors="white", linewidths=0.4, zorder=2)
ax3b.text(0.97, 0.04, f"r = {r_sp:.2f}", transform=ax3b.transAxes,
          ha="right", va="bottom", fontsize=8, color="#666666")
ax3b.set_xlim(0, lim3b); ax3b.set_ylim(0, lim3b)
ax3b.set_xlabel("STAT species detected", fontsize=9)
ax3b.set_ylabel("Pathogens-only DB species detected", fontsize=9)
ax3b.set_title("B   Species detected", fontsize=10, fontweight="bold", pad=5, loc="left")
ax3b.tick_params(labelsize=8)

fig3.legend(handles=sp_patches, loc="upper left", fontsize=8, frameon=False,
            title="Top pathogen (DB)", title_fontsize=8.5,
            bbox_to_anchor=(0.745, 0.95))

print(f"\nstat_vs_pathogens.png  (n={len(stat_runs)})")
print(f"  STAT euk%:        mean={sx.mean():.1f}%  median={np.median(sx):.1f}%")
print(f"  DB path%:         mean={sy.mean():.1f}%  median={np.median(sy):.1f}%")
print(f"  Genus match:      {n_match}/{len(stat_runs)}")
print(f"  STAT species/run: mean={sn_stat.mean():.1f}  range={int(sn_stat.min())}-{int(sn_stat.max())}")
print(f"  DB species/run:   mean={sn_db.mean():.1f}    range={int(sn_db.min())}-{int(sn_db.max())}")

out3 = OUT_DIR / "stat_vs_pathogens.png"
fig3.savefig(out3, dpi=180, bbox_inches="tight")
plt.close(fig3)
print(f"  Saved {out3}")

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
## Kraken2 database comparison: masked+hosts vs unmasked vs pathogens-only vs STAT

### Figure 1: Masked+hosts vs unmasked (`masked_vs_unmasked.png`)

![Masked vs unmasked Kraken2 DB comparison](masked_vs_unmasked.png)

### Figure 2: Masked+hosts vs pathogens-only (`masked_vs_pathogens.png`)

![Masked+hosts vs pathogens-only](masked_vs_pathogens.png)

### Figure 3: STAT euk% vs pathogens-only Kraken2 DB (`stat_vs_pathogens.png`)

![STAT vs pathogens-only DB](stat_vs_pathogens.png)

### Background and motivation

NCBI STAT pre-computed k-mer taxonomy was used in the upstream pipeline to screen
~593k SRA runs for co-infection signal. Screening microbe-as-library (MAL) runs —
runs where the sequenced organism is itself a PHI-base plant pathogen — revealed
that a large proportion returned zero detected eukaryotic pathogen reads in STAT,
particularly for wheat stripe rust (*Puccinia striiformis* f. sp. *tritici*, PST),
while Kraken2 and kallisto independently detected *P. striiformis* at 10–68% of
reads. STAT's reference k-mer database inadequately covers Basidiomycota plant
pathogens.

To replace STAT, a custom Kraken2 database was built from the CDS of PHI-base
eukaryotic pathogens (fungi + oomycetes). Three versions were evaluated across 32–41
high-confidence field RNA-seq runs (biosample-representative, same-genus secondary
pathogens excluded):

- **Unmasked** — pathogen + host CDS used as-is.
- **Masked** — bidirectional BBDuk masking (`k=35`): pathogens masked against shared
  pathogen k-mers + all host CDS; hosts masked against all pathogen CDS.
- **Pathogens-only** — pathogen CDS only, masked against k-mers shared between any
  two pathogens. No host sequences.

### Figure 1: Why host sequences were removed (n=41)

**Panels A–B (Host species).** Even the masked database detects a mean of 4.4 host
species per run — biologically impossible (one host per run). Masking cannot
eliminate this noise because k-mer similarity between related plant genomes is
intrinsic to plant genome evolution.

**Panels C–D (Pathogen species).** Masked is more conservative (mean 18.9% vs
31.9%) and more specific: spurious cross-genus hits in the unmasked results are
eliminated. Species counts correlate well (*r* = 0.95), confirming masking does not
suppress genuine detections.

### Figure 2: Effect of removing host sequences (n=32)

**Panels A–B (Host species).** Pathogens-only assigns zero reads to host species
(0.0% vs 4.5% for masked). Host noise is completely resolved.

**Panels C–D (Pathogen species).** Pathogen detection is consistent between
masked+hosts and pathogens-only (mean 17.9% vs 18.3%; same dominant species per
run). Pathogens-only recovers slightly more species (9.6 vs 8.2) as k-mers no
longer compete with host sequences.

### Figure 3: STAT vs pathogens-only Kraken2 DB (n=32)

Points above the diagonal indicate the Kraken2 DB detects more than STAT; points
below indicate STAT reports more. Two distinct patterns emerge by pathogen group:

- **Basidiomycota rusts (*Puccinia* spp.):** STAT strongly underestimates (STAT
  2–9% vs DB 10–41% for PST runs). STAT's k-mer reference inadequately covers
  rust fungi, confirming the original motivation for this pipeline.
- **Ascomycota (*Zymoseptoria tritici*, *Puccinia graminis*):** STAT tends to
  overestimate vs the DB (STAT 37–41% vs DB 7–13% for some *Z. tritici* runs),
  likely due to broad k-mer matches to non-diagnostic sequences in the STAT
  reference.
- **Overall correlation:** *r* = {r_val:.2f} across all 32 runs, with the divergence
  explained almost entirely by pathogen taxonomy. Runs agree well for *Sclerotinia*
  and *Monilinia*.

### Key finding

The pathogens-only Kraken2 DB outperforms STAT for Basidiomycota pathogens,
eliminates host noise entirely, and is more specific than the unmasked DB for
Ascomycota pathogens. It is the recommended tool for full production screening
of ~8,243 HC biosample-representative runs.
""".format(r_val=float(np.corrcoef(sx, sy)[0, 1]))

for old in (OUT_DIR / "masked_vs_unmasked_caption.txt",
            OUT_DIR / "masked_vs_unmasked_caption.md"):
    if old.exists():
        old.unlink()

(OUT_DIR / "README.md").write_text(readme)
print(f"Saved {OUT_DIR / 'README.md'}")
