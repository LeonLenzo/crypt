#!/usr/bin/env python3
"""
scripts/review_excluded_hosts.py — review MAL runs excluded from the guild network
due to unresolved / broad-clade host assignment in STAT.

Outputs:
  output/04_filter_meta/data/excluded_hosts.tsv   one row per BioProject
  output/04_filter_meta/excluded_hosts.html        browsable HTML tool

Run from crypt/:
  python scripts/review_excluded_hosts.py
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

CRYPT_TSV  = Path("output/02_filter_runs/data/crypt.tsv")
META_TSV   = Path("output/04_filter_meta/data/bioproject_meta.tsv")
OUT_DIR    = Path("output/04_filter_meta")
OUT_TSV    = OUT_DIR / "data" / "excluded_hosts.tsv"
OUT_HTML   = OUT_DIR / "excluded_hosts.html"

BROAD_CLADE_NAMES = {
    "Viridiplantae", "Mesangiospermae", "eudicotyledons", "rosids",
    "Euphyllophyta", "Pentapetalae", "Poales", "asterids", "Gunneridae",
    "Magnoliopsida", "Lamiales", "BOP clade", "PACMAD clade",
    "Streptophyta",
}


def _host_category(h: str) -> str:
    if h == "Viridiplantae":
        return "unresolved"
    if h in BROAD_CLADE_NAMES:
        return "broad_clade"
    parts = h.split()
    if len(parts) >= 2 and parts[1][0].islower():
        return "species"
    return "tribe_genus"


def main() -> None:
    # ── Load BioProject metadata ──────────────────────────────────────────────
    bp_meta: dict[str, dict] = {}
    with open(META_TSV, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bp_meta[r["BioProject"]] = r

    field_bps = {bp for bp, m in bp_meta.items()
                 if m.get("llm_study_setting") == "field"}

    # ── Find excluded runs: MAL, biosample_rep, field, broad/unresolved host ──
    # Also collect total field MAL biosample_rep runs for context.
    bp_excluded: dict[str, Counter] = defaultdict(Counter)  # bp → host → count
    bp_primary:  dict[str, Counter] = defaultdict(Counter)  # bp → primary_pathogen → count
    bp_total_runs: dict[str, int] = defaultdict(int)

    with open(CRYPT_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("mode") != "mal":
                continue
            if row.get("biosample_representative") != "True":
                continue
            if row.get("BioProject") not in field_bps:
                continue
            bp = row["BioProject"]
            bp_total_runs[bp] += 1
            host = row.get("host", "").strip()
            cat  = _host_category(host)
            if cat in ("unresolved", "broad_clade"):
                bp_excluded[bp][host] += 1
                prim = row.get("primary_pathogen", "").strip()
                if prim:
                    bp_primary[bp][prim] += 1

    print(f"BioProjects with excluded runs: {len(bp_excluded)}")
    total_excl = sum(sum(c.values()) for c in bp_excluded.values())
    print(f"Total excluded runs: {total_excl}")

    # ── Write TSV ─────────────────────────────────────────────────────────────
    rows: list[dict] = []
    for bp in sorted(bp_excluded, key=lambda b: -sum(bp_excluded[b].values())):
        excl_counts = bp_excluded[bp]
        n_excl = sum(excl_counts.values())
        n_tot  = bp_total_runs.get(bp, 0)
        host_dist = "; ".join(f"{h}: {n}" for h, n in excl_counts.most_common())
        top_primary = bp_primary[bp].most_common(1)[0][0] if bp_primary[bp] else ""
        meta = bp_meta.get(bp, {})
        rows.append({
            "BioProject":       bp,
            "n_excluded":       n_excl,
            "n_total_mal":      n_tot,
            "pct_excluded":     f"{n_excl/max(n_tot,1)*100:.0f}",
            "host_distribution": host_dist,
            "top_primary":      top_primary,
            "llm_study_setting": meta.get("llm_study_setting", ""),
            "llm_treatment":    meta.get("llm_treatment", ""),
            "llm_named_pathogen": meta.get("llm_named_pathogen", ""),
            "title":            meta.get("title", ""),
            "primary_pmid":     meta.get("primary_pmid", ""),
        })

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Written: {OUT_TSV}")

    # ── Write HTML ────────────────────────────────────────────────────────────
    rows_js = json.dumps(rows, ensure_ascii=False)
    bp_meta_js = json.dumps(
        {bp: {k: bp_meta[bp].get(k, "") for k in
              ("title", "description", "abstract", "methods_text",
               "llm_treatment", "llm_study_setting", "llm_named_pathogen",
               "llm_rationale", "llm_confidence", "primary_pmid",
               "primary_publication")}
         for bp in bp_excluded},
        ensure_ascii=False
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Excluded host runs — crypt pipeline</title>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 0; background: #f8f9fa; color: #222; }}
  h1   {{ font-size: 1.2rem; padding: 12px 16px; margin: 0;
          background: #2c3e50; color: #fff; }}
  #info {{ padding: 6px 16px; background: #ecf0f1; font-size: 0.85rem; color: #555; }}
  #controls {{ padding: 8px 16px; background: #fff; border-bottom: 1px solid #ddd;
               display: flex; gap: 8px; align-items: center; flex-wrap: wrap; }}
  #controls input {{ padding: 4px 8px; border: 1px solid #ccc; border-radius: 4px;
                      font-size: 0.85rem; width: 280px; }}
  #controls label {{ font-size: 0.83rem; color: #555; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.83rem; }}
  th    {{ background: #34495e; color: #fff; padding: 6px 8px; text-align: left;
           cursor: pointer; user-select: none; position: sticky; top: 0; }}
  th:hover {{ background: #2c3e50; }}
  td    {{ padding: 5px 8px; border-bottom: 1px solid #eee; vertical-align: top; }}
  tr:hover td {{ background: #fef9e7; cursor: pointer; }}
  .n-excl {{ font-weight: bold; color: #c0392b; }}
  .pct    {{ color: #7f8c8d; font-size: 0.78rem; }}
  .host-dist {{ max-width: 340px; color: #555; }}
  .title  {{ max-width: 320px; }}
  .pmid a {{ color: #2980b9; text-decoration: none; font-size: 0.78rem; }}

  /* modal */
  #overlay {{ display:none; position:fixed; inset:0; background:rgba(0,0,0,.45);
              z-index:100; }}
  #modal   {{ position:fixed; top:5vh; left:50%; transform:translateX(-50%);
              width:min(900px,94vw); max-height:88vh; overflow-y:auto;
              background:#fff; border-radius:8px; padding:24px; z-index:101;
              box-shadow:0 8px 32px rgba(0,0,0,.3); }}
  #modal h2 {{ margin:0 0 4px; font-size:1.05rem; }}
  #modal .bp-id {{ font-size:0.82rem; color:#888; margin-bottom:12px; }}
  .section {{ margin-bottom:14px; }}
  .section h3 {{ font-size:0.85rem; font-weight:600; color:#555;
                 border-bottom:1px solid #eee; padding-bottom:3px;
                 margin:0 0 6px; }}
  .section p, .section pre {{ font-size:0.82rem; margin:0; white-space:pre-wrap;
                               word-break:break-word; color:#333; }}
  .tag {{ display:inline-block; padding:2px 7px; border-radius:10px;
          font-size:0.75rem; font-weight:600; margin:0 3px 3px 0; }}
  .tag-treatment {{ background:#dceefb; color:#1a6fa0; }}
  .tag-setting   {{ background:#d5f5e3; color:#1a7a42; }}
  .tag-conf      {{ background:#fef3cd; color:#856404; }}
  .host-row {{ display:flex; gap:6px; flex-wrap:wrap; margin-top:6px; }}
  .host-chip {{ background:#fdebd0; color:#784212; border-radius:4px;
                padding:2px 8px; font-size:0.78rem; }}
  #close-btn {{ float:right; background:none; border:none; font-size:1.3rem;
                cursor:pointer; color:#888; padding:0; }}
  #close-btn:hover {{ color:#333; }}
</style>
</head>
<body>
<h1>Excluded host runs — MAL field biosample_representative</h1>
<div id="info">Runs excluded from the guild network due to unresolved or broad-clade host assignment in NCBI STAT.<br>
Use publication metadata to identify the likely host plant and inform future LLM host extraction.</div>
<div id="controls">
  <input id="search" placeholder="Search title, BioProject, pathogen, host…" oninput="applyFilters()">
  <label>Min excluded runs: <input id="min-excl" type="number" value="1" min="1" style="width:60px"
    oninput="applyFilters()"></label>
  <span id="count-display" style="font-size:0.82rem;color:#555;margin-left:8px;"></span>
</div>
<div style="overflow-x:auto">
<table id="bp-table">
  <thead>
    <tr>
      <th onclick="sortBy('BioProject')">BioProject</th>
      <th onclick="sortBy('n_excluded')">Excl. runs</th>
      <th onclick="sortBy('pct_excluded')">% of total</th>
      <th onclick="sortBy('top_primary')">Top primary pathogen</th>
      <th>Host distribution</th>
      <th onclick="sortBy('title')">Title</th>
      <th>PMID</th>
    </tr>
  </thead>
  <tbody id="tbody"></tbody>
</table>
</div>

<div id="overlay" onclick="closeModal()"></div>
<div id="modal" style="display:none">
  <button id="close-btn" onclick="closeModal()">&#x2715;</button>
  <h2 id="m-title"></h2>
  <div class="bp-id" id="m-bp"></div>
  <div class="section">
    <h3>LLM classification</h3>
    <div id="m-tags"></div>
    <p id="m-rationale"></p>
  </div>
  <div class="section">
    <h3>Excluded host STAT values</h3>
    <div class="host-row" id="m-hosts"></div>
  </div>
  <div class="section" id="m-abstract-sec">
    <h3>Abstract</h3>
    <p id="m-abstract"></p>
  </div>
  <div class="section" id="m-methods-sec">
    <h3>Methods (PMC extract)</h3>
    <pre id="m-methods"></pre>
  </div>
  <div class="section" id="m-desc-sec">
    <h3>BioProject description</h3>
    <p id="m-desc"></p>
  </div>
</div>

<script>
const ROWS    = {rows_js};
const BP_META = {bp_meta_js};

let sortKey = 'n_excluded', sortAsc = false;

function sortBy(key) {{
  if (sortKey === key) sortAsc = !sortAsc;
  else {{ sortKey = key; sortAsc = false; }}
  applyFilters();
}}

function applyFilters() {{
  const q       = document.getElementById('search').value.toLowerCase();
  const minExcl = parseInt(document.getElementById('min-excl').value) || 1;
  let visible = ROWS.filter(r =>
    r.n_excluded >= minExcl &&
    (!q || r.BioProject.toLowerCase().includes(q) ||
           r.title.toLowerCase().includes(q) ||
           r.top_primary.toLowerCase().includes(q) ||
           r.host_distribution.toLowerCase().includes(q))
  );
  visible.sort((a, b) => {{
    let av = a[sortKey], bv = b[sortKey];
    if (!isNaN(av) && !isNaN(bv)) {{ av = +av; bv = +bv; }}
    if (av < bv) return sortAsc ? -1 : 1;
    if (av > bv) return sortAsc ? 1 : -1;
    return 0;
  }});
  const tbody = document.getElementById('tbody');
  tbody.innerHTML = visible.map(r => `
    <tr onclick="openModal('${{r.BioProject}}')">
      <td><code>${{r.BioProject}}</code></td>
      <td class="n-excl">${{r.n_excluded}}</td>
      <td class="pct">${{r.pct_excluded}}%</td>
      <td>${{r.top_primary}}</td>
      <td class="host-dist">${{r.host_distribution}}</td>
      <td class="title">${{r.title}}</td>
      <td class="pmid">${{r.primary_pmid
        ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${{r.primary_pmid}}/" target="_blank">${{r.primary_pmid}}</a>`
        : '—'}}</td>
    </tr>`).join('');
  document.getElementById('count-display').textContent =
    `${{visible.length}} BioProject${{visible.length !== 1 ? 's' : ''}}`;
}}

function openModal(bp) {{
  const row  = ROWS.find(r => r.BioProject === bp);
  const meta = BP_META[bp] || {{}};
  document.getElementById('m-title').textContent = meta.title || bp;
  document.getElementById('m-bp').textContent    = bp +
    (meta.primary_pmid ? '  ·  PMID ' + meta.primary_pmid : '');
  // tags
  const tags = [
    meta.llm_treatment    ? `<span class="tag tag-treatment">${{meta.llm_treatment}}</span>` : '',
    meta.llm_study_setting? `<span class="tag tag-setting">${{meta.llm_study_setting}}</span>` : '',
    meta.llm_confidence   ? `<span class="tag tag-conf">conf: ${{meta.llm_confidence}}</span>` : '',
  ].join('');
  document.getElementById('m-tags').innerHTML = tags;
  document.getElementById('m-rationale').textContent = meta.llm_rationale || '';
  // hosts
  const hostHtml = (row.host_distribution || '').split('; ').map(s =>
    `<span class="host-chip">${{s}}</span>`).join('');
  document.getElementById('m-hosts').innerHTML = hostHtml;
  // text sections
  const abstract = meta.abstract || '';
  document.getElementById('m-abstract').textContent = abstract;
  document.getElementById('m-abstract-sec').style.display = abstract ? '' : 'none';
  const methods = meta.methods_text || '';
  document.getElementById('m-methods').textContent = methods;
  document.getElementById('m-methods-sec').style.display = methods ? '' : 'none';
  const desc = meta.description || '';
  document.getElementById('m-desc').textContent = desc;
  document.getElementById('m-desc-sec').style.display = desc ? '' : 'none';

  document.getElementById('overlay').style.display = '';
  document.getElementById('modal').style.display   = '';
}}

function closeModal() {{
  document.getElementById('overlay').style.display = 'none';
  document.getElementById('modal').style.display   = 'none';
}}

document.addEventListener('keydown', e => {{ if (e.key === 'Escape') closeModal(); }});
applyFilters();
</script>
</body>
</html>"""

    with open(OUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Written: {OUT_HTML}")


if __name__ == "__main__":
    main()
