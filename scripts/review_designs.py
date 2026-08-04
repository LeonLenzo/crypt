#!/usr/bin/env python3
"""
review_designs.py — generate design-review artefacts from bioproject_meta.tsv.

Outputs:
  output/04_filter_meta/data/subset_{design}.tsv   one per study_design category
  output/04_filter_meta/review_designs.html          interactive annotation tool

The HTML file is fully self-contained (no network required).  Open it in any
browser.  List view shows one compact row per BioProject; click any row to open
a detail modal with full abstract, methods text, and annotation controls.
Annotations are saved to localStorage as you work; Export downloads a TSV of
manual overrides.

Run from crypt/:
  python scripts/review_designs.py
"""

import csv
import json
from collections import Counter
from pathlib import Path

META_TSV   = Path("output/04_filter_meta/data/bioproject_meta.tsv")
CACHE_JSON = Path("output/03_fetch_meta/data/find_cache.json")
ALIGN_TSV  = Path("output/analysis/primary_alignment.tsv")
OUT_DIR    = Path("output/04_filter_meta")

TREATMENT_TYPES = ["coinf_experiment", "combined_stress", "abiotic_stress", "host_study", "single", "surveillance", "unclear"]
STUDY_SETTINGS  = ["field", "lab", "mixed", "unclear", "no_data"]

TREATMENT_COLOURS = {
    "coinf_experiment": "#c0392b",
    "combined_stress":  "#e74c3c",
    "abiotic_stress":   "#1abc9c",
    "host_study":       "#8e44ad",
    "single":           "#2980b9",
    "surveillance":     "#d35400",
    "unclear":          "#e67e22",
}
SETTING_COLOURS = {
    "field":   "#27ae60",
    "lab":     "#8e44ad",
    "mixed":   "#2980b9",
    "unclear": "#e67e22",
    "no_data": "#95a5a6",
}

# ── Subset TSVs ────────────────────────────────────────────────────────────────

def write_subsets(rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    axes = [
        ("treatment",     TREATMENT_TYPES, "treat"),
        ("study_setting", STUDY_SETTINGS,  "set"),
    ]
    for col, values, prefix in axes:
        by_val: dict[str, list[dict]] = {v: [] for v in values}
        for r in rows:
            v = r.get(col, "")
            if v in by_val:
                by_val[v].append(r)
        for val, subset in by_val.items():
            path = OUT_DIR / "data" / f"subset_{prefix}_{val}.tsv"
            with open(path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                        extrasaction="ignore")
                writer.writeheader()
                writer.writerows(subset)
            print(f"  {path.name:<40} {len(subset):>4} rows")


# ── HTML generation ────────────────────────────────────────────────────────────

def build_html(rows: list[dict], cache: dict, alignment: dict) -> str:  # noqa: C901
    js_rows = []
    for r in rows:
        bp      = r["BioProject"]
        cached  = cache.get(bp, {})
        align   = alignment.get(bp, {})
        abstract     = (cached.get("abstract")     or "").strip()
        methods_text = (cached.get("methods_text") or "").strip()
        pmcid        = (cached.get("pmcid")        or "").strip()
        pub_title    = ""
        pub = r.get("primary_publication", "")
        if pub and "] " in pub:
            pub_title = pub.split("] ", 1)[1]
        js_rows.append({
            "bp":          bp,
            "title":       r.get("title", ""),
            "description": (r.get("description") or "").strip(),
            "modes":       r.get("modes", ""),
            "n_coinf":     r.get("n_coinf", 0),
            "n_single":    r.get("n_single", 0),
            "coinf_rate":  r.get("coinf_rate", ""),
            "treatment":   r.get("treatment", ""),
            "setting":     r.get("study_setting", ""),
            "treat_kws":   r.get("treatment_keywords", ""),
            "set_kws":     r.get("setting_keywords", ""),
            "pmid":        r.get("primary_pmid", ""),
            "pmcid":       pmcid,
            "pub_date":    r.get("primary_pub_date", ""),
            "pub_title":   pub_title,
            "abstract":    abstract,
            "methods":     methods_text[:6000],
            "bp_date":     r.get("bp_submission_date", ""),
            "primaries":   r.get("primaries", ""),
            "secondaries": r.get("secondaries", ""),
            "sra_top":          align.get("sra_top", ""),
            "stat_top":         align.get("stat_tops", [""])[0] if align.get("stat_tops") else "",
            "meta_pathogen":    align.get("meta_pathogen", ""),
            "align_flags":      align.get("flags", {}),
            "stat_tops":        align.get("stat_tops", []),
            "align_n_conflict": align.get("n_conflict", 0),
            "llm_treatment":    r.get("llm_treatment", ""),
            "llm_setting":      r.get("llm_study_setting", ""),
            "llm_pathogen":     r.get("llm_named_pathogen", ""),
            "llm_host":         "",
            "llm_stat_present": r.get("llm_stat_present", ""),
            "llm_stat_dominant":r.get("llm_stat_dominant", ""),
            "llm_confidence":   r.get("llm_confidence", ""),
            "llm_rationale":    r.get("llm_rationale", ""),
        })

    js_data     = json.dumps(js_rows, ensure_ascii=False).replace('</', '<\\/')
    treat_colours = json.dumps(TREATMENT_COLOURS)
    set_colours   = json.dumps(SETTING_COLOURS)
    treat_types   = json.dumps(TREATMENT_TYPES)
    set_types     = json.dumps(STUDY_SETTINGS)
    total       = len(rows)

    template = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>BioProject Design Review</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;font-size:13px;background:#f4f4f4;color:#222}
#hdr{position:sticky;top:0;z-index:20;background:#fff;border-bottom:1px solid #ddd;padding:10px 14px}
#hdr h1{font-size:15px;font-weight:700;margin-bottom:7px}
.filter-row{display:flex;gap:5px;flex-wrap:wrap;align-items:center;margin-bottom:4px}
.filter-label{font-size:10px;font-weight:700;color:#888;text-transform:uppercase;margin-right:2px;white-space:nowrap}
.fbtn{border:2px solid #bbb;background:#fff;border-radius:4px;padding:3px 9px;cursor:pointer;font-size:11px;font-weight:600;transition:.15s}
.fbtn:hover{opacity:.85}
.fbtn.active{color:#fff;border-color:transparent}
#hdr-row2{display:flex;align-items:center;gap:12px;margin-top:6px}
#stats{font-size:11px;color:#777;flex:1}
.hdr-btn{background:#2c3e50;color:#fff;border:none;border-radius:4px;padding:5px 11px;cursor:pointer;font-size:11px}
table{width:100%;border-collapse:collapse;background:#fff}
th{position:sticky;top:52px;background:#f0f0f0;text-align:left;padding:5px 8px;font-size:11px;font-weight:700;color:#555;border-bottom:1px solid #ccc;white-space:nowrap}
.fth{position:sticky;top:88px;background:#e8e8e8;padding:3px 4px;border-bottom:2px solid #bbb}
.fth select,.fth input{width:100%;font-size:10px;padding:2px 3px;border:1px solid #ccc;border-radius:3px;background:#fff;box-sizing:border-box}
.fth input{min-width:60px}
#col-panel{display:none;position:absolute;top:52px;right:14px;z-index:30;background:#fff;border:1px solid #ccc;border-radius:6px;box-shadow:0 4px 16px rgba(0,0,0,.15);padding:10px 14px;min-width:200px}
#col-panel h4{font-size:11px;color:#888;text-transform:uppercase;margin-bottom:8px;font-weight:700}
.col-cb{display:flex;align-items:center;gap:6px;font-size:12px;padding:2px 0;cursor:pointer;user-select:none}
.col-cb input{cursor:pointer}
td{padding:5px 8px;vertical-align:middle;border-bottom:1px solid #eee;font-size:12px}
tr.data-row:hover td{background:#eef4ff;cursor:pointer}
tr.data-row.reviewed td{background:#f0fff4}
tr.data-row.reviewed:hover td{background:#d6f5e3}
.bp-a{font-weight:700;color:#2980b9;text-decoration:none;font-size:11px}
.bp-a:hover{text-decoration:underline}
.badge{display:inline-block;padding:2px 7px;border-radius:10px;font-size:10px;font-weight:700;color:#fff;white-space:nowrap}
.dot-yes{color:#27ae60;font-weight:700}
.dot-no{color:#ccc}
.rate-hi{color:#c0392b;font-weight:700}
select.ann-sel{font-size:10px;padding:1px 2px;border:1px solid #ccc;border-radius:3px;width:115px;display:block;margin-bottom:2px}
#overlay{display:none;position:fixed;inset:0;background:rgba(0,0,0,.45);z-index:100;overflow-y:auto;padding:24px 16px}
#modal{background:#fff;border-radius:8px;max-width:860px;margin:0 auto;box-shadow:0 8px 40px rgba(0,0,0,.25)}
#modal-hdr{display:flex;align-items:center;gap:8px;padding:12px 16px;border-bottom:1px solid #ddd;background:#f8f8f8;border-radius:8px 8px 0 0}
#modal-hdr .nav-btn{border:1px solid #ccc;background:#fff;border-radius:4px;padding:4px 10px;cursor:pointer;font-size:12px}
#modal-hdr .nav-btn:disabled{opacity:.35;cursor:default}
#modal-counter{font-size:12px;color:#888;margin:0 4px;flex:1}
#modal-close{border:none;background:#e74c3c;color:#fff;border-radius:4px;padding:4px 10px;cursor:pointer;font-size:12px}
#modal-body{padding:16px}
.m-meta{display:grid;grid-template-columns:1fr 1fr;gap:6px 16px;margin-bottom:14px;font-size:12px}
.m-meta .lbl{color:#888;font-size:10px;font-weight:700;text-transform:uppercase;margin-bottom:1px}
.m-meta .val{color:#222}
.m-section{margin-bottom:14px}
.m-section h3{font-size:11px;font-weight:700;color:#555;text-transform:uppercase;letter-spacing:.5px;margin-bottom:6px;padding-bottom:4px;border-bottom:1px solid #eee}
.m-text{font-size:12px;line-height:1.6;color:#333;white-space:pre-wrap;word-break:break-word;max-height:220px;overflow-y:auto;padding:8px;background:#fafafa;border:1px solid #eee;border-radius:4px}
.m-absent{font-size:12px;color:#bbb;font-style:italic}
.ann-row{display:flex;align-items:flex-start;gap:16px;margin-bottom:10px;flex-wrap:wrap}
.ann-row label{font-size:11px;font-weight:700;color:#555;display:block;margin-bottom:4px}
select.ann-sel-lg{font-size:12px;padding:4px 6px;border:1px solid #ccc;border-radius:4px;width:170px}
textarea.ann-notes{font-size:12px;padding:6px;border:1px solid #ccc;border-radius:4px;min-height:60px;resize:vertical;width:100%}
input.ann-doi{font-size:12px;padding:5px 7px;border:1px solid #ccc;border-radius:4px;width:100%}
.pathogen-grid{display:grid;grid-template-columns:1fr 2fr;gap:6px 16px;font-size:12px}
.sec-tag{display:inline-block;background:#eef;border:1px solid #c8d;border-radius:3px;padding:1px 6px;margin:2px 2px 2px 0;font-size:11px;color:#444}
th.sortable{cursor:pointer;user-select:none}
th.sortable:hover{background:#e4e4e4}
.sort-ind{font-size:9px;opacity:.55;margin-left:3px}
</style>
</head>
<body>

<div id="hdr">
  <div id="hdr-row1" style="display:flex;align-items:center;gap:10px;margin-bottom:5px">
    <h1 style="margin:0">BioProject Design Review &mdash; __TOTAL__ BioProjects</h1>
    <input id="text-search" type="text" placeholder="&#128269; Search title / BioProject…"
           style="font-size:12px;padding:4px 8px;border:1px solid #ccc;border-radius:4px;width:260px"
           oninput="setColFilter('text',this.value)">
  </div>
  <div id="hdr-row2">
    <span id="stats"></span>
    <button class="hdr-btn" onclick="toggleColPanel()">&#9776; Columns</button>
    <button class="hdr-btn" onclick="exportTSV()">&#11015; Export annotations (TSV)</button>
  </div>
</div>

<div id="col-panel">
  <h4>Show / hide columns</h4>
</div>
<style id="col-vis-style"></style>

<table>
  <thead>
    <tr>
      <th class="sortable col-bp"        data-col="bp"             onclick="sortBy('bp')">BioProject <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-title"     data-col="title"          onclick="sortBy('title')">Title <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-n-coinf"    data-col="n_coinf"  style="text-align:center" onclick="sortBy('n_coinf')">n coinf <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-n-single"   data-col="n_single" style="text-align:center" onclick="sortBy('n_single')">n single <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-rate"      data-col="rate"    style="text-align:center" onclick="sortBy('rate')">Rate <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-treatment" data-col="treatment"      onclick="sortBy('treatment')">Treatment <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-setting"   data-col="setting"        onclick="sortBy('setting')">Setting <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-abs"       data-col="abs" style="text-align:center"  onclick="sortBy('abs')">Abs <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-meth"      data-col="meth" style="text-align:center" onclick="sortBy('meth')">Meth <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-pmid"      data-col="pmid" style="text-align:center" onclick="sortBy('pmid')">PMID <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-llm-treatment" data-col="llm_treatment" onclick="sortBy('llm_treatment')" title="LLM treatment classification">LLM treat <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-llm-setting"   data-col="llm_setting"   onclick="sortBy('llm_setting')"   title="LLM setting classification">LLM set <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-llm-pathogen"  data-col="llm_pathogen"  onclick="sortBy('llm_pathogen')"  title="LLM-named pathogen">LLM pathogen <span class="sort-ind">&#8645;</span></th>
      <th class="col-llm-host" title="LLM-named host (pending)">LLM host</th>
      <th class="sortable col-sra"       data-col="sra_top"        onclick="sortBy('sra_top')" title="Most common SRA library organism">SRA library <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-stat"      data-col="stat_top"       onclick="sortBy('stat_top')" title="Most common STAT top pathogen">STAT top <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-meta"      data-col="meta_pathogen"  onclick="sortBy('meta_pathogen')" title="Pathogen found in BioProject text">Metadata pathogen <span class="sort-ind">&#8645;</span></th>
      <th class="sortable col-mismatch"  data-col="align_n_conflict" onclick="sortBy('align_n_conflict')" title="Alignment conflict flag">Mismatch <span class="sort-ind">&#8645;</span></th>
      <th class="col-manual">Manual override</th>
    </tr>
    <tr id="filter-row">
      <th class="fth col-bp"><input placeholder="accession…" oninput="setColFilter('bp',this.value)"></th>
      <th class="fth col-title"><input placeholder="title keyword…" oninput="setColFilter('title',this.value)"></th>
      <th class="fth col-n-coinf"></th>
      <th class="fth col-n-single"></th>
      <th class="fth col-rate"></th>
      <th class="fth col-treatment"><select id="cf-treatment" onchange="setColFilter('treatment',this.value)"><option value="all">all treatments</option></select></th>
      <th class="fth col-setting"><select id="cf-setting" onchange="setColFilter('setting',this.value)"><option value="all">all settings</option></select></th>
      <th class="fth col-abs"><select onchange="setColFilter('abs',this.value)"><option value="all">abs?</option><option value="yes">yes</option><option value="no">no</option></select></th>
      <th class="fth col-meth"><select onchange="setColFilter('meth',this.value)"><option value="all">meth?</option><option value="yes">yes</option><option value="no">no</option></select></th>
      <th class="fth col-pmid"><select onchange="setColFilter('pmid',this.value)"><option value="all">pmid?</option><option value="yes">yes</option><option value="no">no</option></select></th>
      <th class="fth col-llm-treatment"><select id="cf-llm-treatment" onchange="setColFilter('llm_treatment',this.value)"><option value="all">all LLM treat</option></select></th>
      <th class="fth col-llm-setting"><select id="cf-llm-setting" onchange="setColFilter('llm_setting',this.value)"><option value="all">all LLM set</option></select></th>
      <th class="fth col-llm-pathogen"></th>
      <th class="fth col-llm-host"></th>
      <th class="fth col-sra"><input placeholder="organism…" oninput="setColFilter('organism',this.value)"></th>
      <th class="fth col-stat"></th>
      <th class="fth col-meta"></th>
      <th class="fth col-mismatch"><select onchange="setColFilter('mismatch',this.value)"><option value="all">all</option><option value="conflict">conflict</option><option value="meta">meta diff</option><option value="ok">ok</option><option value="undet">undet</option></select></th>
      <th class="fth col-manual"></th>
    </tr>
  </thead>
  <tbody id="tbody"></tbody>
</table>

<div id="overlay" onclick="overlayClick(event)">
  <div id="modal">
    <div id="modal-hdr">
      <button class="nav-btn" id="prev-btn" onclick="navigateDetail(-1)">&#8592; Prev</button>
      <span id="modal-counter"></span>
      <button class="nav-btn" id="next-btn" onclick="navigateDetail(1)">Next &#8594;</button>
      <button id="modal-close" onclick="closeModal()">&#x2715; Close</button>
    </div>
    <div id="modal-body"></div>
  </div>
</div>

<script>
const DATA             = __DATA__;
const TREAT_COLOURS    = __TREAT_COLOURS__;
const SET_COLOURS      = __SET_COLOURS__;
const TREATMENT_TYPES  = __TREAT_TYPES__;
const STUDY_SETTINGS   = __SET_TYPES__;
const MANUAL_TREATMENTS = ['', ...TREATMENT_TYPES];
const MANUAL_SETTINGS   = ['', ...STUDY_SETTINGS];

let filteredData   = DATA.slice();
let columnFilters  = { treatment:'all', setting:'all', llm_treatment:'all', llm_setting:'all',
                       abs:'all', meth:'all', pmid:'all', mismatch:'all', text:'', bp:'', title:'', organism:'' };
let detailIdx      = 0;
let sortCol        = null;
let sortDir        = 1;

const SORT_KEYS = {
  bp:            r => r.bp,
  title:         r => r.title.toLowerCase(),
  modes:         r => r.modes,
  n_coinf:       r => +r.n_coinf,
  n_single:      r => +r.n_single,
  rate:          r => r.coinf_rate !== '' ? +r.coinf_rate : -1,
  treatment:     r => r.treatment,
  setting:       r => r.setting,
  llm_treatment: r => r.llm_treatment,
  llm_setting:   r => r.llm_setting,
  llm_pathogen:  r => (r.llm_pathogen || '').toLowerCase(),
  abs:           r => r.abstract ? 1 : 0,
  meth:          r => r.methods  ? 1 : 0,
  pmid:          r => r.pmid     ? 1 : 0,
};

/* ── Annotations ── */
function getAnns() {
  try { return JSON.parse(localStorage.getItem('bp_anns') || '{}'); }
  catch { return {}; }
}
function saveAnn(bp, field, val) {
  const a = getAnns();
  if (!a[bp]) a[bp] = {};
  a[bp][field] = val;
  localStorage.setItem('bp_anns', JSON.stringify(a));
  updateStats();
  const row = document.getElementById('row-' + bp);
  if (row) {
    const ann = a[bp];
    const reviewed = ann.manual_treatment || ann.manual_setting || ann.doi || ann.notes;
    row.className = 'data-row' + (reviewed ? ' reviewed' : '');
    const sel = row.querySelector('select[data-field="' + field + '"]');
    if (sel) sel.value = val;
  }
}
function listSelChange(el)       { saveAnn(el.dataset.bp, el.dataset.field, el.value); }
function detTreatSelChange(el)   { saveAnn(el.dataset.bp, 'manual_treatment', el.value); }
function detSetSelChange(el)     { saveAnn(el.dataset.bp, 'manual_setting',   el.value); }
function detNotesChange(el)    { saveAnn(el.dataset.bp, 'notes', el.value); }
function detDoiChange(el)      { saveAnn(el.dataset.bp, 'doi',   el.value); }

/* ── Sort ── */
function applySort() {
  if (!sortCol) return;
  const key = SORT_KEYS[sortCol];
  filteredData.sort((a, b) => {
    const va = key(a), vb = key(b);
    if (va < vb) return -sortDir;
    if (va > vb) return  sortDir;
    return 0;
  });
}
function sortBy(col) {
  sortDir = sortCol === col ? -sortDir : 1;
  sortCol = col;
  document.querySelectorAll('th.sortable .sort-ind').forEach(el => {
    el.textContent = el.parentElement.dataset.col === col
      ? (sortDir === 1 ? '▲' : '▼') : '⇅';
  });
  applySort();
  renderList();
}

/* ── Column filters ── */
function setColFilter(key, val) {
  columnFilters[key] = val;
  applyAllFilters();
}
function applyAllFilters() {
  const f = columnFilters;
  filteredData = DATA.filter(r => {
    if (f.treatment    !== 'all' && r.treatment    !== f.treatment)    return false;
    if (f.setting      !== 'all' && r.setting      !== f.setting)      return false;
    if (f.llm_treatment !== 'all' && r.llm_treatment !== f.llm_treatment) return false;
    if (f.llm_setting  !== 'all' && r.llm_setting  !== f.llm_setting)  return false;
    if (f.abs  !== 'all' && (f.abs  === 'yes') !== !!r.abstract) return false;
    if (f.meth !== 'all' && (f.meth === 'yes') !== !!r.methods)  return false;
    if (f.pmid !== 'all' && (f.pmid === 'yes') !== !!r.pmid)     return false;
    if (f.mismatch !== 'all') {
      const fl = r.align_flags || {};
      const nC = (fl.stat_conflict||0)+(fl.both_conflict||0)+(fl.stat_vs_meta_conflict||0);
      const nM = fl.meta_conflict || 0;
      const nA = fl.all_match || 0;
      if (f.mismatch === 'conflict' && nC === 0) return false;
      if (f.mismatch === 'meta'     && nM === 0) return false;
      if (f.mismatch === 'ok'       && nA === 0) return false;
      if (f.mismatch === 'undet'    && (fl.stat_undetected||0) === 0) return false;
    }
    if (f.bp && !r.bp.toLowerCase().includes(f.bp.toLowerCase())) return false;
    if (f.title) {
      const q = f.title.toLowerCase();
      if (!r.title.toLowerCase().includes(q) && !r.description.toLowerCase().includes(q)) return false;
    }
    if (f.text) {
      const q = f.text.toLowerCase();
      if (!r.bp.toLowerCase().includes(q) && !r.title.toLowerCase().includes(q)) return false;
    }
    if (f.organism) {
      const q = f.organism.toLowerCase();
      const haystack = (r.sra_top + ' ' + r.stat_top + ' ' + r.meta_pathogen + ' ' + r.llm_pathogen).toLowerCase();
      if (!haystack.includes(q)) return false;
    }
    return true;
  });
  applySort();
  renderList();
}

function updateStats() {
  const anns = getAnns();
  const rev = filteredData.filter(r => {
    const a = anns[r.bp];
    return a && (a.manual_treatment || a.manual_setting || a.doi || a.notes);
  }).length;
  document.getElementById('stats').textContent =
    filteredData.length + ' BioProjects shown · ' + rev + ' manually reviewed';
}

/* ── List view ── */
function esc(s) {
  return String(s || '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
function badge(d, colours) {
  const c = colours[d] || '#aaa';
  return '<span class="badge" style="background:' + c + '">' + esc(d) + '</span>';
}
function dot(present) {
  return present ? '<span class="dot-yes" title="present">&#10003;</span>'
                 : '<span class="dot-no" title="absent">&mdash;</span>';
}
function sraStyle(modes) {
  if (modes === 'hal') return {color:'#27ae60', label:'host plant'};
  if (modes === 'mal') return {color:'#c0392b', label:'pathogen'};
  return {color:'#2980b9', label:'both modes'};
}
function alignBadge(r) {
  const f  = r.align_flags || {};
  const nC = (f.stat_conflict||0) + (f.both_conflict||0) + (f.stat_vs_meta_conflict||0);
  const nM = (f.meta_conflict||0);
  const nA = (f.all_match||0);
  const nU = (f.stat_undetected||0);
  if (nC > 0) return '<span style="background:#e74c3c;color:#fff;border-radius:3px;padding:1px 5px;font-size:11px" title="STAT detects different organism">' + nC + ' conflict</span>';
  if (nM > 0) return '<span style="background:#e67e22;color:#fff;border-radius:3px;padding:1px 5px;font-size:11px" title="metadata mentions different pathogen">' + nM + ' meta</span>';
  if (nA > 0) return '<span style="background:#27ae60;color:#fff;border-radius:3px;padding:1px 5px;font-size:11px">' + nA + '&#x2713;</span>';
  return '<span style="color:#aaa;font-size:11px">' + nU + ' undet.</span>';
}

function renderList() {
  const anns  = getAnns();
  const tbody = document.getElementById('tbody');
  tbody.innerHTML = filteredData.map((r, i) => {
    const ann = anns[r.bp] || {};
    const reviewed = ann.manual_treatment || ann.manual_setting || ann.doi || ann.notes;
    const cls  = reviewed ? 'data-row reviewed' : 'data-row';
    const rate = r.coinf_rate !== ''
      ? '<span class="' + (parseFloat(r.coinf_rate) > 0.1 ? 'rate-hi' : '') + '">'
        + (parseFloat(r.coinf_rate)*100).toFixed(0) + '%</span>'
      : '&mdash;';
    const treatOpts = MANUAL_TREATMENTS.map(d =>
      '<option value="' + d + '"' + (ann.manual_treatment === d ? ' selected' : '') + '>'
      + (d || '— treatment —') + '</option>'
    ).join('');
    const setOpts = MANUAL_SETTINGS.map(d =>
      '<option value="' + d + '"' + (ann.manual_setting === d ? ' selected' : '') + '>'
      + (d || '— setting —') + '</option>'
    ).join('');
    return '<tr class="' + cls + '" id="row-' + r.bp + '" onclick="rowClick(event,' + i + ')">'
      + '<td class="col-bp"><a class="bp-a" href="https://www.ncbi.nlm.nih.gov/bioproject/' + r.bp
        + '" target="_blank" onclick="event.stopPropagation()">' + r.bp + '</a></td>'
      + '<td class="col-title" style="max-width:200px">' + esc(r.title.slice(0,70)) + (r.title.length > 70 ? '…' : '') + '</td>'
      + '<td class="col-n-coinf"  style="text-align:center">' + r.n_coinf  + '</td>'
      + '<td class="col-n-single" style="text-align:center">' + r.n_single + '</td>'
      + '<td class="col-rate" style="text-align:center">' + rate + '</td>'
      + '<td class="col-treatment">' + badge(r.treatment, TREAT_COLOURS) + '</td>'
      + '<td class="col-setting">'   + badge(r.setting,   SET_COLOURS)   + '</td>'
      + '<td class="col-abs"  style="text-align:center">' + dot(!!r.abstract) + '</td>'
      + '<td class="col-meth" style="text-align:center">' + dot(!!r.methods)  + '</td>'
      + '<td class="col-pmid" style="text-align:center">' + dot(!!r.pmid)     + '</td>'
      + (r.llm_treatment
          ? '<td class="col-llm-treatment">' + badge(r.llm_treatment, TREAT_COLOURS)
            + (r.llm_treatment !== r.treatment ? '&nbsp;<span style="color:#e74c3c;font-size:9px">&#8800;</span>' : '') + '</td>'
          : '<td class="col-llm-treatment"><span style="color:#ccc;font-size:11px">&mdash;</span></td>')
      + (r.llm_setting
          ? '<td class="col-llm-setting">' + badge(r.llm_setting, SET_COLOURS)
            + (r.llm_setting !== r.setting ? '&nbsp;<span style="color:#e74c3c;font-size:9px">&#8800;</span>' : '') + '</td>'
          : '<td class="col-llm-setting"><span style="color:#ccc;font-size:11px">&mdash;</span></td>')
      + '<td class="col-llm-pathogen" style="font-size:11px;font-style:italic;max-width:130px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="' + esc(r.llm_pathogen) + '">'
        + (r.llm_pathogen ? esc(r.llm_pathogen) : '<span style="color:#ccc">&mdash;</span>') + '</td>'
      + '<td class="col-llm-host" style="font-size:11px;color:#bbb">&mdash;</td>'
      + (function(){
          const s = sraStyle(r.modes);
          return '<td class="col-sra" style="font-size:11px;max-width:130px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="' + esc(r.sra_top) + ' (' + s.label + ')">'
            + (r.sra_top ? '<span style="color:' + s.color + ';font-style:italic">' + esc(r.sra_top) + '</span>' : '<span style="color:#bbb">—</span>')
            + '</td>';
        })()
      + '<td class="col-stat" style="font-size:11px;max-width:120px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="' + esc(r.stat_top) + '">'
        + (r.stat_top ? esc(r.stat_top) : '<span style="color:#bbb">—</span>') + '</td>'
      + '<td class="col-meta" style="font-size:11px;max-width:120px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="' + esc(r.meta_pathogen) + '">'
        + (r.meta_pathogen ? esc(r.meta_pathogen) : '<span style="color:#bbb">—</span>') + '</td>'
      + '<td class="col-mismatch">' + alignBadge(r) + '</td>'
      + '<td class="col-manual" onclick="event.stopPropagation()">'
        + '<select class="ann-sel" data-bp="' + r.bp + '" data-field="manual_treatment" onchange="listSelChange(this)">' + treatOpts + '</select>'
        + '<select class="ann-sel" data-bp="' + r.bp + '" data-field="manual_setting"   onchange="listSelChange(this)">' + setOpts + '</select>'
        + '</td>'
      + '</tr>';
  }).join('');
  updateStats();
}

function rowClick(event, idx) {
  if (event.target.tagName === 'SELECT' || event.target.tagName === 'OPTION'
      || event.target.tagName === 'A') return;
  openDetail(idx);
}

/* ── Detail modal ── */
function openDetail(idx) {
  detailIdx = idx;
  renderDetail();
  document.getElementById('overlay').style.display = 'block';
  document.body.style.overflow = 'hidden';
}
function closeModal() {
  document.getElementById('overlay').style.display = 'none';
  document.body.style.overflow = '';
}
function overlayClick(e) {
  if (e.target === document.getElementById('overlay')) closeModal();
}
function navigateDetail(dir) {
  const next = detailIdx + dir;
  if (next >= 0 && next < filteredData.length) { detailIdx = next; renderDetail(); }
}
document.addEventListener('keydown', e => {
  if (document.getElementById('overlay').style.display !== 'block') return;
  if (e.key === 'Escape')     closeModal();
  if (e.key === 'ArrowRight') navigateDetail(1);
  if (e.key === 'ArrowLeft')  navigateDetail(-1);
});

function renderDetail() {
  const r    = filteredData[detailIdx];
  const anns = getAnns();
  const ann  = anns[r.bp] || {};

  document.getElementById('modal-counter').textContent =
    (detailIdx + 1) + ' / ' + filteredData.length;
  document.getElementById('prev-btn').disabled = detailIdx === 0;
  document.getElementById('next-btn').disabled = detailIdx === filteredData.length - 1;

  const rate = r.coinf_rate !== '' ? (parseFloat(r.coinf_rate)*100).toFixed(1) + '%' : '—';
  const pmidHtml = r.pmid
    ? '<a href="https://pubmed.ncbi.nlm.nih.gov/' + r.pmid + '/" target="_blank"'
      + ' style="color:#8e44ad">PMID ' + r.pmid + '</a>'
      + (r.pub_date ? '  ·  ' + r.pub_date.slice(0,7) : '')
    : '<span class="m-absent">no PMID found</span>';
  const pmcHtml = r.pmcid
    ? '<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC' + r.pmcid
      + '/" target="_blank" style="color:#27ae60">PMC' + r.pmcid + '</a>'
    : '<span class="m-absent">—</span>';

  const treatOpts = MANUAL_TREATMENTS.map(d =>
    '<option value="' + d + '"' + (ann.manual_treatment === d ? ' selected' : '') + '>'
    + (d || '— not reviewed —') + '</option>'
  ).join('');
  const setOpts = MANUAL_SETTINGS.map(d =>
    '<option value="' + d + '"' + (ann.manual_setting === d ? ' selected' : '') + '>'
    + (d || '— not reviewed —') + '</option>'
  ).join('');

  document.getElementById('modal-body').innerHTML = `
    <div class="m-meta">
      <div>
        <div class="lbl">BioProject &nbsp;&middot;&nbsp; submitted</div>
        <div class="val">
          <a href="https://www.ncbi.nlm.nih.gov/bioproject/${r.bp}" target="_blank"
             style="color:#2980b9;font-weight:700">${r.bp}</a>
          &nbsp;&middot;&nbsp; ${r.bp_date ? r.bp_date.slice(0,7) : '?'}
          &nbsp;&middot;&nbsp; mode: <strong>${esc(r.modes)}</strong>
        </div>
      </div>
      <div>
        <div class="lbl">Co-infection</div>
        <div class="val">${r.n_coinf} coinf / ${r.n_single} single = <strong>${rate}</strong></div>
      </div>
      <div>
        <div class="lbl">PubMed</div>
        <div class="val">${pmidHtml}${r.pub_title ? ' &nbsp;&middot;&nbsp; <em style="color:#555;font-size:11px">' + esc(r.pub_title) + '</em>' : ''}</div>
      </div>
      <div>
        <div class="lbl">PMC full text</div>
        <div class="val">${pmcHtml}</div>
      </div>
      <div style="grid-column:1/-1">
        <div class="lbl">BioProject title</div>
        <div class="val" style="font-size:13px;font-weight:500">${esc(r.title)}</div>
      </div>
    </div>
    ${r.description ? `<div class="m-section">
      <h3>BioProject description</h3>
      <div class="m-text" style="font-size:12px;color:#333">${esc(r.description)}</div>
    </div>` : ''}
    <div class="m-section">
      <h3>LLM classification</h3>
      <div class="pathogen-grid" style="margin-bottom:8px">
        <div>
          <div class="lbl">Treatment</div>
          <div class="val">${r.llm_treatment
            ? badge(r.llm_treatment, TREAT_COLOURS)
              + (r.llm_treatment !== r.treatment
                  ? '&nbsp;<span style="color:#e74c3c;font-size:10px">&#8800;&thinsp;kw:&thinsp;' + esc(r.treatment) + '</span>'
                  : '&nbsp;<span style="color:#27ae60;font-size:10px">&#10003;</span>')
            : '<em style="color:#bbb">—</em>'}</div>
        </div>
        <div>
          <div class="lbl">Setting</div>
          <div class="val">${r.llm_setting
            ? badge(r.llm_setting, SET_COLOURS)
              + (r.llm_setting !== r.setting
                  ? '&nbsp;<span style="color:#e74c3c;font-size:10px">&#8800;&thinsp;kw:&thinsp;' + esc(r.setting) + '</span>'
                  : '&nbsp;<span style="color:#27ae60;font-size:10px">&#10003;</span>')
            : '<em style="color:#bbb">—</em>'}</div>
        </div>
        <div>
          <div class="lbl">Named pathogen</div>
          <div class="val" style="font-size:12px;font-style:italic">${r.llm_pathogen || '<em style="color:#bbb">—</em>'}</div>
        </div>
        <div>
          <div class="lbl">Named host</div>
          <div class="val" style="font-size:12px;font-style:italic;color:#bbb">— (pending)</div>
        </div>
        ${r.llm_pathogen ? `<div>
          <div class="lbl">STAT present / dominant</div>
          <div class="val" style="font-size:12px">${
            (r.llm_stat_present === 'True'  ? '<span style="color:#27ae60">present</span>'   :
             r.llm_stat_present === 'False' ? '<span style="color:#e74c3c">absent</span>'    : '—')
            + ' / ' +
            (r.llm_stat_dominant === 'True'  ? '<span style="color:#27ae60">dominant</span>'      :
             r.llm_stat_dominant === 'False' ? '<span style="color:#e67e22">not dominant</span>'  : '—')
          }</div>
        </div>
        <div>
          <div class="lbl">Confidence</div>
          <div class="val" style="font-size:12px;color:#555">${esc(r.llm_confidence || '—')}</div>
        </div>` : '<div></div><div></div>'}
      </div>
      ${r.llm_rationale ? '<div class="m-text" style="font-size:12px;color:#444;background:#fffdf0;border-color:#e8e0a0">' + esc(r.llm_rationale) + '</div>' : '<div class="m-absent">No LLM rationale</div>'}
    </div>
    <div class="m-section">
      <h3>SRA / STAT evidence</h3>
      <div class="pathogen-grid">
        <div>
          <div class="lbl">SRA library organism</div>
          <div class="val" style="font-size:12px;font-style:italic">${r.sra_top ? esc(r.sra_top) : '<em style="color:#bbb">—</em>'}</div>
        </div>
        <div>
          <div class="lbl">STAT top detections (by run count)</div>
          <div class="val" style="font-size:12px">${r.stat_tops && r.stat_tops.length
            ? r.stat_tops.map(s => '<span class="sec-tag">' + esc(s) + '</span>').join(' ')
            : '<em style="color:#bbb">—</em>'}</div>
        </div>
        <div>
          <div class="lbl">CRYPT primary</div>
          <div class="val" style="font-size:12px">${r.primaries ? esc(r.primaries) : '<em style="color:#bbb">—</em>'}</div>
        </div>
        <div>
          <div class="lbl">CRYPT secondaries</div>
          <div class="val" style="font-size:12px">${r.secondaries
            ? r.secondaries.split(';').map(s => '<span class="sec-tag">' + esc(s.trim()) + '</span>').join(' ')
            : '<em style="color:#bbb">none</em>'}</div>
        </div>
        <div>
          <div class="lbl">Metadata pathogen (text match)</div>
          <div class="val" style="font-size:12px">${r.meta_pathogen ? esc(r.meta_pathogen) : '<em style="color:#bbb">none found</em>'}</div>
        </div>
        <div>
          <div class="lbl">Alignment flags</div>
          <div class="val">${(function(){
            const f = r.align_flags || {};
            const pairs = [
              ['all_match','#27ae60'],['stat_conflict','#e74c3c'],
              ['meta_conflict','#e67e22'],['both_conflict','#c0392b'],
              ['stat_undetected','#95a5a6'],['stat_vs_meta_conflict','#e74c3c'],
              ['no_meta','#bdc3c7']
            ];
            const parts = pairs.filter(([k]) => f[k]).map(([k,c]) =>
              '<span style="background:' + c + ';color:#fff;border-radius:3px;padding:1px 5px;font-size:10px;margin-right:3px">'
              + k + ':' + f[k] + '</span>'
            );
            return parts.length ? parts.join('') : '<em style="color:#bbb">—</em>';
          })()}</div>
        </div>
      </div>
    </div>
    <div class="m-section">
      <h3>Abstract</h3>
      ${r.abstract
        ? '<div class="m-text">' + esc(r.abstract) + '</div>'
        : '<div class="m-absent">No abstract available</div>'}
    </div>
    <div class="m-section">
      <h3>Methods (from PMC full text)</h3>
      ${r.methods
        ? '<div class="m-text">' + esc(r.methods) + '</div>'
        : '<div class="m-absent">Not available (paper not in PMC or no methods section found)</div>'}
    </div>
    <div class="m-section">
      <h3>Manual annotation</h3>
      <div class="ann-row">
        <div>
          <label>Treatment</label>
          <select class="ann-sel-lg" data-bp="${r.bp}" onchange="detTreatSelChange(this)">${treatOpts}</select>
        </div>
        <div>
          <label>Study setting</label>
          <select class="ann-sel-lg" data-bp="${r.bp}" onchange="detSetSelChange(this)">${setOpts}</select>
        </div>
        <div style="flex:1;min-width:200px">
          <label>DOI / URL for methods paper</label>
          <input type="text" class="ann-doi" data-bp="${r.bp}"
                 placeholder="e.g. 10.1111/ppa.12345"
                 value="${esc(ann.doi || '')}"
                 onchange="detDoiChange(this)">
        </div>
      </div>
      <div>
        <label style="font-size:11px;font-weight:700;color:#555;display:block;margin-bottom:4px">Notes</label>
        <textarea class="ann-notes" data-bp="${r.bp}"
                  onchange="detNotesChange(this)">${esc(ann.notes || '')}</textarea>
      </div>
    </div>`;
}

/* ── Populate column filter dropdowns ── */
function buildFilters() {
  function populateSel(id, values, colours, key) {
    const sel = document.getElementById(id);
    if (!sel) return;
    values.forEach(v => {
      const n = DATA.filter(r => r[key] === v).length;
      if (n === 0) return;
      const opt = document.createElement('option');
      opt.value = v;
      opt.textContent = v + ' (' + n + ')';
      sel.appendChild(opt);
    });
  }
  populateSel('cf-treatment',     TREATMENT_TYPES, TREAT_COLOURS, 'treatment');
  populateSel('cf-setting',       STUDY_SETTINGS,  SET_COLOURS,   'setting');
  populateSel('cf-llm-treatment', TREATMENT_TYPES, TREAT_COLOURS, 'llm_treatment');
  populateSel('cf-llm-setting',   STUDY_SETTINGS,  SET_COLOURS,   'llm_setting');
}

/* ── Column visibility ── */
const COL_DEFS = [
  { cls: 'col-bp',           label: 'BioProject' },
  { cls: 'col-title',        label: 'Title' },
  { cls: 'col-n-coinf',      label: 'n coinf' },
  { cls: 'col-n-single',     label: 'n single' },
  { cls: 'col-rate',         label: 'Rate' },
  { cls: 'col-treatment',    label: 'Treatment (kw)' },
  { cls: 'col-setting',      label: 'Setting (kw)' },
  { cls: 'col-abs',          label: 'Abstract' },
  { cls: 'col-meth',         label: 'Methods' },
  { cls: 'col-pmid',         label: 'PMID' },
  { cls: 'col-llm-treatment',label: 'LLM treatment' },
  { cls: 'col-llm-setting',  label: 'LLM setting' },
  { cls: 'col-llm-pathogen', label: 'LLM pathogen' },
  { cls: 'col-llm-host',     label: 'LLM host (pending)' },
  { cls: 'col-sra',          label: 'SRA library' },
  { cls: 'col-stat',         label: 'STAT top' },
  { cls: 'col-meta',         label: 'Metadata pathogen' },
  { cls: 'col-mismatch',     label: 'Mismatch' },
  { cls: 'col-manual',       label: 'Manual override' },
];

let hiddenCols = new Set(JSON.parse(localStorage.getItem('hidden_cols') || '[]'));

function updateColVis() {
  let css = '';
  hiddenCols.forEach(c => { css += '.' + c + '{display:none}'; });
  document.getElementById('col-vis-style').textContent = css;
  localStorage.setItem('hidden_cols', JSON.stringify([...hiddenCols]));
}

function toggleCol(cls) {
  if (hiddenCols.has(cls)) hiddenCols.delete(cls);
  else hiddenCols.add(cls);
  updateColVis();
  const cb = document.querySelector('#col-panel input[data-cls="' + cls + '"]');
  if (cb) cb.checked = !hiddenCols.has(cls);
}

function buildColPanel() {
  const panel = document.getElementById('col-panel');
  COL_DEFS.forEach(({cls, label}) => {
    const lbl = document.createElement('label');
    lbl.className = 'col-cb';
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = !hiddenCols.has(cls);
    cb.dataset.cls = cls;
    cb.onchange = () => toggleCol(cls);
    lbl.appendChild(cb);
    lbl.appendChild(document.createTextNode(label));
    panel.appendChild(lbl);
  });
  updateColVis();
}

function toggleColPanel() {
  const p = document.getElementById('col-panel');
  p.style.display = p.style.display === 'block' ? 'none' : 'block';
}
document.addEventListener('click', e => {
  const p = document.getElementById('col-panel');
  if (p.style.display === 'block' && !p.contains(e.target) && e.target.textContent.trim() !== '☰ Columns')
    p.style.display = 'none';
});

/* ── Export ── */
function exportTSV() {
  const anns = getAnns();
  const bps  = Object.keys(anns).filter(bp =>
    anns[bp].manual_treatment || anns[bp].manual_setting || anns[bp].doi || anns[bp].notes);
  if (!bps.length) { alert('No annotations yet.'); return; }
  const lookup = {};
  DATA.forEach(r => { lookup[r.bp] = { treatment: r.treatment, setting: r.setting }; });
  const lines = ['BioProject\\tauto_treatment\\tauto_setting\\tmanual_treatment\\tmanual_setting\\tdoi\\tnotes',
    ...bps.map(bp => [
      bp,
      (lookup[bp] || {}).treatment || '',
      (lookup[bp] || {}).setting   || '',
      anns[bp].manual_treatment || '',
      anns[bp].manual_setting   || '',
      anns[bp].doi   || '',
      (anns[bp].notes || '').replace(/\\n/g,' ')
    ].join('\\t'))];
  const a = Object.assign(document.createElement('a'), {
    href: URL.createObjectURL(new Blob([lines.join('\\n')], {type:'text/tab-separated-values'})),
    download: 'manual_design_overrides.tsv',
  });
  a.click();
}

buildFilters();
buildColPanel();
renderList();
</script>
</body>
</html>
"""
    return (template
            .replace("__DATA__",          js_data)
            .replace("__TREAT_COLOURS__", treat_colours)
            .replace("__SET_COLOURS__",   set_colours)
            .replace("__TREAT_TYPES__",   treat_types)
            .replace("__SET_TYPES__",     set_types)
            .replace("__TOTAL__",         str(total)))


# ── Entry point ────────────────────────────────────────────────────────────────

def _load_alignment(path: Path) -> dict:
    """Aggregate primary_alignment.tsv to BioProject level."""
    alignment: dict[str, dict] = {}
    if not path.exists():
        return alignment
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            bp   = r["BioProject"]
            flag = r.get("mismatch_flag", "")
            top  = r.get("stat_top_pathogen", "").strip()
            sra  = r.get("sra_organism", "").strip()
            if bp not in alignment:
                alignment[bp] = {
                    "meta_pathogen":    r.get("meta_pathogen", ""),
                    "flags":            Counter(),
                    "stat_top_counter": Counter(),
                    "sra_counter":      Counter(),
                }
            alignment[bp]["flags"][flag] += 1
            if top:
                alignment[bp]["stat_top_counter"][top] += 1
            if sra:
                alignment[bp]["sra_counter"][sra] += 1
    result: dict[str, dict] = {}
    for bp, d in alignment.items():
        flags = dict(d["flags"])
        n_conflict = (flags.get("stat_conflict", 0) + flags.get("both_conflict", 0) +
                      flags.get("stat_vs_meta_conflict", 0))
        result[bp] = {
            "meta_pathogen": d["meta_pathogen"],
            "flags":         flags,
            "stat_tops":     [n for n, _ in d["stat_top_counter"].most_common(3)],
            "sra_top":       d["sra_counter"].most_common(1)[0][0] if d["sra_counter"] else "",
            "n_conflict":    n_conflict,
        }
    return result


def main() -> None:
    if not META_TSV.exists():
        raise SystemExit(f"ERROR: {META_TSV} not found — run 04_filter_meta.py first")

    with open(META_TSV, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    cache     = json.loads(CACHE_JSON.read_text()) if CACHE_JSON.exists() else {}
    alignment = _load_alignment(ALIGN_TSV)

    print(f"Loaded {len(rows)} BioProjects from {META_TSV}")
    print(f"Loaded alignment data for {len(alignment)} BioProjects from {ALIGN_TSV}")

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    print("Writing subset TSVs:")
    write_subsets(rows)

    html_path = OUT_DIR / "review_designs.html"
    html_path.write_text(build_html(rows, cache, alignment), encoding="utf-8")
    print(f"\nReview tool: {html_path}")
    print("Open in browser. Click any row for detail view. Arrows / Esc to navigate/close.")


if __name__ == "__main__":
    main()
