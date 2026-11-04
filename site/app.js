/* app.js -- the crypt browser. Three views over site/data/, drawn client-side.
 *
 * WHY DRAWN AND NOT RENDERED. The batch renderer drew 2,862 PNGs in two hours and a 23-species
 * change to taxon_status.tsv invalidated all of them the same afternoon. So the data ships and
 * the picture is made on demand. Inherited from prep_browser.py, which found that out.
 *
 * THE FLOORS ARE RE-APPLIED HERE, which is the point: a reader can move them and see what a
 * different threshold would have done. The moment one leaves detect.py's own value the header
 * says WHAT-IF, so a number read off a moved slider cannot be mistaken for a pipeline result.
 *
 * THE RE-CALL MUST REPRODUCE detect.py EXACTLY, including the congener rule, which detect.py
 * applies AFTER the three floors: a taxon whose only membership route is sharing a genus with a
 * pathogen is not reportable. Dropping it called 8 runs differently. That was found by diffing
 * this viewer against the pipeline, not by reading the code, so leave the flag in.
 */
'use strict';

/* ── the palette, which belongs to the figures ───────────────────────────────────────────── */
const TIER_COL = {host:'#3e9e3b', pathogen_crop:'#8c0404', pathogen_other:'#e77e7e',
                  resident_crop:'#22127e', resident_other:'#add0e4', no_record:'#7a7a7a'};
const TIER_LAB = {host:'host (plant)', pathogen_crop:'host pathogen',
                  pathogen_other:'non-host pathogen',
                  resident_crop:'resident, recorded on this crop',
                  resident_other:'resident, recorded on other plants',
                  no_record:'resident, no plant host record'};
const TIER_HEAD = {host:'Host', pathogen_crop:'Host pathogen', pathogen_other:'Non-host pathogen',
                   resident_crop:'Recorded on this crop',
                   resident_other:'Recorded on other plants', no_record:'No plant host record'};
/* sequential, for magnitude. One hue, light to dark, validated. */
const RAMP = ['#dba0a0','#cf7b7b','#bf5757','#a53636','#8c0404'];
const SITE_COL = '#5b6b7a';

/* ── decision-row field offsets. The rows are arrays to keep runs.json small. ─────────────── */
const SP=0, TIER=1, TROPH=2, EVID=3, DECL=4, COV=5, SAT=6, RED=7, DIST=8, READS=9, CONG=10;

const S = {meta:null, sites:null, species:null, evidence:null, runs:null, world:null,
           view:'map', cov:0.02, sat:0.05, dist:1e4, runMode:'rows', mapMode:'extent',
           rankedAll:false, sel:{site:null, species:null, run:null}, map:{k:1, tx:0, ty:0}};

const $ = s => document.querySelector(s);
const el = (t, c, h) => { const e = document.createElement(t); if (c) e.className = c;
                          if (h !== undefined) e.innerHTML = h; return e; };
const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c =>
  ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const n0 = v => v == null ? '--' : Number(v).toLocaleString('en', {maximumFractionDigits:0});
const pct = v => v == null ? '--' : (v * 100 < 0.01 ? (v * 100).toExponential(1)
                                                    : (v * 100).toFixed(2)) + '%';
const sig = (v, d=3) => v == null ? '--' : Number(v).toPrecision(d).replace(/\.?0+$/, '');
const ital = s => `<i>${esc(s)}</i>`;

/* ── a shared log scale, because every measure here spans orders of magnitude ─────────────── */
function logScale(lo, hi, w) {
  const a = Math.log10(lo), b = Math.log10(hi);
  return v => v == null || v <= 0 ? 0 : Math.max(0, Math.min(w, (Math.log10(v) - a) / (b - a) * w));
}

/* ── tooltip ─────────────────────────────────────────────────────────────────────────────── */
const tip = {
  show(html, ev) {
    const t = $('#tip'); t.innerHTML = html; t.classList.add('on');
    const r = t.getBoundingClientRect();
    t.style.left = Math.min(ev.clientX + 14, innerWidth - r.width - 8) + 'px';
    t.style.top = Math.min(ev.clientY + 14, innerHeight - r.height - 8) + 'px';
  },
  hide() { $('#tip').classList.remove('on'); },
};
const kv = pairs => '<div class="kv">' +
  pairs.filter(p => p[1] !== undefined && p[1] !== null && p[1] !== '')
       .map(p => `<span>${esc(p[0])}</span><span>${p[1]}</span>`).join('') + '</div>';

/* ══ THE CALL, re-applied client-side ═════════════════════════════════════════════════════ */
function called(r) {
  return r[DIST] >= S.dist && r[COV] >= S.cov && r[SAT] >= S.sat && !r[CONG];
}
const whatIf = () => S.cov !== S.meta.floors.cov || S.sat !== S.meta.floors.sat
                                                || S.dist !== S.meta.floors.distinct;

/* ══ MAP ══════════════════════════════════════════════════════════════════════════════════ */
const MAP_W = 1000, LAT_N = 84, LAT_S = -58;
const MAP_H = Math.round(MAP_W * (LAT_N - LAT_S) / 360);
const mx = lon => (lon + 180) / 360 * MAP_W;
const my = lat => (LAT_N - lat) / (LAT_N - LAT_S) * MAP_H;

function siteBurden(s) {           /* undeclared species per classified library at this site */
  return s.n_classified ? s.n_undec / s.n_classified : null;
}

function buildMap() {
  const stage = $('#map-stage');
  stage.innerHTML = '';

  const head = el('div');
  head.append(el('h2', 't', 'Where the cohort was collected'));
  head.append(el('p', 's', `${n0(S.sites.length)} distinct collection sites across `
    + `${n0(S.meta.counts.cohort_runs)} curated runs. Marker area is libraries at the site; `
    + `placement is the curated coordinate, never the submitting institution.`));
  stage.append(head);

  const wrap = el('div'); wrap.id = 'mapwrap';
  wrap.innerHTML = `<svg id="mapsvg" viewBox="0 0 ${MAP_W} ${MAP_H}" role="img"
      aria-label="World map of collection sites"><g id="mapg"></g></svg>
    <div class="mapctl"><button id="zin" title="Zoom in">+</button>
      <button id="zout" title="Zoom out">&minus;</button>
      <button id="zrst" title="Reset">&#8634;</button></div>`;
  stage.append(wrap);

  const g = wrap.querySelector('#mapg');
  const coast = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  coast.setAttribute('class', 'coast');
  coast.setAttribute('d', S.world.map(ring =>
    'M' + ring.map(p => `${mx(p[0]).toFixed(1)},${my(p[1]).toFixed(1)}`).join('L') + 'Z').join(''));
  g.append(coast);

  const maxN = Math.max(...S.sites.map(s => s.n_runs));
  const burdens = S.sites.map(siteBurden).filter(v => v != null).sort((a, b) => a - b);
  const brk = [0.2, 0.4, 0.6, 0.8].map(q => burdens[Math.floor(q * burdens.length)] || 0);
  const rampOf = v => v == null ? '#dedede' : RAMP[brk.filter(b => v >= b).length];

  const marks = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  marks.id = 'marks'; g.append(marks);
  S.sites.forEach((s, i) => {
    const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    c.setAttribute('class', 'site'); c.dataset.i = i;
    c.setAttribute('cx', mx(s.lon).toFixed(2)); c.setAttribute('cy', my(s.lat).toFixed(2));
    c.dataset.r = (2.6 + 12 * Math.sqrt(s.n_runs / maxN)).toFixed(2);
    marks.append(c);
  });
  paintMap(rampOf);
  wireMap(rampOf);

  /* keys: a size key for the magnitude channel, a ramp key only in burden mode */
  const keys = el('div');
  keys.innerHTML = `<div class="sizekey" id="sizekey"></div><div class="ramp" id="rampkey"></div>`;
  stage.append(keys);
  const sk = keys.querySelector('#sizekey');
  const demo = [10, 200, maxN];
  sk.innerHTML = '<span>Libraries</span>' + demo.map(n => {
    const r = 2.6 + 12 * Math.sqrt(n / maxN);
    return `<svg width="${2 * r + 4}" height="${2 * r + 4}"><circle cx="${r + 2}" cy="${r + 2}"
      r="${r}" fill="${SITE_COL}" fill-opacity=".5" stroke="#fff"/></svg><span>${n0(n)}</span>`;
  }).join('');
  keys.querySelector('#rampkey').innerHTML = '<span>Undeclared species per library</span>'
    + '<span class="steps">' + RAMP.map(c => `<i style="background:${c}"></i>`).join('')
    + '</span><span>low to high</span>';
  syncMapKeys();

  stage.append(mapTable());
}

function syncMapKeys() {
  const b = S.view === 'map' && S.mapMode === 'burden';
  const rk = $('#rampkey'); if (rk) rk.style.display = b ? 'flex' : 'none';
}

function paintMap(rampOf) {
  const crop = $('#m-crop') ? $('#m-crop').value : '';
  const k = S.map.k;
  document.querySelectorAll('#marks .site').forEach(c => {
    const s = S.sites[+c.dataset.i];
    const hit = !crop || s.crop.some(x => x[0] === crop);
    c.style.display = hit ? '' : 'none';
    if (!hit) return;
    c.setAttribute('r', (+c.dataset.r / k).toFixed(2));
    c.setAttribute('stroke-width', (1 / k).toFixed(2));
    if (S.mapMode === 'burden') {
      c.setAttribute('fill', rampOf(siteBurden(s)));
      c.setAttribute('fill-opacity', '.85');
    } else {
      c.setAttribute('fill', SITE_COL);
      c.setAttribute('fill-opacity', '.5');
    }
    c.classList.toggle('sel', S.sel.site === +c.dataset.i);
  });
}

function wireMap(rampOf) {
  const svg = $('#mapsvg'), g = $('#mapg');
  const apply = () => {
    g.setAttribute('transform', `translate(${S.map.tx},${S.map.ty}) scale(${S.map.k})`);
    paintMap(rampOf);
  };
  const zoomAt = (f, px, py) => {
    const k2 = Math.max(1, Math.min(14, S.map.k * f));
    const r = k2 / S.map.k;
    S.map.tx = px - (px - S.map.tx) * r;
    S.map.ty = py - (py - S.map.ty) * r;
    S.map.k = k2;
    /* keep the map inside its frame: no panning off into white space */
    S.map.tx = Math.min(0, Math.max(MAP_W * (1 - S.map.k), S.map.tx));
    S.map.ty = Math.min(0, Math.max(MAP_H * (1 - S.map.k), S.map.ty));
    apply();
  };
  const toLocal = ev => {
    const r = svg.getBoundingClientRect();
    return [(ev.clientX - r.left) / r.width * MAP_W, (ev.clientY - r.top) / r.height * MAP_H];
  };
  svg.addEventListener('wheel', ev => {
    ev.preventDefault();
    const [px, py] = toLocal(ev);
    zoomAt(ev.deltaY < 0 ? 1.18 : 1 / 1.18, px, py);
  }, {passive:false});
  $('#zin').onclick = () => zoomAt(1.5, MAP_W / 2, MAP_H / 2);
  $('#zout').onclick = () => zoomAt(1 / 1.5, MAP_W / 2, MAP_H / 2);
  $('#zrst').onclick = () => { S.map = {k:1, tx:0, ty:0}; apply(); };

  let drag = null;
  svg.addEventListener('pointerdown', ev => {
    drag = {x:ev.clientX, y:ev.clientY, tx:S.map.tx, ty:S.map.ty};
    svg.classList.add('drag'); svg.setPointerCapture(ev.pointerId);
  });
  svg.addEventListener('pointermove', ev => {
    if (!drag) return;
    const r = svg.getBoundingClientRect(), f = MAP_W / r.width;
    S.map.tx = Math.min(0, Math.max(MAP_W * (1 - S.map.k), drag.tx + (ev.clientX - drag.x) * f));
    S.map.ty = Math.min(0, Math.max(MAP_H * (1 - S.map.k), drag.ty + (ev.clientY - drag.y) * f));
    apply();
  });
  svg.addEventListener('pointerup', () => { drag = null; svg.classList.remove('drag'); });

  svg.addEventListener('mousemove', ev => {
    const t = ev.target.closest('.site');
    if (!t) return tip.hide();
    const s = S.sites[+t.dataset.i], b = siteBurden(s);
    tip.show(`<b>${esc(s.loc || s.country)}</b>` + kv([
      ['libraries', n0(s.n_runs)], ['classified', n0(s.n_classified)],
      ['projects', n0(s.n_proj)],
      ['crops', s.crop.map(c => c[0]).join(', ')],
      ['years', s.year.length ? s.year[0][0] + (s.year.length > 1
        ? '-' + s.year[s.year.length - 1][0] : '') : ''],
      ['co-infected', s.n_classified ? n0(s.n_coinf) + '/' + n0(s.n_classified) : ''],
      ['undeclared/lib', b == null ? '' : b.toFixed(2)],
      ['precision', s.prec],
    ]) + '<div style="margin-top:5px;color:#b9b9b2">click for the species list</div>', ev);
  });
  svg.addEventListener('mouseleave', tip.hide);
  svg.addEventListener('click', ev => {
    const t = ev.target.closest('.site');
    if (!t) return;
    S.sel.site = +t.dataset.i;
    paintMap(rampOf);
    showSite(S.sites[S.sel.site]);
  });
}

function showSite(s) {
  const box = $('#site-detail');
  const b = siteBurden(s);
  box.innerHTML = `<div class="lbl">Selected site</div>
    <div style="font-weight:650;margin:0 0 2px">${esc(s.loc || s.country)}</div>
    <div class="note" style="margin:0 0 9px">${esc(s.prec)} &middot; ${s.lat}, ${s.lon}</div>
    <div class="stats" style="margin:0 0 10px">
      <div class="stat"><div class="k">Libraries</div><div class="v">${n0(s.n_runs)}</div></div>
      <div class="stat"><div class="k">Classified</div><div class="v">${n0(s.n_classified)}</div></div>
      <div class="stat"><div class="k">Undecl/lib</div><div class="v">${b == null ? '--' : b.toFixed(2)}</div></div>
    </div>
    <div class="lbl">Species detected here</div>` + (s.sp.length
      ? '<div class="ranked">' + s.sp.map(([sp, n]) => {
          const w = Math.round(100 * n / s.sp[0][1]);
          return `<div class="row" style="grid-template-columns:1fr 54px 40px">
            <div class="nm">${ital(sp)}</div>
            <div class="track"><div class="fill" style="width:${w}%;background:${SITE_COL}"></div></div>
            <div class="num">${n0(n)}</div></div>`;
        }).join('') + '</div>'
      : '<div class="missing">No classified libraries from this site yet.</div>');
}

function mapTable() {
  const d = el('details', 'tbl');
  d.innerHTML = `<summary>Table view &mdash; all ${n0(S.sites.length)} sites</summary>
    <div class="scrollx"><table><thead><tr><th>Site</th><th>Country</th><th>Precision</th>
    <th>Libraries</th><th>Classified</th><th>Projects</th><th>Undecl/lib</th><th>Crops</th>
    </tr></thead><tbody>` + S.sites.map(s => {
      const b = siteBurden(s);
      return `<tr><td>${esc(s.loc)}</td><td>${esc(s.country)}</td><td>${esc(s.prec)}</td>
        <td class="n">${n0(s.n_runs)}</td><td class="n">${n0(s.n_classified)}</td>
        <td class="n">${n0(s.n_proj)}</td><td class="n">${b == null ? '--' : b.toFixed(2)}</td>
        <td>${esc(s.crop.map(c => c[0]).join(', '))}</td></tr>`;
    }).join('') + '</tbody></table></div>';
  return d;
}

/* ══ SPECIES ══════════════════════════════════════════════════════════════════════════════ */
function speciesList() {
  const q = ($('#sp-q').value || '').toLowerCase();
  const only = $('#sp-path').checked;
  const king = $('#sp-king').value;
  const ul = $('#sp-pick'); ul.innerHTML = '';
  /* Fungi by default. The host plant appears in detections.tsv because detect.py reports it --
     a reader checks which plant the library came from before reading anything else -- but it
     is not a finding, and ranking by libraries would otherwise open this view on wheat. */
  const hits = S.species.filter(s => (!q || s.sp.toLowerCase().includes(q))
    && (!only || s.pathogen) && (!king || s.kingdom === king));
  $('#sp-count').textContent = `${n0(hits.length)} of ${n0(S.species.length)} species`;
  hits.forEach(s => {
      const li = el('li');
      li.innerHTML = `<i>${esc(s.sp)}</i><span class="n">${n0(s.n_runs)}</span>`;
      li.setAttribute('aria-selected', S.sel.species === s.sp);
    li.onclick = () => { S.sel.species = s.sp; speciesList(); showSpecies(s); };
    ul.append(li);
  });
}

function barRows(pairs, total, colour) {
  if (!pairs.length) return '<div class="missing">No record.</div>';
  const max = pairs[0][1];
  return '<div class="ranked">' + pairs.slice(0, 14).map(([k, n]) =>
    `<div class="row" style="grid-template-columns:minmax(90px,1fr) minmax(70px,1.4fr) 50px 48px">
      <div class="nm">${esc(k)}</div>
      <div class="track"><div class="fill" style="width:${Math.round(100 * n / max)}%;
        background:${colour}"></div></div>
      <div class="num">${n0(n)}</div>
      <div class="num dim">${total ? (100 * n / total).toFixed(0) + '%' : ''}</div>
    </div>`).join('') + '</div>'
    + (pairs.length > 14 ? `<div class="note">+ ${pairs.length - 14} more, in the table view</div>` : '');
}

function showSpecies(s) {
  const st = $('#sp-stage');
  const q = s.cov_q, undeclPct = s.n_runs ? 100 * s.n_undeclared / s.n_runs : 0;
  st.innerHTML = `
    <h2 class="t"><i>${esc(s.sp)}</i></h2>
    <p class="s">taxid ${esc(s.taxid)} &middot; ${esc(s.kingdom)}
      ${s.pathogen ? '&middot; <b>pathogen</b>' : ''} &middot; membership route
      <span style="font-family:var(--mono)">${esc(s.route)}</span></p>
    ${s.guild ? `<p class="note">FUNGuild: ${esc(s.guild)}</p>` : ''}

    <div class="stats">
      <div class="stat"><div class="k">Libraries</div><div class="v">${n0(s.n_runs)}</div></div>
      <div class="stat"><div class="k">BioProjects</div><div class="v">${n0(s.n_proj)}</div></div>
      <div class="stat"><div class="k">Undeclared</div><div class="v">${n0(s.n_undeclared)}
        <small>${undeclPct.toFixed(0)}%</small></div></div>
      <div class="stat"><div class="k">Declared</div><div class="v">${n0(s.n_declared)}</div></div>
      <div class="stat"><div class="k">db clade k-mers</div><div class="v">${n0(s.clade)}</div></div>
    </div>

    <h3 style="font-size:13px;margin:22px 0 2px">Coverage where detected</h3>
    <p class="note" style="margin:0 0 8px">Fraction of this species' k-mer space seen, log scale.
      Box is the interquartile range, the rule is the median, whiskers reach the minimum and the
      99th percentile. The dashed rule is the ${(S.meta.floors.cov * 100).toFixed(0)}% floor every
      one of these cleared.</p>
    <div id="sp-cov"></div>

    <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:22px;margin-top:22px">
      <div><h3 style="font-size:13px;margin:0 0 6px">Crop</h3>
        ${barRows(s.crop, s.n_runs, SITE_COL)}</div>
      <div><h3 style="font-size:13px;margin:0 0 6px">Country</h3>
        ${barRows(s.country, s.n_runs, SITE_COL)}</div>
    </div>

    <h3 style="font-size:13px;margin:22px 0 6px">Collection year</h3>
    ${barRows(s.year, s.n_runs, SITE_COL)}

    <details class="tbl" open><summary>Libraries detecting this species
      (${n0(Math.min(s.runs.length, 400))}${s.n_runs > s.runs.length
        ? ' of ' + n0(s.n_runs) + ', capped' : ''})</summary>
      <div class="scrollx"><p class="note">Click a run to open it in the run view.</p>
      ${s.runs.map(r => `<button class="runlink" data-run="${esc(r)}"
        style="font:12px var(--mono);border:1px solid var(--line);background:#fff;border-radius:3px;
        padding:2px 6px;margin:0 3px 3px 0;cursor:pointer">${esc(r)}</button>`).join('')}</div>
    </details>`;

  covStrip($('#sp-cov'), q);
  st.querySelectorAll('.runlink').forEach(b =>
    b.onclick = () => { location.hash = '#/runs/' + b.dataset.run; });
}

/* five-number summary on one log axis. Not a chart type for its own sake: the question is
   "when this is found, how strongly", and a distribution answers it where a mean would not. */
function covStrip(host, q) {
  if (!q || q[0] == null) { host.innerHTML = '<div class="missing">No coverage recorded.</div>'; return; }
  const W = 560, H = 44, PAD = 8, lo = 1e-3, hi = 1;
  const x = logScale(lo, hi, W - 2 * PAD);
  const X = v => PAD + x(Math.max(lo, Math.min(hi, v)));
  const [mn, p25, p50, p75, p99] = q;
  const ticks = [1e-3, 1e-2, 1e-1, 1];
  host.innerHTML = `<svg viewBox="0 0 ${W} ${H}" style="max-width:${W}px;width:100%;height:auto">
    ${ticks.map(t => `<line class="gridline" x1="${X(t)}" x2="${X(t)}" y1="4" y2="24"/>`).join('')}
    <line class="axis" x1="${PAD}" x2="${W - PAD}" y1="24" y2="24"/>
    <line class="floorline" x1="${X(S.meta.floors.cov)}" x2="${X(S.meta.floors.cov)}" y1="2" y2="26"/>
    <line x1="${X(mn)}" x2="${X(p99)}" y1="14" y2="14" stroke="${TIER_COL.pathogen_other}"
      stroke-width="2"/>
    <rect x="${X(p25)}" y="7" width="${Math.max(2, X(p75) - X(p25))}" height="14" rx="3"
      fill="${TIER_COL.pathogen_other}" fill-opacity=".45" stroke="${TIER_COL.pathogen_crop}"
      stroke-width="0"/>
    <line x1="${X(p50)}" x2="${X(p50)}" y1="5" y2="23" stroke="${TIER_COL.pathogen_crop}"
      stroke-width="2"/>
    ${ticks.map(t => `<text class="axlab" x="${X(t)}" y="38" text-anchor="middle">${
      t >= 1 ? '100%' : (t * 100) + '%'}</text>`).join('')}
    <text class="axlab" x="${X(p50)}" y="4" text-anchor="middle">${pct(p50)}</text>
  </svg>`;
}

/* ══ RUNS ═════════════════════════════════════════════════════════════════════════════════ */
function runList() {
  const q = ($('#r-q').value || '').toLowerCase();
  const crop = $('#r-crop').value, proj = $('#r-proj').value;
  const ul = $('#r-pick'); ul.innerHTML = '';
  const hits = S.runs.filter(r =>
    (!crop || r.crop === crop) && (!proj || r.p === proj) &&
    (!q || r.r.toLowerCase().includes(q) || (r.loc || '').toLowerCase().includes(q)));
  $('#r-count').textContent = `${n0(hits.length)} of ${n0(S.runs.length)} runs`;
  hits.slice(0, 500).forEach(r => {
    const li = el('li');
    const n = r.d.filter(called).length;
    li.innerHTML = `<span style="font-family:var(--mono);font-size:12px">${esc(r.r)}</span>
      <span class="n">${n}</span>`;
    li.setAttribute('aria-selected', S.sel.run === r.r);
    li.onclick = () => { location.hash = '#/runs/' + r.r; };
    ul.append(li);
  });
  if (hits.length > 500) ul.append(el('li', '', `<span class="n">+ ${n0(hits.length - 500)} more; narrow the filter</span>`));
}

function showRun(run) {
  S.sel.run = run.r;
  const st = $('#r-stage');
  const calls = run.d.filter(called);
  const undecl = calls.filter(r => !r[DECL] && r[TIER] !== 0);

  st.innerHTML = `
    <div style="display:flex;align-items:flex-start;gap:14px;flex-wrap:wrap">
      <div style="min-width:0">
        <h2 class="t"><button id="r-title" style="all:unset;cursor:pointer;
          border-bottom:1px dashed var(--ink-3)">${esc(run.r)}</button></h2>
        <p class="s">${esc(run.p)} &middot; ${esc(run.host)} &middot; ${esc(run.tis)}<br>
          ${esc(run.loc || 'location not curated')}${run.year ? ' &middot; ' + esc(run.year) : ''}
          ${run.prec ? ' &middot; <span style="font-family:var(--mono);font-size:11.5px">'
            + esc(run.prec) + '</span>' : ''}</p>
        <p class="note">Click the run id for the curated evidence behind every one of those values.</p>
      </div>
      <div style="margin-left:auto;display:flex;flex-direction:column;gap:5px;flex:0 0 auto">
        <div class="seg" id="r-mode" style="width:196px">
          <button data-m="rows" aria-pressed="${S.runMode === 'rows'}">Ranked</button>
          <button data-m="scatter" aria-pressed="${S.runMode === 'scatter'}">Floor space</button>
        </div>
        <div class="seg" id="r-depth" style="width:196px;${
          S.runMode === 'rows' ? '' : 'visibility:hidden'}">
          <button data-a="0" aria-pressed="${!S.rankedAll}">Called</button>
          <button data-a="1" aria-pressed="${S.rankedAll}">All candidates</button>
        </div>
      </div>
    </div>
    <div class="whatif ${whatIf() ? 'on' : ''}">WHAT-IF: a floor has been moved off the
      pipeline's value. These calls are not pipeline results.</div>

    <div class="stats">
      <div class="stat"><div class="k">Called</div><div class="v">${n0(calls.length)}</div></div>
      <div class="stat"><div class="k">Undeclared</div><div class="v">${n0(undecl.length)}</div></div>
      <div class="stat"><div class="k">Declared by study</div>
        <div class="v" style="font-size:13px;font-weight:500;line-height:1.35">${
          run.dp ? ital(run.dp) : '<span style="color:var(--ink-3)">none</span>'}</div></div>
    </div>

    <div class="comp">
      <div class="lbl">What the library is made of &mdash; classified read mass by tier</div>
      <div class="bar" id="r-comp"></div>
      <div class="note" id="r-comp-note"></div>
    </div>

    <div id="r-chart"></div>`;

  compStrip(run);
  drawRunChart(run);
  $('#r-title').onclick = () => openEvidence(run);
  st.querySelectorAll('#r-mode button').forEach(b => b.onclick = () => {
    S.runMode = b.dataset.m; showRun(run);
  });
  st.querySelectorAll('#r-depth button').forEach(b => b.onclick = () => {
    S.rankedAll = b.dataset.a === '1'; showRun(run);
  });
}

function compStrip(run) {
  const tiers = S.meta.tiers, total = run.mass.reduce((a, b) => a + b, 0);
  const bar = $('#r-comp');
  if (!total) { bar.innerHTML = ''; $('#r-comp-note').textContent = 'No classified reads.'; return; }
  bar.innerHTML = tiers.map((t, i) => run.mass[i]
    ? `<i style="width:${(100 * run.mass[i] / total).toFixed(2)}%;background:${TIER_COL[t]}"
         data-t="${t}" data-n="${run.mass[i]}"></i>` : '').join('');
  $('#r-comp-note').textContent = n0(total) + ' classified reads. '
    + tiers.map((t, i) => run.mass[i] ? `${TIER_HEAD[t].toLowerCase()} `
        + (100 * run.mass[i] / total).toFixed(1) + '%' : null).filter(Boolean).join(', ') + '.';
  bar.onmousemove = ev => {
    const t = ev.target.closest('i');
    if (!t) return tip.hide();
    tip.show(`<b>${esc(TIER_LAB[t.dataset.t])}</b>` + kv([
      ['reads', n0(+t.dataset.n)],
      ['share', (100 * t.dataset.n / total).toFixed(2) + '%']]), ev);
  };
  bar.onmouseleave = tip.hide;
}

function drawRunChart(run) {
  const host = $('#r-chart');
  host.innerHTML = legendHTML()
    + (S.runMode === 'rows' ? rankedHTML(run) : '<div class="scatter" id="sc"></div>')
    + runTable(run);
  if (S.runMode === 'scatter') scatter($('#sc'), run);
  else wireRanked(host, run);
}

function legendHTML() {
  return '<div class="legend">' + S.meta.tiers.map(t =>
    `<span><i class="sw" style="background:${TIER_COL[t]}"></i>${esc(TIER_LAB[t])}</span>`).join('')
    + '<span><i class="sw" style="background:none;border:1.6px solid var(--ink)"></i>'
    + 'declared by the study</span></div>';
}

/* ── the ranked rows. The default, because the browsing question is "what is in this sample
      and should I believe it", which is a ranked list with a pass mark, not a 2-D space. ──── */
const COVLO = 1e-5, SATLO = 1e-3;
function rankedHTML(run) {
  const rows = run.d.slice().sort((a, b) => b[COV] - a[COV]);
  const byTier = new Map();
  rows.forEach(r => {
    const t = S.meta.tiers[r[TIER]];
    if (!byTier.has(t)) byTier.set(t, []);
    byTier.get(t).push(r);
  });
  const xc = logScale(COVLO, 1, 100), xs = logScale(SATLO, 10, 100);

  let h = `<div class="ranked"><div class="hd"><div></div><div>Species</div>
    <div>Coverage</div><div style="text-align:right">value</div><div>Saturation</div>
    <div style="text-align:right">value</div><div style="text-align:right">Distinct</div></div>`;

  for (const t of S.meta.tiers) {
    const rs = byTier.get(t);
    if (!rs || !rs.length) continue;
    const nc = rs.filter(called).length;
    const show = S.rankedAll ? rs : rs.filter(called);
    if (!show.length) continue;
    h += `<div class="grp"><i class="sw" style="background:${TIER_COL[t]}"></i>
      ${esc(TIER_HEAD[t])}<span class="c">${nc} called of ${rs.length}</span></div>`;
    h += show.map((r, i) => {
      const ok = called(r), col = TIER_COL[t];
      const sp = S.meta.trophic[r[TROPH]];
      return `<div class="row ${ok ? '' : 'out'}" data-sp="${esc(r[SP])}"
          data-k="${S.meta.tiers.indexOf(t)}-${i}">
        <div class="dec">${r[DECL] ? '&#9679;' : ''}</div>
        <div class="nm">${ital(r[SP])}${sp ? `<em>${esc(sp)}</em>` : ''}${
          r[CONG] ? '<em>congener only</em>' : ''}</div>
        <div class="track"><div class="fill" style="width:${xc(r[COV]).toFixed(1)}%;
            background:${col};opacity:${ok ? 1 : .5}"></div>
          <div class="fl" style="left:${xc(S.cov).toFixed(1)}%"></div></div>
        <div class="num ${ok ? '' : 'dim'}">${pct(r[COV])}</div>
        <div class="track"><div class="fill" style="width:${xs(r[SAT]).toFixed(1)}%;
            background:${col};opacity:${ok ? 1 : .5}"></div>
          <div class="fl" style="left:${xs(S.sat).toFixed(1)}%"></div></div>
        <div class="num ${ok ? '' : 'dim'}">${sig(r[SAT], 2)}</div>
        <div class="num ${r[DIST] >= S.dist ? '' : 'dim'}">${n0(r[DIST])}</div>
      </div>`;
    }).join('');
  }
  h += '</div>';
  const nShown = S.rankedAll ? rows.length : rows.filter(called).length;
  if (!nShown) h = h.replace('</div>',
    '<div class="missing" style="margin:10px 0">Nothing was called in this library at the '
    + 'floors in force. Switch to <b>All candidates</b> to see what was rejected and by '
    + 'which floor.</div></div>');
  /* axis captions under the two tracks, so the log scale is stated rather than implied */
  h += `<div class="axcap"><div></div><div></div>
    <div class="r">${[1e-5, 1e-3, 1e-1, 1].map(v =>
      `<span style="left:${logScale(COVLO, 1, 100)(v)}%">${v >= 1 ? '100%' : pct(v)}</span>`).join('')}</div>
    <div></div>
    <div class="r">${[1e-3, 1e-1, 1, 10].map(v =>
      `<span style="left:${logScale(SATLO, 10, 100)(v)}%">${v}</span>`).join('')}</div>
    <div></div><div></div></div>
    <p class="note">Both tracks are log scales; the dashed rule is the floor in force.
      ${S.rankedAll ? '' : `Showing the ${n0(nShown)} called rows; each heading gives the
      group's total, and <b>All candidates</b> shows the rest with the floor each one missed. `}
      A further ${n0(run.kn)} candidates fell below the ${n0(S.meta.floors.distinct)}-k-mer
      absolute floor and are never listed &mdash; they carry no identity worth reading, and
      their shape is the faint cloud in the floor-space view.</p>`;
  return h;
}

function wireRanked(host, run) {
  const idx = new Map();
  run.d.forEach(r => idx.set(r[SP], r));
  host.querySelectorAll('.ranked .row[data-sp]').forEach(row => {
    row.onmousemove = ev => {
      const r = idx.get(row.dataset.sp);
      if (!r) return;
      const ok = called(r);
      tip.show(`<b><i>${esc(r[SP])}</i></b>` + kv([
        ['coverage', pct(r[COV])], ['saturation', sig(r[SAT], 3)],
        ['redundancy', sig(r[RED], 3)], ['distinct k-mers', n0(r[DIST])],
        ['reads', n0(r[READS])],
        ['tier', TIER_LAB[S.meta.tiers[r[TIER]]]],
        ['trophic', S.meta.trophic[r[TROPH]]],
        ['evidence', S.meta.evid[r[EVID]]],
        ['declared', r[DECL] ? 'yes' : 'no'],
        ['call', ok ? 'YES' : (r[CONG] ? 'no (congener only)'
          : r[DIST] < S.dist ? 'no (absolute floor)'
          : r[COV] < S.cov ? 'no (coverage floor)' : 'no (saturation floor)')],
      ]), ev);
    };
    row.onmouseleave = tip.hide;
  });
}

/* ── the floor-space scatter: kept, stripped. Colour carries tier and nothing else; trophic
      moved to the row list and the tooltip, labels moved to hover. ──────────────────────── */
function scatter(host, run) {
  const W = 620, H = 430, L = 54, R = 12, T = 12, B = 44;
  const xs = [1e-4, 1e4], ys = [1e-5, 1];
  const px = logScale(xs[0], xs[1], W - L - R), py = logScale(ys[0], ys[1], H - T - B);
  const X = v => L + px(Math.max(xs[0], Math.min(xs[1], v)));
  const Y = v => H - B - py(Math.max(ys[0], Math.min(ys[1], v)));
  const xt = [1e-4, 1e-2, 1, 1e2, 1e4], yt = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];

  /* the floor curve: coverage = sat_floor x (1 - exp(-redundancy)). The shape test IS this
     curve, so drawing it is the only way the rejections read as a rule rather than a scatter. */
  const curve = [];
  for (let i = 0; i <= 140; i++) {
    const red = Math.pow(10, Math.log10(xs[0]) + i / 140 * (Math.log10(xs[1]) - Math.log10(xs[0])));
    curve.push(`${X(red).toFixed(1)},${Y(Math.max(ys[0], S.sat * (1 - Math.exp(-red)))).toFixed(1)}`);
  }

  const pts = run.d.slice().sort((a, b) => a[COV] - b[COV]).map(r => {
    const ok = called(r), col = TIER_COL[S.meta.tiers[r[TIER]]];
    return `<circle class="mk${r[DECL] ? '' : ''}" cx="${X(r[RED]).toFixed(1)}"
      cy="${Y(r[COV]).toFixed(1)}" r="${r[DECL] ? 7 : 5}"
      fill="${ok ? col : '#fff'}" stroke="${ok ? 'var(--surface)' : col}"
      stroke-width="${ok ? 2 : 1.6}" data-sp="${esc(r[SP])}"/>`
      + (r[DECL] ? `<circle class="mk hollow" cx="${X(r[RED]).toFixed(1)}"
          cy="${Y(r[COV]).toFixed(1)}" r="10.5" stroke="${col}"/>` : '');
  }).join('');

  host.innerHTML = `<svg viewBox="0 0 ${W} ${H}" role="img"
      aria-label="Coverage against redundancy, with the detection floors">
    ${yt.map(t => `<line class="gridline" x1="${L}" x2="${W - R}" y1="${Y(t)}" y2="${Y(t)}"/>`).join('')}
    ${xt.map(t => `<line class="gridline" x1="${X(t)}" x2="${X(t)}" y1="${T}" y2="${H - B}"/>`).join('')}
    <line class="axis" x1="${L}" x2="${W - R}" y1="${H - B}" y2="${H - B}"/>
    <line class="axis" x1="${L}" x2="${L}" y1="${T}" y2="${H - B}"/>
    ${run.k.map(p => `<circle class="cloudpt" cx="${X(p[0]).toFixed(1)}"
      cy="${Y(p[1]).toFixed(1)}" r="1.6"/>`).join('')}
    <polyline class="floorline" points="${curve.join(' ')}"/>
    <line class="floorline" x1="${L}" x2="${W - R}" y1="${Y(S.cov)}" y2="${Y(S.cov)}"/>
    ${pts}
    ${yt.map(t => `<text class="axlab" x="${L - 7}" y="${Y(t) + 4}" text-anchor="end">${
      t >= 1 ? '100%' : pct(t)}</text>`).join('')}
    ${xt.map(t => `<text class="axlab" x="${X(t)}" y="${H - B + 16}" text-anchor="middle">${
      t >= 1 ? n0(t) : sig(t, 1)}</text>`).join('')}
    <text class="axttl" x="${(L + W - R) / 2}" y="${H - 8}" text-anchor="middle">Redundancy</text>
    <text class="axttl" transform="translate(13,${(T + H - B) / 2}) rotate(-90)"
      text-anchor="middle">Coverage</text>
  </svg>
  <p class="note">Each point is one candidate species. The curved rule is the saturation floor,
    coverage = ${S.sat} &times; (1 &minus; e<sup>&minus;redundancy</sup>); the straight rule is the
    ${pct(S.cov)} coverage floor. Filled points were called, hollow were not; a ring marks a
    species the study itself declared. The faint cloud is ${n0(run.kn)} candidates below the
    absolute floor, downsampled to ${n0(S.meta.cloud_n)} for the shape.</p>`;

  const idx = new Map(run.d.map(r => [r[SP], r]));
  host.querySelector('svg').addEventListener('mousemove', ev => {
    const t = ev.target.closest('circle[data-sp]');
    if (!t) return tip.hide();
    const r = idx.get(t.dataset.sp);
    tip.show(`<b><i>${esc(r[SP])}</i></b>` + kv([
      ['coverage', pct(r[COV])], ['saturation', sig(r[SAT], 3)],
      ['redundancy', sig(r[RED], 3)], ['distinct', n0(r[DIST])],
      ['tier', TIER_LAB[S.meta.tiers[r[TIER]]]],
      ['declared', r[DECL] ? 'yes' : 'no'],
      ['call', called(r) ? 'YES' : 'no']]), ev);
  });
  host.querySelector('svg').addEventListener('mouseleave', tip.hide);
}

function runTable(run) {
  const rows = run.d.slice().sort((a, b) => b[COV] - a[COV]);
  return `<details class="tbl"><summary>Table view &mdash; ${n0(rows.length)} candidate rows</summary>
    <div class="scrollx"><table><thead><tr><th>Species</th><th>Tier</th><th>Trophic</th>
    <th>Declared</th><th>Coverage</th><th>Saturation</th><th>Redundancy</th><th>Distinct</th>
    <th>Reads</th><th>Call</th></tr></thead><tbody>` + rows.map(r =>
      `<tr><td><i>${esc(r[SP])}</i></td><td>${esc(TIER_HEAD[S.meta.tiers[r[TIER]]])}</td>
       <td>${esc(S.meta.trophic[r[TROPH]])}</td><td>${r[DECL] ? 'yes' : ''}</td>
       <td class="n">${pct(r[COV])}</td><td class="n">${sig(r[SAT], 3)}</td>
       <td class="n">${sig(r[RED], 3)}</td><td class="n">${n0(r[DIST])}</td>
       <td class="n">${n0(r[READS])}</td><td>${called(r) ? 'yes' : ''}</td></tr>`).join('')
    + '</tbody></table></div></details>';
}

/* ══ EVIDENCE DRAWER ══════════════════════════════════════════════════════════════════════ */
const FIELD_LAB = {setting:'Setting', selection:'Sampling selection', tissue:'Tissue',
                   location:'Location', year:'Collection year'};

function openEvidence(run) {
  const ev = S.evidence;
  const parts = [];
  const vals = {setting:run.set, selection:run.sel, tissue:run.tis,
                location:run.loc, year:run.year};

  for (const [field, rule] of Object.entries(run.rules || {})) {
    const claims = ev.by_rule[rule] || [];
    const label = FIELD_LAB[field] || field;
    if (!claims.length) {
      parts.push(`<div class="claim"><div class="f">${esc(label)}</div>
        <div class="v">${esc(vals[field] || '--')}</div>
        <div class="missing">No quotation is registered for rule
          <span style="font-family:var(--mono)">${esc(rule)}</span>. The value comes from the
          archive record rather than the manuscript.</div></div>`);
      continue;
    }
    claims.forEach(c => parts.push(`<div class="claim">
      <div class="f">${esc(label)} &middot; rule ${esc(rule)}</div>
      <div class="v">${esc(vals[field] || c.field)}</div>
      <blockquote>${esc(c.quote)}</blockquote>
      <div class="src">${c.doi ? `<a href="https://doi.org/${encodeURIComponent(c.doi)}"
        target="_blank" rel="noopener">${esc(c.doi)}</a>` : 'no DOI recorded'}
        ${c.locus ? ' &middot; ' + esc(c.locus) : ''}
        ${c.n ? ' &middot; asserted for ' + n0(c.n) + ' runs' : ''}</div></div>`));
  }

  const proj = (ev.by_project[run.p] || []).filter(c =>
    !Object.values(run.rules || {}).includes(c.id));
  $('#drawer .dh h3').innerHTML = `Evidence for ${esc(run.r)}`;
  $('#drawer .db').innerHTML = `
    <p class="note" style="margin:0 0 12px">Every curated value below is backed by a quotation
      from the study that produced the sample. ${esc(run.p)} &middot;
      ${n0((ev.by_project[run.p] || []).length)} registered claims for this BioProject.</p>
    ${parts.join('') || '<div class="missing">No curated rules are recorded for this run.</div>'}
    ${proj.length ? `<h4 style="font-size:12.5px;margin:16px 0 7px">Other claims for
      ${esc(run.p)}</h4>` + proj.map(c => `<div class="claim">
        <div class="f">${esc(c.field)}</div><blockquote>${esc(c.quote)}</blockquote>
        <div class="src">${c.doi ? esc(c.doi) : 'no DOI'}${c.locus ? ' &middot; '
          + esc(c.locus) : ''}</div></div>`).join('') : ''}`;
  $('#drawer').classList.add('on');
  $('#scrim').classList.add('on');
}
function closeEvidence() {
  $('#drawer').classList.remove('on');
  $('#scrim').classList.remove('on');
}

/* ══ ROUTER & BOOT ════════════════════════════════════════════════════════════════════════ */
function setView(v) {
  S.view = v;
  ['map', 'species', 'runs'].forEach(x => {
    $('#v-' + x).classList.toggle('on', x === v);
    $('#nav-' + x).setAttribute('aria-selected', x === v);
  });
  syncMapKeys();
}

async function route() {
  const h = location.hash.replace(/^#\/?/, '');
  const [v, arg] = h.split('/');
  if (v === 'species') {
    setView('species');
    if (arg) {
      const s = S.species.find(x => x.sp === decodeURIComponent(arg));
      if (s) { S.sel.species = s.sp; speciesList(); showSpecies(s); }
    }
  } else if (v === 'runs') {
    setView('runs');
    await ensureRuns();
    if (arg) {
      const r = S.runs.find(x => x.r === arg);
      if (r) { showRun(r); runList(); }
    }
  } else {
    setView('map');
  }
}

let runsP = null;
function ensureRuns() {
  if (S.runs) return Promise.resolve();
  if (runsP) return runsP;
  $('#r-stage').innerHTML = '<p class="s">Loading the run panels (11 MB, once per visit)...</p>';
  runsP = fetch('data/runs.json').then(r => r.json()).then(j => {
    S.runs = j.runs;
    const projs = [...new Set(S.runs.map(r => r.p))].sort();
    $('#r-proj').innerHTML = '<option value="">all projects</option>'
      + projs.map(p => `<option>${p}</option>`).join('');
    const crops = [...new Set(S.runs.map(r => r.crop).filter(Boolean))].sort();
    $('#r-crop').innerHTML = '<option value="">all crops</option>'
      + crops.map(c => `<option>${c}</option>`).join('');
    runList();
    $('#r-stage').innerHTML = '<p class="s">Pick a run on the left.</p>';
  });
  return runsP;
}

function floors() {
  const set = () => {
    S.cov = +$('#f-cov').value / 10000;
    S.sat = +$('#f-sat').value / 1000;
    S.dist = +$('#f-dist').value;
    $('#f-cov-v').textContent = pct(S.cov);
    $('#f-sat-v').textContent = S.sat.toFixed(3);
    $('#f-dist-v').textContent = n0(S.dist);
    if (S.runs && S.sel.run) {
      const r = S.runs.find(x => x.r === S.sel.run);
      if (r) showRun(r);
      runList();
    }
  };
  ['f-cov', 'f-sat', 'f-dist'].forEach(id => $('#' + id).addEventListener('input', set));
  $('#f-reset').onclick = () => {
    $('#f-cov').value = S.meta.floors.cov * 10000;
    $('#f-sat').value = S.meta.floors.sat * 1000;
    $('#f-dist').value = S.meta.floors.distinct;
    set();
  };
  set();
}

async function boot() {
  const [meta, sites, species, evidence, world] = await Promise.all(
    ['meta.json', 'sites.json', 'species.json', 'evidence.json', 'world.json']
      .map(f => fetch('data/' + f).then(r => r.json())));
  Object.assign(S, {meta, sites, species, evidence, world});
  S.cov = meta.floors.cov; S.sat = meta.floors.sat; S.dist = meta.floors.distinct;

  $('#hdr-sub').textContent = `${n0(meta.counts.cohort_runs)} curated runs &middot; `
    .replace('&middot;', '·') + `${n0(meta.counts.projects)} BioProjects · `
    + `${n0(meta.counts.species)} species called · ${n0(meta.counts.sites)} sites`;

  ['map', 'species', 'runs'].forEach(v =>
    $('#nav-' + v).onclick = () => { location.hash = '#/' + v; });

  /* map controls */
  const crops = [...new Set(S.sites.flatMap(s => s.crop.map(c => c[0])))].sort();
  $('#m-crop').innerHTML = '<option value="">all crops</option>'
    + crops.map(c => `<option>${c}</option>`).join('');
  buildMap();
  $('#m-crop').onchange = () => { S.sel.site = null; $('#site-detail').innerHTML = ''; buildMap(); };
  document.querySelectorAll('#m-mode button').forEach(b => b.onclick = () => {
    S.mapMode = b.dataset.m;
    document.querySelectorAll('#m-mode button').forEach(x =>
      x.setAttribute('aria-pressed', x.dataset.m === S.mapMode));
    buildMap();
  });
  const un = Object.values(meta.unplaced).reduce((a, b) => a + b, 0);
  $('#m-unplaced').textContent = un
    ? `${n0(un)} runs have no curated coordinate and are not on the map.` : '';

  /* species controls */
  const kings = [...new Set(S.species.map(s => s.kingdom))].sort();
  $('#sp-king').innerHTML = kings.map(k =>
    `<option${k === 'Fungi' ? ' selected' : ''}>${k}</option>`).join('')
    + '<option value="">all kingdoms</option>';
  speciesList();
  ['sp-q', 'sp-king'].forEach(id => $('#' + id).addEventListener('input', speciesList));
  $('#sp-king').addEventListener('change', speciesList);
  $('#sp-path').addEventListener('change', speciesList);

  /* run controls */
  floors();
  ['r-q', 'r-crop', 'r-proj'].forEach(id =>
    $('#' + id).addEventListener('input', () => S.runs && runList()));

  $('#scrim').onclick = closeEvidence;
  $('#drawer .dh button').onclick = closeEvidence;
  addEventListener('keydown', e => { if (e.key === 'Escape') { closeEvidence(); tip.hide(); } });

  addEventListener('hashchange', route);
  route();
}
boot();
