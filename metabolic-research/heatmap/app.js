/* Longevity Biomarker Heat Map — vanilla JS, single view.
   Loads the bundled data layer (window.BIOMARKER_DATA) or fetches the JSON. */

'use strict';

const OUTCOMES = [
  { key: 'all_cause_mortality', label: 'All-Cause\nMortality', anchor: true },
  { key: 'cvd',                 label: 'Cardiovascular\nDisease' },
  { key: 'cancer',              label: 'Cancer' },
  { key: 'dementia',            label: 'Dementia /\nCognitive Decline' },
  { key: 'frailty',             label: 'Frailty /\nSarcopenia' }
];

const MOD_RATINGS = {
  high:     { label: 'High',     color: '#1a7a3c', rank: 1.00 },
  moderate: { label: 'Moderate', color: '#7cae4f', rank: 0.66 },
  low:      { label: 'Low',      color: '#d3cf8f', rank: 0.33 },
  fixed:    { label: 'Fixed',    color: '#bdbdbd', rank: 0.00 }
};

const LN2 = Math.log(2); // magnitude scale ceiling

const state = { sort: 'domain', filter: '', collapsed: new Set() };
let DATA = null;

/* ---- math / color helpers --------------------------------------------- */
function magnitude(hr) {
  if (hr == null || hr <= 0) return null;
  return Math.abs(Math.log(hr));
}
function clamp01(t) { return Math.max(0, Math.min(1, t)); }
function lerp(a, b, t) { return a + (b - a) * t; }

function ramp(stops, t) {
  t = clamp01(t);
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [p0, c0] = stops[i - 1], [p1, c1] = stops[i];
      const f = (t - p0) / (p1 - p0);
      const c = c0.map((v, j) => Math.round(lerp(v, c1[j], f)));
      return `rgb(${c[0]},${c[1]},${c[2]})`;
    }
  }
  const last = stops[stops.length - 1][1];
  return `rgb(${last[0]},${last[1]},${last[2]})`;
}
const RISK_RAMP = [[0, [255, 255, 255]], [0.5, [253, 187, 132]], [1, [178, 24, 43]]];
const PROT_RAMP = [[0, [255, 255, 255]], [0.5, [146, 197, 222]], [1, [33, 102, 172]]];

/* The comparable figure the magnitude scale uses: HR per +1 SD.
   null for categorical / unconvertible / qualitative / no-data cells. */
function comparableHr(cell) {
  return (cell && typeof cell.hr_per_sd === 'number') ? cell.hr_per_sd : null;
}
function effectColor(cell) {
  const h = comparableHr(cell);
  if (h == null) return null;
  const t = clamp01(magnitude(h) / LN2);
  return ramp(cell.direction === 'protective' ? PROT_RAMP : RISK_RAMP, t);
}
function textColorFor(bg) {
  const m = /rgb\((\d+),(\d+),(\d+)\)/.exec(bg);
  if (!m) return '#1c1c1c';
  const lum = (0.299 * +m[1] + 0.587 * +m[2] + 0.114 * +m[3]) / 255;
  return lum > 0.62 ? '#1c1c1c' : '#ffffff';
}
function tierClass(tier) {
  return tier === 'A' ? 'tier-A' : tier === 'B' ? 'tier-B' : 'tier-C';
}

/* Cell states:
   'comparable'    — has a per-SD HR; sits on the comparable magnitude scale
   'categorical'   — categorical exposure; native HR shown, off the scale
   'unconvertible' — continuous but not standardizable; native HR shown, flagged
   'qualitative'   — evidence documented but no usable pooled HR
   'none'          — no evidence at all */
function cellState(cell) {
  if (!cell) return 'none';
  if (typeof cell.hr_per_sd === 'number') return 'comparable';
  if (typeof cell.hr === 'number') {
    return cell.comparability === 'categorical' ? 'categorical' : 'unconvertible';
  }
  if (cell.qualitative) return 'qualitative';
  return 'none';
}

/* ---- biomarker-level metrics ------------------------------------------ */
function predictiveStrength(bm) {
  let best = 0;
  OUTCOMES.forEach(o => {
    const c = bm.outcomes && bm.outcomes[o.key];
    const m = c ? magnitude(comparableHr(c)) : null;
    if (m != null && m > best) best = m;
  });
  return best;
}
function modRank(bm) {
  const r = bm.modifiability && bm.modifiability.rating;
  return MOD_RATINGS[r] ? MOD_RATINGS[r].rank : 0;
}
function priorityScore(bm) {
  return clamp01(predictiveStrength(bm) / LN2) * modRank(bm);
}
function isPriority(bm) {
  const r = bm.modifiability && bm.modifiability.rating;
  return r === 'high' && predictiveStrength(bm) >= 0.22;
}

/* ---- data shaping ------------------------------------------------------ */
function domainList() {
  const defined = (DATA.meta.domains || []).slice().sort((a, b) => a.order - b.order);
  const seen = new Set(defined.map(d => d.id));
  const extra = [];
  DATA.biomarkers.forEach(bm => {
    if (!seen.has(bm.domain)) { seen.add(bm.domain); extra.push({ id: bm.domain, name: bm.domain, order: 999 }); }
  });
  return defined.concat(extra);
}
function biomarkersForDomain(domainId) {
  let list = DATA.biomarkers.filter(bm => bm.domain === domainId);
  if (state.filter) {
    const q = state.filter.toLowerCase();
    list = list.filter(bm => {
      const hay = (bm.name + ' ' + (bm.aliases || []).join(' ')).toLowerCase();
      return hay.includes(q);
    });
  }
  if (state.sort === 'predictive') {
    list = list.slice().sort((a, b) => predictiveStrength(b) - predictiveStrength(a));
  } else if (state.sort === 'priority') {
    list = list.slice().sort((a, b) => priorityScore(b) - priorityScore(a));
  }
  return list;
}

/* ---- rendering --------------------------------------------------------- */
const $map = document.getElementById('map');

function el(tag, cls, html) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html != null) e.innerHTML = html;
  return e;
}

function labelCell(bm) {
  const td = el('td', 'bm-label');
  const name = el('span', 'bm-label__name', escapeHtml(bm.name));
  name.addEventListener('click', () => openModal(bm));
  td.appendChild(name);
  if (isPriority(bm)) td.appendChild(el('span', 'bm-label__tag', '★ PRIORITY'));
  return td;
}

function domainRow(domain, span, count) {
  const tr = el('tr', 'domain-row');
  const td = el('td');
  td.colSpan = span;
  const collapsed = state.collapsed.has(domain.id);
  const caret = collapsed ? '▸' : '▾';
  td.innerHTML = `<span class="domain-toggle">${caret} ${escapeHtml(domain.name)}` +
    `<span class="domain-count">${count} biomarker${count === 1 ? '' : 's'}</span></span>`;
  td.querySelector('.domain-toggle').addEventListener('click', () => {
    if (collapsed) state.collapsed.delete(domain.id);
    else state.collapsed.add(domain.id);
    render();
  });
  tr.appendChild(td);
  return tr;
}

function dataCell(bm, outcomeKey) {
  const cell = bm.outcomes ? bm.outcomes[outcomeKey] : null;
  const td = el('td', 'cell');
  const inner = el('div', 'cell__inner');
  const st = cellState(cell);
  if (st === 'none') {
    td.classList.add('cell--nodata');
    inner.textContent = '—';
    td.appendChild(inner);
    return td;
  }
  if (st === 'qualitative') {
    td.classList.add('cell--qual');
    inner.classList.add('has-data');
    inner.innerHTML = `<span class="cell__qual">${escapeHtml(cell.evidence_tier || '•')}</span>`;
    attachTip(inner, () => tipForCell(bm, outcomeKey, cell));
    td.appendChild(inner);
    return td;
  }
  if (st === 'categorical' || st === 'unconvertible') {
    td.classList.add(st === 'categorical' ? 'cell--cat' : 'cell--unconv');
    inner.classList.add('has-data', tierClass(cell.evidence_tier));
    inner.innerHTML = `<span class="cell__hr">${fmtHr(cell.hr)}</span>` +
      `<span class="cell__sub">${st === 'categorical' ? 'categ.' : 'n/c'}</span>`;
    attachTip(inner, () => tipForCell(bm, outcomeKey, cell));
    td.appendChild(inner);
    return td;
  }
  const bg = effectColor(cell);
  inner.classList.add('has-data', tierClass(cell.evidence_tier));
  inner.style.background = bg;
  inner.style.color = textColorFor(bg);
  inner.innerHTML = `<span class="cell__hr">${fmtHr(cell.hr_per_sd)}</span>`;
  attachTip(inner, () => tipForCell(bm, outcomeKey, cell));
  td.appendChild(inner);
  return td;
}

function modCell(bm) {
  const m = bm.modifiability || {};
  const r = MOD_RATINGS[m.rating];
  const td = el('td', 'cell mod-cell');
  const inner = el('div', 'cell__inner');
  if (!r) {
    td.classList.add('cell--nodata');
    inner.textContent = '—';
  } else {
    inner.classList.add('has-data');
    const dot = el('span', 'mod-dot');
    dot.style.background = r.color;
    inner.appendChild(dot);
    inner.appendChild(el('span', null, r.label));
    attachTip(inner, () => tipForMod(bm, m));
  }
  td.appendChild(inner);
  return td;
}

function render() {
  const table = el('table', 'hmtable');
  const head = el('tr');
  head.appendChild(el('th', 'col-head col-head--label', 'Biomarker'));
  OUTCOMES.forEach(o => {
    const th = el('th', 'col-head' + (o.anchor ? ' col-acm' : ''),
      escapeHtml(o.label).replace(/\n/g, '<br>'));
    head.appendChild(th);
  });
  head.appendChild(el('th', 'col-head', 'Modifiability'));
  const thead = el('thead'); thead.appendChild(head); table.appendChild(thead);

  const tbody = el('tbody');
  const span = OUTCOMES.length + 2;
  domainList().forEach(domain => {
    const bms = biomarkersForDomain(domain.id);
    if (!bms.length) return;
    tbody.appendChild(domainRow(domain, span, bms.length));
    if (state.collapsed.has(domain.id)) return;
    bms.forEach(bm => {
      const tr = el('tr', 'bm-row');
      if (isPriority(bm)) tr.classList.add('is-priority');
      tr.appendChild(labelCell(bm));
      OUTCOMES.forEach(o => {
        const td = dataCell(bm, o.key);
        if (o.anchor) td.classList.add('col-acm');
        tr.appendChild(td);
      });
      tr.appendChild(modCell(bm));
      tbody.appendChild(tr);
    });
  });
  table.appendChild(tbody);
  $map.innerHTML = '';
  $map.appendChild(table);
  renderLegend();
}

/* ---- tooltip ----------------------------------------------------------- */
const $tip = document.getElementById('tooltip');

function attachTip(node, builder) {
  node.addEventListener('mouseenter', () => {
    $tip.innerHTML = builder();
    $tip.hidden = false;
  });
  node.addEventListener('mousemove', moveTip);
  node.addEventListener('mouseleave', () => { $tip.hidden = true; });
}
function moveTip(e) {
  if ($tip.hidden) return;
  const pad = 14;
  let x = e.clientX + pad, y = e.clientY + pad;
  const r = $tip.getBoundingClientRect();
  if (x + r.width > window.innerWidth - 8) x = e.clientX - r.width - pad;
  if (y + r.height > window.innerHeight - 8) y = e.clientY - r.height - pad;
  $tip.style.left = Math.max(8, x) + 'px';
  $tip.style.top = Math.max(8, y) + 'px';
}

function tipForCell(bm, outcomeKey, cell) {
  const oLabel = OUTCOMES.find(o => o.key === outcomeKey).label.replace(/\n/g, ' ');
  const ci = (cell.ci_low != null && cell.ci_high != null)
    ? `${fmtHr(cell.ci_low)}–${fmtHr(cell.ci_high)}` : 'not reported';
  const flag = cell.verification_status || 'unverified';
  let html = `<h5>${escapeHtml(bm.name)} &mdash; ${escapeHtml(oLabel)}</h5><dl>`;
  if (typeof cell.hr_per_sd === 'number') {
    const ciSd = cell.hr_per_sd_ci
      ? `${fmtHr(cell.hr_per_sd_ci[0])}–${fmtHr(cell.hr_per_sd_ci[1])}` : 'n/a';
    html += row('HR per SD', `<b>${fmtHr(cell.hr_per_sd)}</b> (${ciSd}) ` +
      `<span style="color:#a9b4bd">${escapeHtml(cell.direction || '')}</span>`);
    html += row('Native', `${fmtHr(cell.hr)} (${ci}) · ${escapeHtml(cell.hr_unit_detail || cell.hr_metric || '')}`);
    html += row('Standardized', escapeHtml((cell.comparability || '') +
      (cell.hr_per_sd_basis ? ' — ' + cell.hr_per_sd_basis : '')));
  } else if (typeof cell.hr === 'number') {
    html += row('Native HR', `${fmtHr(cell.hr)} (${ci}) ` +
      `<span style="color:#a9b4bd">(${escapeHtml(cell.direction || '')})</span>`);
    html += row('Metric', escapeHtml(cell.hr_unit_detail || cell.hr_metric || '—'));
    html += row('Standardized', `<em>not per-SD comparable — ${escapeHtml(cell.comparability || 'n/a')}</em>`);
  } else {
    html += row('HR', '<em>evidence documented; no usable pooled HR</em>');
    html += row('Metric', escapeHtml(cell.hr_unit_detail || cell.hr_metric || '—'));
  }
  html += row('Evidence', `Tier ${escapeHtml(cell.evidence_tier || '?')}`);
  html += row('Study type', escapeHtml(cell.study_type || '—'));
  html += row('Sample', cell.n != null ? escapeHtml(String(cell.n)) : '—');
  html += row('Follow-up', cell.followup_years != null ? cell.followup_years + ' y' : '—');
  html += '</dl>';
  if (cell.study) html += `<div class="tt-cite">${escapeHtml(cell.study)}</div>`;
  if (cell.citation) html += `<div class="tt-cite">${escapeHtml(cell.citation)}</div>`;
  if (cell.notes) html += `<div class="tt-cite"><em>${escapeHtml(cell.notes)}</em></div>`;
  html += `<span class="tt-flag flag-${flag}">${flag}</span>`;
  return html;
}
function tipForMod(bm, m) {
  let html = `<h5>${escapeHtml(bm.name)} &mdash; Modifiability</h5><dl>`;
  html += row('Rating', escapeHtml((MOD_RATINGS[m.rating] || {}).label || m.rating || '—'));
  html += row('Intervention', escapeHtml(m.intervention || '—'));
  html += row('Effect size', escapeHtml(m.effect_size || '—'));
  html += row('Timeframe', escapeHtml(m.timeframe || '—'));
  html += row('Evidence', escapeHtml(m.evidence_type || '—'));
  html += '</dl>';
  if (m.citation) html += `<div class="tt-cite">${escapeHtml(m.citation)}</div>`;
  return html;
}
function row(k, v) { return `<dt>${k}</dt><dd>${v}</dd>`; }

/* ---- modal: biomarker profile ----------------------------------------- */
const $modal = document.getElementById('modal');

function openModal(bm) {
  const p = bm.profile || {};
  const meas = p.measurement || {};
  const mod = bm.modifiability || {};
  const domain = (DATA.meta.domains || []).find(d => d.id === bm.domain);
  let html = `<div class="modal__card">
    <button class="modal__close" aria-label="Close">&times;</button>
    <h2>${escapeHtml(bm.name)}</h2>
    <p class="modal__domain">${escapeHtml(domain ? domain.name : bm.domain)}</p>`;

  html += section('Mechanism', `<p>${escapeHtml(p.mechanism || 'Not yet documented.')}</p>`);

  html += `<section><h3>Measurement</h3><dl class="modal__grid">
    ${dRow('Method', meas.method)}
    ${dRow('Sample type', meas.sample_type)}
    ${dRow('Approx. cost', meas.approx_cost_usd ? '$' + escapeHtml(meas.approx_cost_usd) : null)}
    ${dRow('Signal quality', p.signal_quality
      ? escapeHtml(p.signal_quality) + (p.signal_quality_note ? ' — ' + escapeHtml(p.signal_quality_note) : '')
      : null)}
    ${dRow('Testing cadence', p.testing_cadence)}
  </dl></section>`;

  if (p.confounders && p.confounders.length) {
    html += `<section><h3>Major confounders</h3><div class="modal__chips">` +
      p.confounders.map(c => `<span class="modal__chip">${escapeHtml(c)}</span>`).join('') +
      `</div></section>`;
  }

  html += `<section><h3>Modifiability</h3><dl class="modal__grid">
    ${dRow('Rating', (MOD_RATINGS[mod.rating] || {}).label)}
    ${dRow('Intervention', mod.intervention)}
    ${dRow('Effect size', mod.effect_size)}
    ${dRow('Timeframe', mod.timeframe)}
    ${dRow('Evidence type', mod.evidence_type)}
    ${dRow('Citation', mod.citation)}
  </dl></section>`;

  const oRows = OUTCOMES.map(o => {
    const c = bm.outcomes && bm.outcomes[o.key];
    const lbl = o.label.replace(/\n/g, ' ');
    const st = cellState(c);
    if (st === 'none') return dRow(lbl, 'no usable evidence');
    if (st === 'qualitative') return dRow(lbl,
      'documented; effect size not quantified — ' + (c.study || c.notes || ''));
    if (st === 'comparable') {
      const ci = c.hr_per_sd_ci ? ` (${fmtHr(c.hr_per_sd_ci[0])}–${fmtHr(c.hr_per_sd_ci[1])})` : '';
      return dRow(lbl, `HR/SD ${fmtHr(c.hr_per_sd)}${ci}, tier ${c.evidence_tier} — ${c.study || ''}`);
    }
    const ci = (c.ci_low != null) ? ` (${fmtHr(c.ci_low)}–${fmtHr(c.ci_high)})` : '';
    return dRow(lbl, `native HR ${fmtHr(c.hr)}${ci} [${c.comparability}], tier ${c.evidence_tier} — ${c.study || ''}`);
  }).join('');
  html += `<section><h3>Outcome associations</h3><dl class="modal__grid">${oRows}</dl></section>`;

  html += `</div>`;
  $modal.innerHTML = html;
  $modal.hidden = false;
  $modal.querySelector('.modal__close').addEventListener('click', closeModal);
  $modal.addEventListener('click', e => { if (e.target === $modal) closeModal(); });
}
function closeModal() { $modal.hidden = true; $modal.innerHTML = ''; }
function section(title, body) { return `<section><h3>${title}</h3>${body}</section>`; }
function dRow(k, v) {
  if (v == null || v === '') return '';
  return `<dt>${escapeHtml(k)}</dt><dd>${escapeHtml(String(v))}</dd>`;
}

/* ---- legend ------------------------------------------------------------ */
function renderLegend() {
  const $l = document.getElementById('legend');
  const gradient = rampStops => {
    const stops = [];
    for (let i = 0; i <= 6; i++) stops.push(ramp(rampStops, i / 6) + ' ' + Math.round(i / 6 * 100) + '%');
    return stops.join(',');
  };
  let html = `<div class="legend__block">
    <h4>Effect magnitude (HR per +1 SD)</h4>
    <div class="legend__ramp">
      <span>weak</span>
      <span class="legend__swatches" style="width:120px;height:14px;
        background:linear-gradient(90deg,${gradient(PROT_RAMP)})"></span>
      <span>protective</span>
    </div>
    <div class="legend__ramp" style="margin-top:5px">
      <span>weak</span>
      <span class="legend__swatches" style="width:120px;height:14px;
        background:linear-gradient(90deg,${gradient(RISK_RAMP)})"></span>
      <span>risk</span>
    </div>
  </div>`;

  html += `<div class="legend__block"><h4>Evidence tier</h4>
    <div class="legend__row"><span class="legend__chip tier-A" style="background:#fff"></span>A &mdash; multiple MA / MA + MR</div>
    <div class="legend__row"><span class="legend__chip tier-B" style="background:#fff"></span>B &mdash; single MA / strong cohort</div>
    <div class="legend__row"><span class="legend__chip tier-C" style="background:#fff"></span>C &mdash; limited / observational</div>
  </div>`;

  html += `<div class="legend__block"><h4>Off the per-SD scale</h4>
    <div class="legend__row"><span class="legend__chip legend__chip--cat"></span>categorical exposure (native HR shown)</div>
    <div class="legend__row"><span class="legend__chip legend__chip--unconv"></span>not standardizable (flagged)</div>
    <div class="legend__row"><span class="legend__chip legend__chip--qual"></span>documented, no usable pooled HR</div>
    <div class="legend__row"><span class="legend__chip cell--nodata" style="background:#f0f0ee"></span>no usable evidence</div>
  </div>`;

  html += `<div class="legend__block"><h4>Modifiability</h4>` +
    Object.keys(MOD_RATINGS).map(k =>
      `<div class="legend__row"><span class="mod-dot" style="background:${MOD_RATINGS[k].color}"></span>${MOD_RATINGS[k].label}</div>`
    ).join('') +
    `<div class="legend__row" style="margin-top:4px;color:var(--priority)">★ priority = high predictive &amp; high modifiable</div>
  </div>`;

  $l.innerHTML = html;
}

/* ---- formatting -------------------------------------------------------- */
function fmtHr(v) { return (Math.round(v * 100) / 100).toFixed(2); }
function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

/* ---- boot -------------------------------------------------------------- */
function wireControls() {
  document.getElementById('sort-select').addEventListener('change', e => {
    state.sort = e.target.value; render();
  });
  let t;
  document.getElementById('filter-input').addEventListener('input', e => {
    clearTimeout(t);
    t = setTimeout(() => { state.filter = e.target.value.trim(); render(); }, 140);
  });
  window.addEventListener('keydown', e => { if (e.key === 'Escape') closeModal(); });
}

function showMeta() {
  const n = DATA.biomarkers.length;
  const cells = DATA.biomarkers.reduce((acc, bm) => {
    OUTCOMES.forEach(o => { if (bm.outcomes && bm.outcomes[o.key]) acc++; });
    return acc;
  }, 0);
  const domains = new Set(DATA.biomarkers.map(b => b.domain)).size;
  document.getElementById('meta-line').textContent =
    `${n} biomarkers · ${domains} domains · ${cells} evidence cells`;
  document.getElementById('footnote-version').textContent =
    `Data layer v${DATA.meta.version || '?'} · updated ${DATA.meta.last_updated || 'n/a'}. ` +
    `Source: data/biomarkers.json`;
}

function showLoadError() {
  $map.innerHTML = `<div class="loaderr">
    <h3>Could not load the data layer</h3>
    <p>The heat map expects <code style="display:inline">heatmap/data.bundle.js</code>
       (a generated mirror of <code style="display:inline">data/biomarkers.json</code>).
       If it is missing, regenerate it:</p>
    <code>node tools/standardize.js</code>
    <p>Run that from the project root, then reload this page.</p>
  </div>`;
}

function init(json) {
  DATA = json;
  DATA.meta = DATA.meta || {};
  DATA.biomarkers = DATA.biomarkers || [];
  showMeta();
  render();
}

function boot() {
  wireControls();
  if (window.BIOMARKER_DATA) { init(window.BIOMARKER_DATA); return; }
  fetch('../data/biomarkers.json', { cache: 'no-store' })
    .then(r => { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then(init)
    .catch(showLoadError);
}
document.addEventListener('DOMContentLoaded', boot);
