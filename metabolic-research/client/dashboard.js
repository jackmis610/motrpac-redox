/* Personal biomarker dashboard — scores your measurements against the
   verified evidence map. Vanilla JS; loads refdata.js + measurements.js. */

'use strict';

const OUTCOMES = [
  { key: 'all_cause_mortality', label: 'all-cause mortality' },
  { key: 'cvd', label: 'cardiovascular disease' },
  { key: 'cancer', label: 'cancer' },
  { key: 'dementia', label: 'dementia' },
  { key: 'frailty', label: 'frailty' }
];
const MOD_RANK = { high: 1.0, moderate: 0.66, low: 0.33, fixed: 0 };

let BM = {}, DOMAINS = {};

function ln(x) { return Math.log(x); }
function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }
function fmt(x, n) { return (typeof x === 'number') ? x.toFixed(n == null ? 2 : n) : '—'; }
function esc(s) {
  return String(s == null ? '' : s).replace(/[&<>"]/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

/* ---- per-biomarker evidence helpers ----------------------------------- */
function comparableHr(cell) {
  return (cell && typeof cell.hr_per_sd === 'number') ? cell.hr_per_sd : null;
}
function predictiveStrength(bm) {
  let best = 0;
  OUTCOMES.forEach(o => {
    const h = comparableHr(bm.outcomes && bm.outcomes[o.key]);
    if (h != null) best = Math.max(best, Math.abs(ln(h)));
  });
  return best;
}
function modRank(bm) {
  const r = bm.modifiability && bm.modifiability.rating;
  return MOD_RANK[r] != null ? MOD_RANK[r] : 0;
}
/* the outcome used for modeled hazard: all-cause mortality if comparable,
   else the comparable outcome with the largest effect */
function primaryOutcome(bm) {
  const acm = bm.outcomes && bm.outcomes.all_cause_mortality;
  if (comparableHr(acm) != null) return OUTCOMES[0];
  let best = null, bestMag = 0;
  OUTCOMES.forEach(o => {
    const h = comparableHr(bm.outcomes && bm.outcomes[o.key]);
    if (h != null && Math.abs(ln(h)) > bestMag) { bestMag = Math.abs(ln(h)); best = o; }
  });
  return best;
}

/* ---- measurement helpers ---------------------------------------------- */
function history(id) {
  return (window.MEASUREMENTS || [])
    .filter(p => p.values && p.values[id] != null)
    .map(p => ({ date: p.date, value: p.values[id], source: p.source }))
    .sort((a, b) => a.date < b.date ? -1 : 1);
}
function trackedIds() {
  const s = new Set();
  (window.MEASUREMENTS || []).forEach(p => Object.keys(p.values || {}).forEach(id => s.add(id)));
  return [...s];
}

/* ---- scoring ----------------------------------------------------------- */
function toScale(ref, v) { return ref.scale === 'log' ? ln(v) : v; }
function zOf(ref, v) { return (toScale(ref, v) - ref.population_mean) / ref.population_sd; }

/* signed SD distance of a value above the optimal band (marker-value units) */
function dzAboveOptimal(ref, v) {
  const z = zOf(ref, v);
  const o = ref.optimal || {};
  let zlo = (o.low != null) ? zOf(ref, o.low) : -Infinity;
  let zhi = (o.high != null) ? zOf(ref, o.high) : Infinity;
  if (zlo > zhi) { const t = zlo; zlo = zhi; zhi = t; }
  if (z > zhi && isFinite(zhi)) return z - zhi;
  if (z < zlo && isFinite(zlo)) return z - zlo;
  return 0;
}

function score(id) {
  const bm = BM[id];
  const ref = (window.REFERENCE_RANGES.ranges || {})[id];
  const hist = history(id);
  if (!bm) return { id, state: 'unknown', hist };
  if (!ref || !ref.applicable || ref.population_mean == null || ref.population_sd == null)
    return { id, bm, ref, state: 'tracked', hist };

  const latest = hist[hist.length - 1];
  const prev = hist.length > 1 ? hist[hist.length - 2] : null;
  const po = primaryOutcome(bm);
  const hr = po ? comparableHr(bm.outcomes[po.key]) : null;

  function modeled(v) {
    if (hr == null) return null;
    return Math.pow(hr, dzAboveOptimal(ref, v));
  }
  const userZ = zOf(ref, latest.value);
  const dz = dzAboveOptimal(ref, latest.value);
  const modHR = modeled(latest.value);
  const prevModHR = prev ? modeled(prev.value) : null;

  const pred = predictiveStrength(bm);
  const priority = clamp(Math.abs(dz), 0, 2) / 2 * pred * modRank(bm);
  const isPriority = Math.abs(dz) >= 0.5 && pred >= 0.15 && modRank(bm) >= 0.66;

  return {
    id, bm, ref, hist, latest, prev, state: 'scored',
    primaryOutcome: po, hr, userZ, dz, modHR, prevModHR,
    prevZ: prev ? zOf(ref, prev.value) : null,
    priority, isPriority,
    statusClass: Math.abs(dz) < 0.001 ? 'flat' : Math.abs(dz) < 1 ? 'warn' : 'bad'
  };
}

/* ---- rendering --------------------------------------------------------- */
const $view = document.getElementById('view');
let STATE = { view: 'priorities' };

function el(tag, cls, html) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html != null) e.innerHTML = html;
  return e;
}

function sdBar(s) {
  const pct = z => clamp((z + 3) / 6 * 100, 0, 100);
  const ref = s.ref, o = ref.optimal || {};
  let zlo = o.low != null ? zOf(ref, o.low) : -3;
  let zhi = o.high != null ? zOf(ref, o.high) : 3;
  if (zlo > zhi) { const t = zlo; zlo = zhi; zhi = t; }
  const bandL = pct(zlo), bandR = pct(zhi);
  let h = '<div class="bar"><div class="bar__track"></div>';
  h += `<div class="bar__optimal" style="left:${bandL}%;width:${Math.max(bandR - bandL, 1.5)}%"></div>`;
  if (s.prev) h += `<div class="bar__prev" style="left:calc(${pct(s.prevZ)}% - 4px)"></div>`;
  h += `<div class="bar__you" style="left:${pct(s.userZ)}%"></div>`;
  h += '</div><div class="bar__scale"><span>&minus;3 SD</span><span>population mean</span><span>+3 SD</span></div>';
  return h;
}

function statusText(s) {
  if (Math.abs(s.dz) < 0.001) return '<span class="status status--good">within your optimal range</span>';
  const dir = s.dz > 0 ? 'above' : 'below';
  return `<span class="status status--${s.statusClass}">${fmt(Math.abs(s.dz), 1)} SD ${dir} optimal</span>`;
}

function modeledText(s) {
  if (s.modHR == null) return '<span class="modeled">no per-SD outcome model for this marker</span>';
  let h = `<span class="modeled">modeled ${esc(s.primaryOutcome.label)} hazard ` +
    `<b>${fmt(s.modHR)}&times;</b> vs your optimal`;
  if (s.prevModHR != null) {
    const better = s.modHR < s.prevModHR - 0.005;
    const worse = s.modHR > s.prevModHR + 0.005;
    const cls = better ? 'better' : worse ? 'worse' : 'flat';
    const arrow = better ? '▼' : worse ? '▲' : '→';
    h += ` &nbsp;<span class="delta delta--${cls}">${arrow} from ${fmt(s.prevModHR)}&times; ` +
      `(${esc(s.prev.date)})</span>`;
  }
  return h + '</span>';
}

function card(s) {
  const c = el('div', 'card' + (s.isPriority ? ' is-priority' : ''));
  const mod = s.bm.modifiability || {};
  c.innerHTML =
    `<div class="card__name">${esc(s.bm.name)}` +
      (s.isPriority ? '<span class="tag">★ PRIORITY</span>' : '') + '</div>' +
    `<div class="card__right">${sdBar(s)}</div>` +
    `<div class="card__sub">` +
      `<span class="card__val"><b>${esc(fmtVal(s))}</b> ${esc(s.ref.units || '')}` +
      ` &middot; ${esc(s.latest.date)}</span> &nbsp; ${statusText(s)}<br>` +
      `${modeledText(s)}<br>` +
      `<span class="card__sub">Lever: ${esc(mod.intervention || '—')}</span>` +
    `</div>`;
  c.addEventListener('click', () => openModal(s));
  return c;
}
function fmtVal(s) {
  const v = s.latest.value;
  return (Math.abs(v) >= 100 || Number.isInteger(v)) ? String(v) : v.toFixed(1);
}

function renderPriorities(scored) {
  const note = el('p', 'sectionnote',
    'Your biomarkers ranked by <strong>priority</strong> — predictive weight &times; ' +
    'modifiability &times; how far you sit from your optimal range. Markers you are ' +
    'already optimal on fall to the bottom; ★ flags high-leverage targets.');
  $view.appendChild(note);
  const sorted = scored.slice().sort((a, b) => b.priority - a.priority);
  if (!sorted.length) { $view.appendChild(el('p', 'sectionnote', 'No scored biomarkers yet.')); return; }
  sorted.forEach(s => $view.appendChild(card(s)));
}

function renderPanel(scored, tracked) {
  $view.appendChild(el('p', 'sectionnote',
    'Every measured biomarker with a reference range, grouped by domain.'));
  const order = (window.BIOMARKER_DATA.meta.domains || []).slice().sort((a, b) => a.order - b.order);
  order.forEach(dom => {
    const inDom = scored.filter(s => s.bm.domain === dom.id);
    if (!inDom.length) return;
    $view.appendChild(el('div', 'domain-head', esc(dom.name)));
    inDom.sort((a, b) => b.priority - a.priority).forEach(s => $view.appendChild(card(s)));
  });
}

function renderTracked(tracked) {
  $view.appendChild(el('p', 'sectionnote',
    'Biomarkers you submitted that are recorded but not scored — either not yet ' +
    'in the evidence model, or not a single quantitative value (genotypes, ' +
    'screening tests, age-acceleration scores).'));
  if (!tracked.length) { $view.appendChild(el('p', 'sectionnote', 'None.')); return; }
  const list = el('div', 'tracked-list');
  tracked.forEach(s => {
    const name = s.bm ? s.bm.name : s.id;
    const why = s.state === 'unknown' ? 'not in the 135-biomarker model'
      : 'no per-SD reference range';
    const last = s.hist[s.hist.length - 1];
    list.appendChild(el('div', null,
      `<strong>${esc(name)}</strong> — ${esc(String(last.value))} ` +
      `<span style="color:var(--ink-faint)">(${esc(why)})</span>`));
  });
  $view.appendChild(list);
}

function render() {
  document.querySelectorAll('.tab').forEach(t =>
    t.classList.toggle('is-active', t.dataset.view === STATE.view));
  $view.innerHTML = '';

  const ids = trackedIds();
  if (!ids.length) {
    $view.appendChild(el('div', 'empty',
      '<strong>No measurements yet.</strong><br>Add your biomarker panels to ' +
      '<code>client/measurements.js</code> and reload. The dashboard works from ' +
      'a single panel; a second panel adds the trajectory view.'));
    return;
  }
  const scored = [], tracked = [];
  ids.forEach(id => {
    const s = score(id);
    if (s.state === 'scored') scored.push(s);
    else tracked.push(s);
  });

  if ((window.MEASUREMENTS || []).every(p => /example/i.test(p.source || ''))) {
    $view.appendChild(el('div', 'empty',
      'Showing <strong>example data</strong> — replace the panels in ' +
      '<code>client/measurements.js</code> with your own.'));
  }

  if (STATE.view === 'priorities') renderPriorities(scored);
  else if (STATE.view === 'panel') renderPanel(scored, tracked);
  else renderTracked(tracked);
}

/* ---- modal ------------------------------------------------------------- */
const $modal = document.getElementById('modal');
function openModal(s) {
  const bm = s.bm, p = bm.profile || {}, mod = bm.modifiability || {};
  let h = `<div class="modal__card"><button class="modal__close">&times;</button>` +
    `<h2>${esc(bm.name)}</h2>`;
  h += `<h3>Your readings</h3><table class="hist"><tr><th>Date</th><th>Value</th>` +
    `<th>SD position</th><th>Modeled hazard vs optimal</th></tr>`;
  s.hist.forEach(r => {
    const z = zOf(s.ref, r.value);
    const mh = (s.hr != null) ? Math.pow(s.hr, dzAboveOptimal(s.ref, r.value)) : null;
    h += `<tr><td>${esc(r.date)}</td><td>${esc(String(r.value))} ${esc(s.ref.units || '')}</td>` +
      `<td>${z >= 0 ? '+' : ''}${fmt(z, 2)} SD</td>` +
      `<td>${mh != null ? fmt(mh) + '×' : '—'}</td></tr>`;
  });
  h += `</table>`;
  h += `<h3>Why it matters</h3><p>${esc(p.mechanism || '')}</p>`;
  h += `<h3>Evidence — verified per-SD hazard ratios</h3><table><tr><th>Outcome</th>` +
    `<th>HR / SD</th><th>Tier</th></tr>`;
  OUTCOMES.forEach(o => {
    const c = bm.outcomes && bm.outcomes[o.key];
    if (!c) return;
    const hps = comparableHr(c);
    h += `<tr><td>${esc(o.label)}</td>` +
      `<td>${hps != null ? fmt(hps) : (typeof c.hr === 'number' ? fmt(c.hr) + ' (native)' : '—')}</td>` +
      `<td>${esc(c.evidence_tier || '—')}</td></tr>`;
  });
  h += `</table>`;
  h += `<h3>Optimal range</h3><p>${esc(s.ref.units || '')}: ` +
    `${esc(String((s.ref.optimal || {}).low))}–${esc(String((s.ref.optimal || {}).high))}. ` +
    `${esc((s.ref.optimal || {}).rationale || '')}</p>`;
  h += `<h3>Lever</h3><p>${esc(mod.intervention || '—')} — <em>${esc(mod.effect_size || '')}</em>, ` +
    `${esc(mod.timeframe || '')}.</p>`;
  h += `<p style="color:var(--ink-faint);font-size:11.5px;margin-top:12px">` +
    `Population reference: ${esc(s.ref.population_source || 'n/a')} · ` +
    `verification: ${esc(s.ref.verification_status || 'n/a')}</p>`;
  h += `</div>`;
  $modal.innerHTML = h;
  $modal.hidden = false;
  $modal.querySelector('.modal__close').onclick = () => { $modal.hidden = true; };
  $modal.onclick = e => { if (e.target === $modal) $modal.hidden = true; };
}

/* ---- boot -------------------------------------------------------------- */
function boot() {
  if (!window.BIOMARKER_DATA || !window.REFERENCE_RANGES) {
    $view.innerHTML = '<div class="empty"><strong>Data not loaded.</strong><br>' +
      'Run <code>node tools/build-client.js</code> to generate client/refdata.js.</div>';
    return;
  }
  (window.BIOMARKER_DATA.biomarkers || []).forEach(b => { BM[b.id] = b; });
  const n = trackedIds().length;
  const panels = (window.MEASUREMENTS || []).length;
  document.getElementById('meta-line').textContent =
    `${n} biomarkers tracked · ${panels} panel${panels === 1 ? '' : 's'} · ` +
    `evidence map v${window.BIOMARKER_DATA.meta.version || '?'}`;
  document.querySelectorAll('.tab').forEach(t =>
    t.addEventListener('click', () => { STATE.view = t.dataset.view; render(); }));
  window.addEventListener('keydown', e => { if (e.key === 'Escape') $modal.hidden = true; });
  render();
}
document.addEventListener('DOMContentLoaded', boot);
