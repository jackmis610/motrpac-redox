#!/usr/bin/env node
/* Internal QA audit of data/biomarkers.json. Checks structural integrity, the
   standardization math, CI/HR consistency, evidence-tier coherence, and
   plausibility. Does NOT verify citations against primary sources — that is a
   separate web-research step.

   Usage:  node tools/audit.js [path/to/biomarkers.json] */

'use strict';
const fs = require('fs');

const FILE = process.argv[2] || 'data/biomarkers.json';
const OUTCOMES = ['all_cause_mortality', 'cvd', 'cancer', 'dementia', 'frailty'];
const SPREAD = { tertile: 2.18, quartile: 2.54, quintile: 2.80, median_split: 1.60, decile: 3.51 };
const METRICS = ['per_sd', 'per_log_sd', 'tertile', 'quartile', 'quintile',
  'decile', 'median_split', 'per_unit', 'per_doubling', 'categorical'];

const data = JSON.parse(fs.readFileSync(FILE, 'utf8'));
const errors = [], warns = [];
const err = m => errors.push(m);
const warn = m => warns.push(m);

const domainIds = new Set((data.meta.domains || []).map(d => d.id));
const ids = new Set(), names = new Set();
const vstat = { verified: 0, approximate: 0, unverified: 0, missing: 0 };
let quant = 0, qual = 0, nul = 0;

for (const b of data.biomarkers) {
  const at = b.id || b.name || '?';
  if (!b.id) err(`${at}: missing id`);
  if (ids.has(b.id)) err(`duplicate id: ${b.id}`);
  ids.add(b.id);
  if (names.has(b.name)) warn(`duplicate name: ${b.name}`);
  names.add(b.name);
  if (!domainIds.has(b.domain)) err(`${at}: domain "${b.domain}" not in meta.domains`);
  for (const k of ['name', 'profile', 'outcomes', 'modifiability'])
    if (!b[k]) err(`${at}: missing ${k}`);

  // comparability tag consistency
  let offScale = false;

  for (const o of OUTCOMES) {
    if (!b.outcomes || !(o in b.outcomes)) { err(`${at}: outcomes missing key ${o}`); continue; }
    const c = b.outcomes[o];
    if (c == null) { nul++; continue; }
    if (typeof c.hr !== 'number') { qual++; continue; }
    quant++;
    const cat = `${b.id}/${o}`;

    // required fields
    for (const k of ['direction', 'hr_metric', 'evidence_tier', 'study', 'citation', 'verification_status'])
      if (c[k] == null || c[k] === '') warn(`${cat}: missing ${k}`);

    if (!METRICS.includes(c.hr_metric)) err(`${cat}: bad hr_metric "${c.hr_metric}"`);
    if (!['risk', 'protective'].includes(c.direction)) err(`${cat}: bad direction "${c.direction}"`);
    if (!['A', 'B', 'C'].includes(c.evidence_tier)) err(`${cat}: bad evidence_tier "${c.evidence_tier}"`);
    if (vstat[c.verification_status] != null) vstat[c.verification_status]++;
    else vstat.missing++;

    // CI sanity
    if (typeof c.ci_low === 'number' && typeof c.ci_high === 'number') {
      if (c.ci_low > c.ci_high) err(`${cat}: ci_low ${c.ci_low} > ci_high ${c.ci_high}`);
      if (c.hr < c.ci_low - 0.02 || c.hr > c.ci_high + 0.02)
        err(`${cat}: hr ${c.hr} outside CI [${c.ci_low}, ${c.ci_high}]`);
      if (c.ci_low <= 0) err(`${cat}: ci_low <= 0`);
    }
    if (c.hr <= 0) err(`${cat}: hr <= 0`);

    // evidence-tier coherence
    if (c.evidence_tier === 'A' &&
        !['meta-analysis', 'mendelian-randomization', 'pooled-cohort'].includes(c.study_type))
      warn(`${cat}: tier A but study_type "${c.study_type}"`);

    // recompute hr_per_sd independently and compare
    let k = null;
    const m = c.hr_metric, cv = c.conversion || {};
    if (m === 'per_sd' || m === 'per_log_sd') k = 1;
    else if (SPREAD[m]) k = 1 / SPREAD[m];
    else if (m === 'per_unit' && cv.unit_value > 0 && typeof cv.population_sd === 'number')
      k = cv.population_sd / cv.unit_value;
    else if (m === 'per_doubling' && typeof cv.log2_sd === 'number') k = cv.log2_sd;

    if (m === 'per_unit' && k == null) {
      if (c.comparability !== 'unconvertible') warn(`${cat}: per_unit lacks conversion inputs but comparability="${c.comparability}"`);
    }
    if (k != null) {
      const expect = Math.round(Math.pow(c.hr, k) * 1e4) / 1e4;
      if (typeof c.hr_per_sd !== 'number')
        err(`${cat}: convertible (k=${k.toFixed(3)}) but hr_per_sd is ${c.hr_per_sd}`);
      else if (Math.abs(c.hr_per_sd - expect) > 0.01)
        err(`${cat}: hr_per_sd ${c.hr_per_sd} != recomputed ${expect}`);
    }
    if (m === 'categorical' && c.hr_per_sd != null)
      err(`${cat}: categorical but hr_per_sd is set`);

    // plausibility
    if (typeof c.hr_per_sd === 'number') {
      const lm = Math.abs(Math.log(c.hr_per_sd));
      if (lm > Math.log(2.5)) warn(`${cat}: extreme per-SD HR ${c.hr_per_sd} — scrutinize`);
    }
    if (typeof c.followup_years === 'number' && (c.followup_years <= 0 || c.followup_years > 60))
      warn(`${cat}: implausible followup_years ${c.followup_years}`);
    if (c.citation && !/\b(19|20)\d\d\b/.test(String(c.citation)))
      warn(`${cat}: citation has no year — "${String(c.citation).slice(0, 50)}"`);

    if (c.comparability === 'categorical' || c.comparability === 'unconvertible') offScale = true;
  }

  if (offScale && !b.comparability_tag) warn(`${at}: has off-scale cells but no comparability_tag`);
  if (!offScale && b.comparability_tag) warn(`${at}: comparability_tag set but no off-scale cells`);
}

console.log(`AUDIT ${FILE}`);
console.log(`  ${data.biomarkers.length} biomarkers · ${quant} quantified · ${qual} qualitative · ${nul} null`);
console.log(`  verification_status: ` + JSON.stringify(vstat));
console.log(`\nERRORS (${errors.length}):`);
errors.forEach(e => console.log('  ✗ ' + e));
console.log(`\nWARNINGS (${warns.length}):`);
warns.forEach(w => console.log('  ! ' + w));
process.exitCode = errors.length ? 1 : 0;
