# Longevity Biomarker Heat Map

A rigorous, citation-backed heat map of longevity biomarkers across multiple
framings, with an intervention-modifiability layer. Built as both a defensible
academic artifact and a clinical coaching tool for a metabolic longevity
assessment protocol.

## Status

**Complete and verified — all 14 domains.** 135 biomarkers, 519 evidence cells
(356 quantified, 163 documented-but-unquantified); 241 cells standardized to a
comparable HR-per-SD scale. Every quantified cell has been checked against
primary sources across two passes — **257 verified, 99 approximate, 0
unverified**; cells that could not be anchored to a real source were nulled, not
retained. The data structure, the interactive heat map, the per-biomarker
profiles, and the synthesis writeup are all in place.

| # | Domain | Biomarkers |
|---|---|---|
| 1 | Insulin Sensitivity / Glucose Metabolism | 16 |
| 2 | Lipids / Cardiovascular Risk | 15 |
| 3 | Mitochondrial / Cardiorespiratory Fitness | 14 |
| 4 | Inflammation | 9 |
| 5 | Biological Age Clocks | 8 |
| 6 | Hormonal Axes | 21 |
| 7 | Body Composition & Anthropometrics | 7 |
| 8 | Cardiovascular & Autonomic Function | 6 |
| 9 | Renal & Hepatic Function | 8 |
| 10 | Cognitive & Neurological | 6 |
| 11 | Sleep & Recovery | 6 |
| 12 | Micronutrient & Metabolic Cofactors | 13 |
| 13 | Cancer Screening & Risk | 3 |
| 14 | Genetic & Pharmacogenomic | 3 |

## Layout

```
metabolic-research/
├── data/
│   ├── biomarkers.json        the data layer — biomarkers, evidence cells, modifiability
│   ├── SCHEMA.md              field-by-field schema documentation
│   └── HR_STANDARDIZATION.md  how hazard ratios are converted to a common per-SD scale
├── tools/
│   ├── standardize.js         deterministic HR → per-SD converter (idempotent)
│   └── build-standalone.js    bundles the heat map into one portable HTML file
├── heatmap/
│   ├── index.html             interactive 3-view heat map
│   ├── styles.css
│   ├── app.js                 vanilla JS, no framework, no build step
│   ├── data.bundle.js         generated — data layer mirrored for file:// use
│   └── heatmap-standalone.html generated — single self-contained file
├── profiles/
│   └── biomarker-profiles.md  Deliverable 1 — per-biomarker profiles
└── writeup/
    └── summary.md             Deliverable 5 — synthesis writeup
```

The data layer (`data/biomarkers.json`) is deliberately separate from the
visualization so biomarkers can be added or re-cited without touching the
heat map code. See `data/SCHEMA.md` for the contract.

## Running the heat map

Two options, both server-free — just open the file in any browser:

- `heatmap/index.html` — the multi-file version (loads `data.bundle.js`).
- `heatmap/heatmap-standalone.html` — a single self-contained file (stylesheet,
  data, and script all inlined). The most portable; good for sharing.

`data/biomarkers.json` remains the editable source of truth. After changing it,
regenerate everything downstream:

```sh
node tools/standardize.js        # per-SD values + heatmap/data.bundle.js
node tools/build-standalone.js   # heatmap/heatmap-standalone.html
```

## The three views

- **View A — All-Cause Mortality.** Biomarkers placed in effect-size columns;
  cell fill encodes evidence tier. The cleanest defensible academic view.
- **View B — Domain-Stratified Outcomes.** HR magnitude across all-cause
  mortality, CVD, cancer, dementia, and frailty; cell border = evidence tier.
- **View C — Outcomes + Modifiability.** View B plus an intervention-
  modifiability column; high-predictive / high-modifiable markers flagged as
  priority targets. The coaching view.

## Evidence standards

- Hazard ratios prioritize meta-analyses, Mendelian randomization, and large
  prospective cohorts (UK Biobank, ARIC, Framingham, MESA, Whitehall II, etc.).
- HR units differ across source studies; every continuous marker is standardized
  to HR per +1 SD by `tools/standardize.js` (method in `HR_STANDARDIZATION.md`).
  Native HR and the conversion basis are retained per cell. Categorical exposures
  are flagged and kept off the comparable scale.
- Each evidence cell carries a `verification_status` flag. **Verify all figures
  against primary sources before clinical use.**

To recompute the standardized values after editing native HRs:

```sh
node tools/standardize.js
```

## Moving this to its own repository

This folder is self-contained:

```sh
cp -r metabolic-research ~/longevity-heatmap
cd ~/longevity-heatmap && git init && git add . && git commit -m "Initial commit"
```
