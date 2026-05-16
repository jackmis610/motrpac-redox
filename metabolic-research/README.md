# Longevity Biomarker Heat Map

A rigorous, citation-backed heat map of longevity biomarkers across multiple
framings, with an intervention-modifiability layer. Built as both a defensible
academic artifact and a clinical coaching tool for a metabolic longevity
assessment protocol.

## Status

**Framework-first build.** The data structure, the full 3-view interactive
visualization, and three fully-researched domains are complete and under
review. Remaining domains are populated in subsequent passes.

| Domain | Status |
|---|---|
| 2 · Lipids / Cardiovascular Risk | researched |
| 3 · Mitochondrial / Cardiorespiratory Fitness | researched |
| 4 · Inflammation | researched |
| 1, 5–14 (glucose, clocks, hormones, body comp, etc.) | pending |

## Layout

```
metabolic-research/
├── data/
│   ├── biomarkers.json        the data layer — biomarkers, evidence cells, modifiability
│   ├── SCHEMA.md              field-by-field schema documentation
│   └── HR_STANDARDIZATION.md  how hazard ratios are converted to a common per-SD scale
├── tools/
│   └── standardize.js         deterministic HR → per-SD converter (idempotent)
├── heatmap/
│   ├── index.html             interactive 3-view heat map
│   ├── styles.css
│   └── app.js                 vanilla JS, no framework, no build step
├── profiles/
│   └── biomarker-profiles.md  Deliverable 1 — per-biomarker profiles
└── writeup/
    └── summary.md             Deliverable 5 — synthesis writeup
```

The data layer (`data/biomarkers.json`) is deliberately separate from the
visualization so biomarkers can be added or re-cited without touching the
heat map code. See `data/SCHEMA.md` for the contract.

## Running the heat map

Open `heatmap/index.html` directly in any browser — no web server needed. The
data layer is bundled into `heatmap/data.bundle.js` so it loads from `file://`.

`data/biomarkers.json` remains the editable source of truth. After changing it,
regenerate the bundle and the per-SD standardized values:

```sh
node tools/standardize.js
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
