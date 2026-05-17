# Longevity Biomarker Heat Map

A rigorous, citation-backed map of how testable biomarkers predict aging-related
outcomes — across 14 domains, on one comparable scale, with an
intervention-modifiability layer and a personal scoring dashboard.

> ⚠️ **This is a research and educational resource — not medical advice and not
> a diagnostic tool.** See [`DISCLAIMER.md`](DISCLAIMER.md) before using it.

## What this is

Most consumer longevity panels report a biomarker against a static "optimal
range." This project does something different: it ties every biomarker to the
**published hazard ratios** linking it to aging-related outcomes, puts those
hazard ratios on **one comparable scale**, and is honest about **evidence
quality** and **where the field oversells**.

- **135 biomarkers**, 14 domains, **519 evidence cells** across five outcomes
  (all-cause mortality, cardiovascular disease, cancer, dementia, frailty).
- Hazard ratios are **standardized to HR per +1 SD** so markers measured in
  different units are directly comparable — see
  [`data/HR_STANDARDIZATION.md`](data/HR_STANDARDIZATION.md).
- Every quantified cell was **checked against primary sources** in two passes:
  **259 verified, 97 approximate, 0 unverified**. Cells that could not be
  anchored to a real source were left blank, not guessed.
- Each cell carries an **A/B/C evidence tier** and a modifiability rating.

## The two tools

Both are server-free — open the HTML file in any browser.

**The heat map** — `heatmap/index.html` (or the single-file
`heatmap/heatmap-standalone.html`). The population evidence map: 135 biomarkers
× 5 outcomes, coloured by standardized effect size, with the all-cause-mortality
column highlighted and a modifiability column. Hover any cell for the source;
click a biomarker for its profile.

**The personal dashboard** — `client/dashboard.html`. Scores *your own*
measurements against the map. Put your biomarker panels into
`client/measurements.js`; each value is placed on the population SD scale,
translated through the verified per-SD hazard ratios into modeled risk, and
triaged by predictive weight × modifiability × your distance from optimal. A
second panel adds a risk trajectory. There is deliberately **no single
composite score** — correlated hazard ratios cannot be multiplied together.

## Repository layout

```
├── data/
│   ├── biomarkers.json        the evidence layer — biomarkers, cells, modifiability
│   ├── reference_ranges.json  population distribution + risk-optimal range per biomarker
│   ├── SCHEMA.md              data-layer schema
│   └── HR_STANDARDIZATION.md  how hazard ratios are converted to one per-SD scale
├── tools/
│   ├── standardize.js         deterministic HR → per-SD converter (idempotent)
│   ├── audit.js               internal QA audit of the data layer
│   ├── build-standalone.js    bundles the heat map into one portable HTML file
│   └── build-client.js        bundles data + reference ranges for the dashboard
├── heatmap/                   the population evidence map
├── client/                    the personal coaching dashboard
├── profiles/biomarker-profiles.md   per-biomarker profiles + modifiability
├── writeup/summary.md         14-domain synthesis
├── DISCLAIMER.md   ·   LICENSE   ·   LICENSE-DATA
```

The data layer is kept separate from the visualizations so biomarkers can be
added or re-cited without touching code.

## Regenerating

After editing `data/biomarkers.json` or `data/reference_ranges.json`:

```sh
node tools/standardize.js       # per-SD values + heatmap/data.bundle.js
node tools/build-standalone.js  # heatmap/heatmap-standalone.html
node tools/build-client.js      # client/refdata.js
node tools/audit.js             # QA check (expect 0 errors)
```

## Evidence standards

- Hazard ratios prioritize meta-analyses, Mendelian randomization, and large
  prospective cohorts (UK Biobank, ARIC, Framingham, MESA, Whitehall II, etc.).
- The native HR and the conversion basis are retained per cell; categorical
  exposures (genotypes, U-shaped markers, clinical strata) are flagged and kept
  off the comparable scale rather than forced onto it.
- `verification_status` is honest. `approximate` and `unverified` figures — and
  the `approximate`/`unverified` reference ranges — are provisional; verify
  against the cited primary sources before any serious use.

## Corrections welcome

This is an evidence base, and evidence bases improve by being checked. If you
find a wrong hazard ratio, a misattributed citation, or a better source, please
open an issue or a pull request — cite the primary source.

## Licensing

- **Source code** (`tools/`, and the scripts/markup in `heatmap/` and
  `client/`): MIT — see [`LICENSE`](LICENSE).
- **Dataset and written content** (`data/`, `profiles/`, `writeup/`): Creative
  Commons Attribution 4.0 (CC BY 4.0) — see [`LICENSE-DATA`](LICENSE-DATA).

Re-use is welcome, including commercially; please attribute.

## Citing

> Mislinski, J. Longevity Biomarker Heat Map: a citation-backed,
> per-SD-standardized map of biomarker–outcome hazard ratios across 14 domains.
> 2026.

---

*Built as both a defensible academic artifact and the evidence engine for a
metabolic longevity assessment practice. Not medical advice — see
[`DISCLAIMER.md`](DISCLAIMER.md).*
