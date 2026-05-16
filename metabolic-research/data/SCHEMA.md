# Data Layer Schema — `biomarkers.json`

The data layer is deliberately kept separate from the visualization so biomarkers
can be added, edited, or re-cited without touching the heat map code.

## Top-level structure

```json
{
  "meta": { ... },
  "biomarkers": [ { ... }, ... ]
}
```

## `meta`

```json
{
  "title": "Longevity Biomarker Heat Map",
  "version": "0.x.0",
  "last_updated": "YYYY-MM-DD",
  "outcomes": ["all_cause_mortality", "cvd", "cancer", "dementia", "frailty"],
  "evidence_tiers": {
    "A": "Multiple meta-analyses, or meta-analysis plus Mendelian randomization",
    "B": "Single meta-analysis, or a strong large prospective cohort",
    "C": "Limited evidence: small, purely observational, or conflicting"
  },
  "domains": [
    { "id": "lipids_cv", "name": "Lipids / Cardiovascular Risk", "order": 2 }
  ]
}
```

## `biomarkers[]`

Each biomarker is one object:

```json
{
  "id": "apob",                       // unique, snake_case
  "name": "ApoB",                     // display name
  "domain": "lipids_cv",              // matches meta.domains[].id
  "aliases": ["apolipoprotein B"],
  "profile": {
    "mechanism": "Why it matters for aging/mortality (1-3 sentences).",
    "measurement": {
      "method": "How it is measured.",
      "sample_type": "venous blood | capillary blood | saliva | urine | hair | in-office device | wearable | imaging | functional test",
      "approx_cost_usd": "15-40"      // string range, USD
    },
    "signal_quality": "clean | moderate | noisy",
    "signal_quality_note": "Brief rationale.",
    "confounders": ["..."],
    "testing_cadence": "e.g. 'baseline + annual', 'quarterly'"
  },
  "outcomes": {
    "all_cause_mortality": { <cell> | null },
    "cvd":                 { <cell> | null },
    "cancer":              { <cell> | null },
    "dementia":            { <cell> | null },
    "frailty":             { <cell> | null }
  },
  "modifiability": {
    "intervention": "Best-evidence intervention (lifestyle / pharma / supplement).",
    "effect_size": "Magnitude of change achievable.",
    "timeframe": "Time to a meaningful change.",
    "evidence_type": "RCT | RCT meta-analysis | observational | mechanistic",
    "citation": "Source.",
    "rating": "high | moderate | low | fixed"
  }
}
```

`null` for an outcome means no usable evidence was found — the heat map renders
it as an explicit "no data" cell rather than implying a null effect.

## `<cell>` — one biomarker × one outcome

```json
{
  "hr": 1.43,
  "hr_metric": "per SD | per quartile | per unit | categorical",
  "hr_unit_detail": "Exact contrast, e.g. 'per 1 SD higher', 'per 5 kg lower', 'Q4 vs Q1'.",
  "direction": "risk | protective",   // risk: higher value -> worse outcome
  "ci_low": 1.35,
  "ci_high": 1.52,
  "evidence_tier": "A | B | C",
  "n": 350000,                        // number or descriptive string
  "followup_years": 12.0,
  "study": "First author Year, Journal, design (short).",
  "study_type": "meta-analysis | mendelian-randomization | prospective-cohort | pooled-cohort | RCT",
  "citation": "Full citation with DOI or PMID where possible.",
  "verification_status": "verified | approximate | unverified",
  "notes": "Caveats: unit conversions, association vs causal, conflicting data."
}
```

### `verification_status`
- `verified` — the cited study was located and the HR / CI / N / follow-up confirmed against it.
- `approximate` — the study is correctly identified but exact CI / N / follow-up could not be confirmed.
- `unverified` — drawn from general domain knowledge; the specific study was not confirmed.

Cells must never carry invented precision. If a number cannot be sourced, leave
the field `null` and explain in `notes`, rather than guessing.

## Color / encoding conventions used by the heat map
- HR magnitude is rendered as `|ln(HR)|` so risk and protective markers of equal
  strength read as equal intensity; `direction` drives the hue.
- `evidence_tier` drives cell border weight / icon (A strongest).
- View A additionally buckets all-cause-mortality HR into effect-size columns.
