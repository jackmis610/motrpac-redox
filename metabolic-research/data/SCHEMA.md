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
  "comparability_tag": "Genotype",    // optional — see below
  "comparability_note": "...",        // optional — see below
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

### `comparability_tag` / `comparability_note` (optional)
Set these only on biomarkers that have one or more `categorical` or
`unconvertible` cells. `comparability_tag` is a short label (≤ ~18 chars, e.g.
`Genotype`, `U-shaped`, `Non-normal`, `Clinical cutpoint`) shown in the heat
map's Standardization column; `comparability_note` is a one-sentence
explanation, shown on hover and in the biomarker profile. Biomarkers whose
cells are all comparable omit both fields.

## `<cell>` — one biomarker × one outcome

```json
{
  "hr": 1.43,                         // NATIVE HR, exactly as the source study reported it
  "hr_metric": "per_sd | per_log_sd | tertile | quartile | quintile | decile | median_split | per_unit | per_doubling | categorical",
  "hr_unit_detail": "Exact contrast, e.g. 'per 1 SD higher', 'per 5 kg lower', 'Q4 vs Q1'.",
  "direction": "risk | protective",   // risk: higher value -> worse outcome
  "ci_low": 1.35,                     // native CI
  "ci_high": 1.52,

  "hr_per_sd": 1.16,                  // STANDARDIZED: HR per +1 SD. null if categorical/unconvertible.
  "hr_per_sd_ci": [1.10, 1.23],       // standardized CI. null if unavailable.
  "comparability": "native | converted | categorical | unconvertible",
  "hr_per_sd_basis": "How hr_per_sd was derived: 'native', or conversion formula + inputs.",
  "conversion": {                     // inputs the standardization script uses; omit/null when not needed
    "unit_value": 10,                 // per_unit: numeric size of the reported unit
    "population_sd": 35,              // per_unit: biomarker SD in the SAME units as unit_value
    "log2_sd": 0.90,                  // per_doubling: SD of log2(marker)
    "sd_source": "Source of the SD (study name or referenced population)."
  },

  "evidence_tier": "A | B | C",
  "n": 350000,                        // number or descriptive string
  "followup_years": 12.0,
  "study": "First author Year, Journal, design (short).",
  "study_type": "meta-analysis | mendelian-randomization | prospective-cohort | pooled-cohort | RCT",
  "citation": "Full citation with DOI or PMID where possible.",
  "verification_status": "verified | approximate | unverified",
  "notes": "Caveats: association vs causal, conflicting data, non-linearity."
}
```

### Standardization fields
`hr` / `ci_low` / `ci_high` always hold the **native** figures as published.
`hr_per_sd` holds the **comparable** figure (HR per +1 SD) and is what the heat
map's magnitude scale uses. See `HR_STANDARDIZATION.md` for the conversion rules.
`hr_per_sd` is computed by a deterministic script from `hr_metric` + `conversion`,
never entered by hand. Quantile metrics need no `conversion` inputs (the SD
spread is fixed); `per_unit` and `per_doubling` require it.

### `verification_status`
- `verified` — the cited study was located and the HR / CI / N / follow-up confirmed against it.
- `approximate` — the study is correctly identified but exact CI / N / follow-up could not be confirmed.
- `unverified` — drawn from general domain knowledge; the specific study was not confirmed.

Cells must never carry invented precision. If a number cannot be sourced, leave
the field `null` and explain in `notes`, rather than guessing.

## Color / encoding conventions used by the heat map
- Magnitude is rendered as `|ln(hr_per_sd)|` so all markers sit on one comparable
  scale; `direction` drives the hue.
- `categorical` and `unconvertible` cells have no `hr_per_sd`; they are shown in a
  separate band and never bucketed on the comparable scale.
- `evidence_tier` drives cell border weight / icon (A strongest).
- View A additionally buckets all-cause-mortality `hr_per_sd` into effect-size columns.
