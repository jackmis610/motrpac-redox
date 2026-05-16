# HR Standardization Methodology

The source literature reports hazard ratios in incompatible units — per standard
deviation, per quartile contrast, per fixed unit, per doubling. To make the heat
map comparable, every continuous-exposure HR is standardized to a common metric:

> **`hr_per_sd` — the hazard ratio per +1 standard deviation of the biomarker.**

Categorical exposures (genotypes, dichotomous tests, disease-stage categories)
cannot be placed on a per-SD scale and are flagged separately, never mixed into
the comparable analysis.

## The `comparability` flag

Every evidence cell carries one of:

| Value | Meaning |
|---|---|
| `native` | The source study reported the HR per SD (or per SD of the log-transformed marker). Used as-is. |
| `converted` | A continuous-exposure HR reported in another unit, converted to per-SD by the rules below. The exact basis is logged in `hr_per_sd_basis`. |
| `categorical` | The exposure is inherently categorical (e.g. ApoE genotype, CAC score strata, "able vs unable"). `hr_per_sd` is `null`; the native HR is retained and shown in a separate band. |
| `unconvertible` | A continuous exposure whose HR could not be converted (missing SD, non-linear with no usable summary). `hr_per_sd` is `null`; flagged and excluded from the comparable main analysis. |

## Conversion rules

`hr_per_sd` is always a monotone power transform of the native HR:
`hr_per_sd = hr ^ k`, and the same exponent `k` is applied to both confidence
limits. This is exact for the log-HR (it rescales the Cox coefficient β by `k`).

### 1. Native per-SD — `k = 1`
Reported per SD, or per SD of the log-transformed marker (standard for skewed
markers such as CRP, IL-6, triglycerides, ferritin). For skewed markers
`hr_per_sd` is therefore "per SD of log marker"; this is the accepted comparable
quantity and is labelled as such.

### 2. Extreme-quantile contrasts — `k = 1 / spread`
A "top vs bottom group" HR spans more than one SD. For an approximately normal
exposure, the mean separation between extreme groups (in SD units) is a fixed
order-statistic of the standard normal:

| Contrast | Mean separation (SD) | `k` |
|---|---|---|
| Above vs below median | 1.60 | 1 / 1.60 |
| Tertile 3 vs 1 | 2.18 | 1 / 2.18 |
| Quartile 4 vs 1 | 2.54 | 1 / 2.54 |
| Quintile 5 vs 1 | 2.80 | 1 / 2.80 |
| Decile 10 vs 1 | 3.51 | 1 / 3.51 |

Derivation (quartiles): the mean of a standard-normal tail above the 75th
percentile is `φ(0.6745) / 0.25 = 1.271`; by symmetry the separation between
top and bottom quartile means is `2 × 1.271 = 2.54`. The other rows use the same
construction at the relevant cut-points.

**Assumption logged per cell:** this presumes the biomarker is approximately
normal (or normal after the transformation the study applied). Native per-SD
estimates are always preferred over converted quantile estimates; the
`comparability` flag preserves that distinction.

### 3. Per fixed unit — `k = population_SD / unit`
An HR reported "per X units" is converted with the biomarker's population SD:
`hr_per_sd = hr ^ (SD / X)`. The SD value and its source (the source study, or a
referenced population such as NHANES / UK Biobank) are logged in
`hr_per_sd_basis`. If no defensible SD is available the cell is `unconvertible`.

### 4. Per doubling — `k = SD_of_log2`
For markers analyzed on a log2 scale, "per doubling" is per 1 unit of log2.
`hr_per_sd = hr ^ (SD of log2 marker)`. The log2 SD and its source are logged.

## Worked example

GlycA, all-cause mortality, reported as quartile 4 vs 1, HR 1.62 (95% CI
1.45–1.81):

```
k          = 1 / 2.54 = 0.3937
hr_per_sd  = 1.62 ^ 0.3937 = 1.207
ci_per_sd  = [1.45 ^ 0.3937, 1.81 ^ 0.3937] = [1.158, 1.258]
comparability   = "converted"
hr_per_sd_basis = "Q4-vs-Q1 contrast; k=1/2.54 (normal-approx extreme-quartile
                   SD separation); native HR 1.62 (1.45-1.81)"
```

## Auditing

The native `hr`, `hr_metric`, and `hr_unit_detail` are always retained alongside
`hr_per_sd`, so any conversion can be checked. `hr_per_sd` values are computed by
a single deterministic script from the logged inputs — not entered by hand — so
the arithmetic is reproducible. The heat map's comparable views (bucketing,
colour intensity) use `hr_per_sd`; tooltips show both the native and the
standardized figure.
