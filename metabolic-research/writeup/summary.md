# Longevity Biomarker Heat Map — Synthesis

A 14-domain, 135-biomarker map of how testable biomarkers predict aging-related
outcomes, with an intervention-modifiability layer. 519 evidence cells: 356
quantified (241 standardized to a comparable HR-per-SD scale, 100 categorical,
15 unconvertible) and 163 documented-but-unquantified.

Every quantified cell has been checked against primary sources, twice: a
verification pass followed by a re-sourcing pass for anything that could not be
confirmed. The result is **257 verified** (source located, HR and CI confirmed)
and **99 approximate** (study identified, exact CI/N not confirmable) — and
**zero unverified**. Cells whose stated figure could not be anchored to a real
primary source were nulled rather than retained: ~60 such cells are now
documented-but-unquantified. No displayed hazard ratio is unsourced.

---

## 1. Which biomarkers carry the most predictive weight

Across all 14 domains, the heaviest, most defensible predictors cluster into six
groups.

**Atherogenic particle burden — causal, Tier A.** ApoB, non-HDL cholesterol, and
remnant cholesterol carry concordant meta-analytic *and* Mendelian-randomization
support for a causal CVD role (per-SD HR ≈ 1.4–1.7). The single most defensible
actionable cluster in the map.

**Blood pressure — causal, Tier A, broad.** Resting and ambulatory blood
pressure are among the strongest CVD predictors here (per-SD HR ≈ 1.3–1.9) and
also track dementia. Causal, cheap, modifiable — they belong beside the lipids.

**Physical function and fitness — Tier A, broadest reach.** VO2max, grip
strength, gait speed, the sit-to-stand test, and appendicular lean mass predict
all-cause mortality, and uniquely also reach frailty and dementia. No blood
marker predicts as many outcomes at once.

**Renal markers — Tier A, under-rated.** Albuminuria and cystatin C-based eGFR
are strong, broad predictors (albuminuria CVD HR/SD ≈ 1.6) and are routinely
under-weighted in consumer longevity panels.

**Inflammation — wide but mostly non-causal.** hsCRP, IL-6, and fibrinogen each
predict mortality, CVD, and dementia (fibrinogen CVD HR/SD ≈ 1.9) — a real,
broad signal, but with one causal exception (§2).

**Integrative second/third-generation biological-age clocks.** GrimAge and
DunedinPACE predict mortality strongly (DunedinPACE HR/SD ≈ 1.65) and survive
lifestyle adjustment — useful as a composite, with caveats (§2).

Glucose metabolism (HbA1c, fasting insulin, HOMA-IR), visceral adipose tissue,
processing speed, polygenic risk scores, and pTau-217 (for dementia
specifically) round out the high-weight set.

## 2. Where the field oversells

- **hsCRP and fibrinogen** are strong predictors but **not causal** (Mendelian
  randomization is null). **IL-6 is the exception** — *IL6R* variants and CANTOS
  support a causal inflammatory pathway. Lowering CRP itself does nothing.
- **Biological-age clocks.** First-generation clocks (Horvath, Hannum) predict
  mortality far more weakly than GrimAge/DunedinPACE and were non-significant in
  several cohorts. **GlycanAge** has essentially no prospective per-SD mortality
  evidence — the clearest over-selling risk. **Telomere length** is measurement-
  noisy with small, inconsistent HRs. Even GrimAge's CVD/cancer signal is
  substantially a built-in smoking surrogate.
- **Micronutrients — the largest observational-vs-trial gap in the dataset.**
  Vitamin D, selenium, omega-3, and B-vitamins show graded observational
  associations but their **supplementation RCTs are overwhelmingly null** for
  hard outcomes (VITAL, D-Health, SELECT, the B-Vitamin Treatment Trialists).
  Several "risk" signals (high B12, high copper) are reverse causation flagging
  occult disease.
- **Homocysteine** predicts events but B-vitamin lowering trials did not reduce
  them — predictive, not actionable.
- **CGM-derived metrics, wearable sleep metrics** have no hard-outcome cohort
  evidence in non-diabetic / general populations — promising, not yet evidential.
- **Hormone "optimization."** Most hormonal axes have U-shaped or strongly
  sex-specific mortality relationships; growth hormone is pulsatile and its
  evidence comes from disease states. A single hormone value is rarely a clean
  per-SD lever.
- **Pharmacogenomic panels** predict drug response, not aging outcomes — valuable
  for medication safety, miscategorised if sold as longevity biomarkers.
- **VO2max** causality: the mortality association is large but observational; MR
  has not established a causal longevity effect. Train for fitness as a *state*,
  not to move a number. **HDL-C** is non-causal and U-shaped. **LDL particle
  size and oxidized LDL** add little over apoB.

## 3. Highest-leverage modifiable markers

The coaching priorities — high predictive weight *and* high modifiability — are,
in descending standardized strength:

| Cluster | Markers | Leverage |
|---|---|---|
| Atherogenic lipids | apoB, non-HDL-C, LDL-C, LDL-P, triglycerides, remnant-C | Statin ± ezetimibe ± PCSK9: 35–60%, 4–6 wk (pharmacological) |
| Blood pressure | resting & ambulatory BP | Lifestyle + antihypertensives; weeks |
| Fitness & strength | VO2max, grip, sit-to-stand, lower-body power | Aerobic + resistance training; 8–12 wk (behavioural) |
| Central adiposity / glucose | visceral adipose tissue, HbA1c, fasting insulin | Weight loss, diet, fitness; months |
| Renal | albuminuria | Responds to BP/glucose control and SGLT2 inhibitors |
| Bone | bone mineral density | Resistance + impact training, modest; months–years |

Note the asymmetry: lipids and BP are highly modifiable *pharmacologically*;
fitness, strength, and adiposity *behaviourally*. hsCRP moves only indirectly
(it follows its upstream drivers) — a false-priority caution: as with
homocysteine, "modifiable" must mean the *outcome* moves, not just the marker.

**Predictive but fixed — stratification, not coaching:** Lp(a), ApoE genotype,
MTHFR, polygenic risk scores, coronary artery calcium. These set the baseline
risk picture and rarely need repeating.

## 4. Gaps worth flagging to clients

- **Cancer prediction is weak** across almost every domain except polygenic risk
  scores and parts of inflammation. Cancer risk rests on screening, not on this
  biomarker map.
- **Frailty evidence is largely cross-sectional or definitional** — several
  functional markers are themselves components of the frailty phenotype.
- **Sleep has no native per-SD evidence at all** — the entire domain is reported
  in clinical strata or U-shaped curves; informative but not standardizable.
- **Emerging blood markers** (pTau-217, NfL, GFAP) are strong for pathology but
  their incident-outcome HR evidence is young; **CGM and wearable metrics** lack
  outcome cohorts entirely.
- **Causal vs observational** — much of the fitness, inflammation, hormonal, and
  micronutrient evidence is observational; the A/B/C tiers and per-cell notes
  flag this and it should be explicit with clients.
- **Comparability** — converted quantile estimates assume approximate normality;
  79 categorical cells (genotypes, U-shaped markers, clinical strata) sit
  off-scale by necessity, flagged in the Standardization column.

## 5. Mapping to an assessment protocol

A defensible protocol tiers biomarkers by *what they are for*:

- **Causal risk panel (act on these):** apoB / non-HDL-C, triglycerides, blood
  pressure, HbA1c + fasting insulin, hsCRP, albuminuria, cystatin C. Cheap,
  causal or near-causal, modifiable.
- **Functional core (the longevity dashboard):** VO2max, grip strength, gait
  speed / sit-to-stand, lower-body power, appendicular lean mass — broad-outcome
  and behaviourally modifiable; the natural extension of a CRF-based practice.
- **One-time genetic / structural stratifiers:** Lp(a), ApoE genotype, polygenic
  risk scores, pharmacogenomic panel, coronary artery calcium, DEXA (visceral
  fat, bone density).
- **Integrative annual check:** a second/third-generation epigenetic clock
  (GrimAge or DunedinPACE) as a composite — interpreted with its caveats.
- **Context / interpret with caution:** hormonal axes, micronutrients, CGM
  metrics, sleep architecture, first-generation clocks, GlycanAge — measure
  where they sharpen a specific question; do not over-weight.

## 6. Suggested testing cadence

| Tier | Markers | Cadence |
|---|---|---|
| Intake (baseline) | Full panel: causal lipids, Lp(a), ApoE, polygenic + pharmacogenomic panels, CAC, DEXA, blood pressure, full functional battery, hsCRP, renal panel, hormonal + micronutrient screen, an epigenetic clock | Once at onboarding |
| Quarterly | apoB / non-HDL, triglycerides, blood pressure, HbA1c, hsCRP, weight / waist, VO2max (if actively training) | Every 3 months |
| Annual | Full functional battery, repeat lipid / renal / inflammatory / glucose panels, epigenetic clock, body composition | Yearly |
| Rarely / never repeat | Lp(a), ApoE, polygenic and pharmacogenomic panels; CAC only on a multi-year interval if clinically indicated | Once |

The quarterly set is deliberately small — only markers that move on a 3-month
timescale and that the client can influence. Fixed stratifiers and slow-moving
measures belong on the annual or one-time cycle.

---

*Generated from `data/biomarkers.json` (v1.0.0). Hazard ratios are standardized
per `data/HR_STANDARDIZATION.md`; every figure is traceable to a cited source in
the data layer and should be verified against the primary literature before
clinical use.*
