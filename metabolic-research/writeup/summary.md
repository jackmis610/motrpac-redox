# Longevity Biomarker Heat Map — Synthesis

**Interim writeup — domains 2–4 of 14.** This covers Lipids / Cardiovascular
Risk, Mitochondrial / Cardiorespiratory Fitness, and Inflammation (38 biomarkers,
80 evidence cells). It will be revised into a full synthesis once the remaining
domains are researched. Conclusions below are scoped to these three domains.

---

## 1. Which biomarkers carry the most predictive weight

Three clusters dominate across these domains.

**Atherogenic particle burden (causal, Tier A).** ApoB, non-HDL cholesterol, and
remnant cholesterol show concordant meta-analytic *and* Mendelian-randomization
support for a causal role in cardiovascular disease (apoB CVD HR ≈ 1.6 per SD;
remnant cholesterol HR ≈ 2.8 in top-vs-bottom contrasts). They are the most
defensible biomarkers in the set: the association is large, causal, and
intervention-confirmed. LDL-C carries the same causal weight but is a noisier
*measurement* of that biology than apoB.

**Imaging of established disease (Tier A, very strong but fixed).** Coronary
artery calcium is the single strongest predictor here (CVD HR ≈ 6.8, all-cause
mortality HR ≈ 2.6 for high vs zero scores). It measures disease that has
*already* occurred — extremely informative for stratification, useless as a
coaching target.

**Physical function and fitness (Tier A, broad).** VO2max, grip strength, and
gait speed predict all-cause mortality with Tier A evidence and unusually broad
reach — the same markers also track frailty and dementia. Their breadth across
outcomes is their distinguishing feature; no lipid or inflammatory marker
predicts as many outcomes at once.

**Inflammation as a wide-but-shallow signal.** hsCRP, IL-6, and fibrinogen each
predict all-cause mortality, CVD, *and* dementia. The signal is real and broad,
but mostly modest per-SD and — with one exception — non-causal (see §2).

> **Units caveat.** Hazard ratios here are not on one scale. Categorical
> extreme-group contrasts (e.g. the sit-to-stand test, lowest vs highest group,
> HR ≈ 3.8–6.0) look dramatically larger than per-SD markers but are not directly
> comparable. The heat map records the exact metric per cell; do not rank
> biomarkers on raw HR without checking it.

## 2. Where the field oversells

- **hsCRP is a predictor, not a cause.** It is cheap and broadly predictive, but
  Mendelian randomization does not support a causal role in CVD or cancer.
  Lowering CRP itself does nothing; CRP improves because the upstream driver
  (adiposity, inactivity, IL-6 signalling) improved. Treat it as a thermometer.
- **IL-6 is the genuine exception.** *IL6R* genetic variants and the CANTOS trial
  support a causal inflammatory pathway. Of the inflammatory markers, IL-6 is the
  one with a mechanistic mandate.
- **Homocysteine.** Robust observational associations with mortality and
  dementia — but B-vitamin trials that successfully lowered homocysteine did not
  reduce events. Predictive, not actionable as a target.
- **Fibrinogen and ferritin.** Fibrinogen predicts but MR does not support
  causality. Ferritin as an "inflammation marker" is badly confounded by iron
  status and the acute-phase response; the signal attenuates after adjustment.
- **VO2max causality.** The VO2max–mortality association is large and consistent
  but observational; two-sample MR has not established a causal effect on
  longevity, and reverse causation (subclinical disease lowering fitness) is real.
  This does not argue against training — it argues against treating a single
  VO2max number as a causal lever rather than a fitness *state*.
- **Lipid sub-fractionation.** LDL particle size and oxidized LDL add little
  independent value over apoB / LDL-P and suffer assay heterogeneity. HDL-C is
  non-causal (MR plus neutral HDL-raising trials); HDL efflux capacity is
  promising but has no standardized commercial assay.
- **The "physiology" markers.** Resting metabolic rate, RER / fat oxidation,
  exercise economy, lactate threshold, ventilatory thresholds, and VO2 kinetics
  have essentially no hard-outcome hazard ratios in aging populations. They are
  valuable for *training prescription*; they should not be sold to clients as
  longevity biomarkers.

## 3. Highest-leverage modifiable markers

The coaching sweet spot is markers that are both predictive *and* responsive to
intervention. From these three domains, the priority targets are:

| Marker | Predictive role | Intervention leverage |
|---|---|---|
| VO2max | Tier A, broad outcomes | +10–25% with structured training, 8–12 wk |
| Grip / lower-body strength | Tier A mortality + frailty | Resistance training, measurable in 8–12 wk |
| ApoB / non-HDL / LDL-P | Tier A causal CVD | Statin ± ezetimibe ± PCSK9: 35–60%, 4–6 wk |
| Triglycerides | Tier B–C, metabolic | Diet, weight loss, fitness — fast-responding |
| hsCRP | Broad predictor | Falls with weight loss, exercise, statins |

Note the asymmetry: lipid markers are highly modifiable *pharmacologically*;
fitness and strength are highly modifiable *behaviourally*. hsCRP is modifiable
only indirectly (it follows its upstream drivers).

**Predictive but fixed — use for stratification, not coaching:** Lp(a) and ApoE
genotype (genetically set), and CAC (established disease that does not regress).
These belong in the baseline risk picture, not the quarterly scorecard.

## 4. Gaps worth flagging to clients

- **Cancer prediction is weak** across all three domains — lipids and fitness
  carry essentially no cancer hazard ratios, and inflammatory markers are modest.
  Cancer risk in this protocol will rest on other domains and on screening.
- **Frailty evidence is largely cross-sectional.** Strength and fitness correlate
  with frailty, but longitudinal *prediction* is weaker, and several functional
  markers are partly definitional components of the frailty phenotype
  (circularity).
- **Causal vs observational gap.** Much of the fitness and inflammation evidence
  is observational. The heat map's A/B/C tiers and per-cell notes flag this; the
  distinction should be made explicit with clients.
- **Unit inconsistency** (per SD vs per quartile vs categorical) limits direct
  cross-marker comparison and is an inherent property of the literature.
- **Newer markers** (suPAR, GlycA, HDL efflux capacity) are mechanistically
  interesting but rest on thin, few-cohort evidence — reasonable to measure, not
  yet reasonable to weight heavily.

## 5. Mapping to an assessment protocol

A defensible protocol tiers biomarkers by what they are *for*:

- **Causal risk panel (act on these):** apoB or non-HDL-C, triglycerides, hsCRP.
  These are cheap, causal or near-causal, and modifiable.
- **One-time genetic / structural stratifiers:** Lp(a) (once in a lifetime),
  ApoE genotype (once), CAC (baseline, age ≈ 40+). These set the risk backdrop
  and rarely need repeating.
- **Functional core (the longevity dashboard):** VO2max, grip strength,
  lower-body strength/power, gait speed, heart-rate recovery, HRV. Broad-outcome,
  behaviourally modifiable, and the natural extension of your existing CRF model.
- **Context / interpret-with-caution:** IL-6, fibrinogen, homocysteine, ferritin,
  NLR, LDL-P, cIMT — measure where it sharpens a specific question; do not
  over-weight.

## 6. Suggested testing cadence

| Tier | Markers | Cadence |
|---|---|---|
| Intake (baseline) | Full panel: causal lipids, Lp(a), ApoE, CAC, hsCRP, full functional battery, cIMT | Once at onboarding |
| Quarterly | apoB / non-HDL, triglycerides, hsCRP, VO2max (if actively training), grip strength | Every 3 months |
| Annual | Full functional battery, repeat lipid + inflammatory panel, IL-6 / fibrinogen if elevated at baseline | Yearly |
| Rarely / never repeat | Lp(a), ApoE genotype | Once; CAC only on a multi-year interval if clinically indicated |

The quarterly set is deliberately small: only markers that move on a 3-month
timescale and that the client can influence. Fixed stratifiers and slow-moving
functional measures belong on the annual or one-time cycle.

---

*Generated for the framework-first review. Verify all hazard ratios against the
primary sources cited in `data/biomarkers.json` before clinical use.*
