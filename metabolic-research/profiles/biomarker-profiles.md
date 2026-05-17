# Biomarker Profiles

Deliverable 1 (profiles) and Deliverable 3 (modifiability layer), generated from `data/biomarkers.json` by `tools/`.
All 14 domains.

Hazard ratios in the outcome line are standardized to **HR per +1 SD** where the
exposure is continuous; categorical exposures show the native HR (see `data/HR_STANDARDIZATION.md`).


## Insulin Sensitivity / Glucose Metabolism

### Fasting glucose

*Also: fasting plasma glucose, FPG, fasting blood sugar*

- **Mechanism.** Fasting glucose reflects hepatic glucose output and basal insulin action; chronic elevation drives advanced glycation, endothelial dysfunction, beta-cell stress and microvascular damage. It is the diagnostic anchor for prediabetes and diabetes, conditions that shorten lifespan and accelerate vascular and cognitive aging.
- **Measurement.** Enzymatic hexokinase or glucose-oxidase assay on plasma after an 8-12 h fast; point-of-care meters available but less accurate. — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** moderate — Single fasting measurements have meaningful day-to-day biological variability (~5-8% CV) and are sensitive to recent diet, stress, illness and sample handling (glycolysis in uncentrifuged tubes).
- **Major confounders.** acute illness/stress hyperglycaemia; recent carbohydrate intake or incomplete fast; delayed sample processing (glycolysis lowers values); glucocorticoids and other hyperglycaemic drugs; reverse causation: weight loss/frailty lowers glucose
- **Testing cadence.** baseline + annual; more frequent in prediabetes
- **Standardization (J-shaped).** Fasting glucose has a J/U-shaped relation with all-cause mortality — both hypoglycaemic and hyperglycaemic ranges raise risk — so no single per-SD slope represents the mortality cell, which uses a non-monotonic flag.
- **Modifiability (high).** Caloric restriction/weight loss, increased physical activity, low-glycaemic diet; metformin and SGLT2/GLP-1 agents pharmacologically. — *effect:* Intensive lifestyle change lowers fasting glucose ~0.3-0.8 mmol/L and cuts progression to diabetes ~58%; metformin lowers fasting glucose ~1-2 mmol/L in dysglycaemia.; *timeframe:* Weeks to a few months; *evidence:* RCT meta-analysis. Diabetes Prevention Program (Knowler 2002) and Finnish DPS; lifestyle and metformin RCTs.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.11 (tier A); Cancer HR/SD 1.04 (tier B); Dementia HR/SD 1.17 (tier B); Frailty native HR 1.25 [categorical, tier C]

### Fasting insulin

*Also: fasting plasma insulin, FPI*

- **Mechanism.** Fasting insulin reflects basal insulin secretion required to maintain euglycaemia; elevated levels signal compensatory hyperinsulinaemia from insulin resistance, which precedes and drives type 2 diabetes, atherogenic dyslipidaemia, hypertension and adverse vascular and tissue aging.
- **Measurement.** Immunoassay (chemiluminescent/ELISA) on fasting plasma or serum; assays are not fully harmonized across platforms. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** noisy — High within-person variability, pulsatile secretion, poor cross-assay standardization and cross-reactivity with proinsulin make single measurements imprecise; right-skewed distribution requires log transformation.
- **Major confounders.** assay/platform differences (no global standardization); recent food intake or incomplete fast; hepatic insulin clearance variation; acute stress and physical activity; exogenous insulin or insulin secretagogues
- **Testing cadence.** baseline + annual if metabolic risk present
- **Modifiability (high).** Weight loss, aerobic plus resistance exercise, carbohydrate restriction; metformin and GLP-1/SGLT2 agents reduce hyperinsulinaemia indirectly. — *effect:* Lifestyle weight loss can lower fasting insulin 25-40%; exercise training reduces it independent of weight; effect is large and rapid.; *timeframe:* Days to weeks; *evidence:* RCT meta-analysis. Exercise and caloric-restriction RCT meta-analyses of insulin sensitivity.
- **Outcome evidence.** All-cause mortality HR/SD 1.13 (tier C); CVD HR/SD 1.18 (tier B); Cancer HR/SD 1.21 (tier C); Dementia HR/SD 1.07 (tier C); Frailty HR/SD 1.13 (tier C)

### HOMA-IR

*Also: homeostatic model assessment of insulin resistance, HOMA index*

- **Mechanism.** HOMA-IR (fasting glucose x fasting insulin / 22.5) estimates whole-body insulin resistance. Insulin resistance is an upstream driver of type 2 diabetes, atherogenic dyslipidaemia, hypertension, NAFLD and accelerated vascular and cognitive aging.
- **Measurement.** Calculated from paired fasting glucose and fasting insulin; inherits the insulin assay's lack of standardization. — *sample:* venous blood; *approx. cost:* $25-70.
- **Signal quality.** noisy — Combines the high variability and poor standardization of fasting insulin with fasting glucose noise; performs poorly in those with beta-cell failure and is platform-dependent.
- **Major confounders.** insulin assay/platform differences; incomplete fast; hepatic vs peripheral insulin resistance not distinguished; beta-cell failure invalidates the model; acute illness and stress
- **Testing cadence.** baseline + annual if metabolic risk present
- **Modifiability (high).** Weight loss, aerobic + resistance exercise, low-glycaemic/Mediterranean diet; metformin, pioglitazone and GLP-1 agents pharmacologically. — *effect:* Intensive lifestyle change improves HOMA-IR 25-50%; exercise improves it independent of weight loss.; *timeframe:* Weeks to a few months; *evidence:* RCT meta-analysis. Lifestyle and exercise RCT meta-analyses of insulin resistance; Diabetes Prevention Program.
- **Outcome evidence.** All-cause mortality HR/SD 1.16 (tier B); CVD HR/SD 1.19 (tier A); Cancer HR/SD 1.11 (tier C); Dementia HR/SD 1.12 (tier C); Frailty HR/SD 1.14 (tier C)

### HbA1c

*Also: glycated haemoglobin, glycohemoglobin, A1c*

- **Mechanism.** HbA1c reflects the average blood glucose over the preceding 8-12 weeks via non-enzymatic glycation of haemoglobin. It integrates chronic glycaemic exposure that drives micro- and macrovascular damage and is the standard marker for diabetes diagnosis and control.
- **Measurement.** HPLC, immunoassay or enzymatic assay; internationally standardized to the IFCC/NGSP reference; no fasting required. — *sample:* venous blood; *approx. cost:* $10-40.
- **Signal quality.** clean — Well-standardized, low within-person variability and no fasting needed; but biased by conditions altering red-cell lifespan (anaemia, haemoglobinopathies, recent transfusion, pregnancy).
- **Major confounders.** anaemia and iron deficiency (can raise HbA1c); haemolysis/haemoglobinopathies (lower HbA1c); chronic kidney or liver disease; recent blood transfusion or erythropoietin; ethnic differences in glycation
- **Testing cadence.** baseline + annual; every 3-6 months if dysglycaemic or treated
- **Modifiability (high).** Weight loss, low-glycaemic diet, physical activity; metformin, GLP-1 agonists, SGLT2 inhibitors and insulin in diabetes. — *effect:* Lifestyle change lowers HbA1c ~0.3-0.5 percentage points in prediabetes; glucose-lowering drugs lower it 0.5-1.5 points in diabetes.; *timeframe:* 3 months (one HbA1c turnover) for a meaningful change; *evidence:* RCT meta-analysis. Diabetes Prevention Program; glucose-lowering RCT meta-analyses (UKPDS, etc.).
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.10 (tier A); Cancer HR/SD 1.04 (tier C); Dementia HR/SD 1.08 (tier B); Frailty HR/SD 1.14 (tier C)

### OGTT 2-hour glucose

*Also: 2-hour postload glucose, 2-h plasma glucose, oral glucose tolerance test 2h glucose*

- **Mechanism.** Two-hour plasma glucose after a 75 g oral glucose load measures the integrated capacity to dispose of a glucose challenge, capturing postprandial/postload hyperglycaemia and muscle insulin sensitivity. Postload glucose detects dysglycaemia missed by fasting tests and predicts diabetes and vascular disease.
- **Measurement.** Plasma glucose 2 h after a standardized 75 g oral glucose load, after an overnight fast; enzymatic assay. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** moderate — Reproducible enzymatic assay, but the test itself has substantial within-person variability (~16-17% CV) driven by gastric emptying, prior diet and activity; burdensome 2-hour protocol.
- **Major confounders.** prior carbohydrate intake and physical activity; gastric emptying rate variation; test-to-test biological variability; acute illness/stress; incomplete glucose-load ingestion or vomiting
- **Testing cadence.** as indicated for dysglycaemia screening; not routine
- **Standardization (Clinical strata).** 2-hour OGTT glucose HRs are reported across diagnostic strata (normal, impaired glucose tolerance, diabetic range), not as a continuous per-SD slope.
- **Modifiability (high).** Weight loss, physical activity (especially post-meal activity), low-glycaemic diet; acarbose and GLP-1 agents blunt postload glucose. — *effect:* Lifestyle intervention substantially lowers 2-h glucose and reverts impaired glucose tolerance in a large fraction of participants.; *timeframe:* Weeks to months; *evidence:* RCT meta-analysis. Diabetes Prevention Program and Finnish DPS (2-h glucose endpoints); STOP-NIDDM (acarbose).
- **Outcome evidence.** All-cause mortality native HR 1.73 [categorical, tier A]; CVD native HR 1.40 [categorical, tier A]; Cancer native HR 1.26 [categorical, tier C]; Dementia native HR 1.30 [categorical, tier C]

### OGTT 2-hour insulin

*Also: 2-hour postload insulin, postchallenge insulin*

- **Mechanism.** Two-hour postload insulin reflects the magnitude of the insulin response sustained after a glucose challenge; elevated postload insulin indicates compensatory hyperinsulinaemia and insulin resistance, while a delayed/sustained peak signals impaired insulin dynamics.
- **Measurement.** Immunoassay for insulin on plasma drawn 2 h into a 75 g OGTT; not standardized across platforms. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** noisy — Combines OGTT protocol variability with poorly standardized, pulsatile insulin measurement; postload insulin interpretation depends on the concurrent glucose, limiting standalone use.
- **Major confounders.** insulin assay/platform differences; OGTT protocol and prior diet/activity; concurrent glucose level needed for interpretation; hepatic insulin clearance variation; gastric emptying rate
- **Testing cadence.** research/specialist use; not routine
- **Modifiability (high).** Weight loss, exercise training, carbohydrate restriction lower the postload insulin response. — *effect:* Lifestyle change can substantially reduce postload insulin secretion as insulin sensitivity improves.; *timeframe:* Weeks to months; *evidence:* observational. Exercise and weight-loss studies reporting OGTT insulin responses.
- **Outcome evidence.** All-cause mortality HR/SD 1.11 (tier C); CVD HR/SD 1.14 (tier C)

### CGM mean glucose

*Also: continuous glucose monitoring mean glucose, average sensor glucose*

- **Mechanism.** CGM mean glucose is the average interstitial glucose over days to weeks, a high-resolution analogue of glycaemic exposure that maps closely onto HbA1c (the glucose management indicator) but captures day-to-day and postprandial glucose that point measures miss.
- **Measurement.** Subcutaneous interstitial glucose sensor (factory-calibrated) worn 10-14 days, sampling every 1-15 min; mean computed over the wear period. — *sample:* wearable; *approx. cost:* $40-90.
- **Signal quality.** moderate — Dense sampling gives a stable mean, but interstitial sensors have ~8-10% MARD vs blood, lag blood glucose, and require sufficient wear time (>=70%) for a representative average.
- **Major confounders.** sensor accuracy/calibration drift (MARD ~8-10%); interstitial-to-blood lag; insufficient wear time; compression artefacts and acetaminophen interference (some sensors); short monitoring window may not represent habitual glycaemia
- **Testing cadence.** intermittent 2-week wears, e.g. quarterly to annually
- **Modifiability (high).** Weight loss, low-glycaemic/lower-carbohydrate diet, post-meal activity; glucose-lowering drugs in diabetes. — *effect:* Diet and activity changes shift CGM mean glucose by ~0.5-1.5 mmol/L within weeks; readily trackable on-device.; *timeframe:* Days to weeks; *evidence:* RCT. CGM-guided lifestyle and dietary intervention trials.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD (documented, no pooled HR)

### CGM glucose variability (CV/SD)

*Also: glycaemic variability, coefficient of variation of glucose, CGM SD*

- **Mechanism.** Glycaemic variability (CV or SD of CGM glucose) captures the amplitude of glucose swings independent of the mean. Large excursions are hypothesized to drive oxidative stress and endothelial dysfunction more than stable hyperglycaemia, and high CV marks hypoglycaemia risk.
- **Measurement.** Computed from a 10-14 day CGM trace as coefficient of variation (SD/mean x 100) or SD of sensor glucose. — *sample:* wearable; *approx. cost:* $40-90.
- **Signal quality.** moderate — Sensitive to sensor noise and to the wear window; CV is more standardized than SD but both depend on adequate wear time and are affected by sensor MARD and compression artefacts.
- **Major confounders.** sensor noise inflates apparent variability; short or low-wear monitoring window; meal timing and composition during wear; physical activity and alcohol; limited population reference ranges in non-diabetic adults
- **Testing cadence.** intermittent 2-week wears
- **Modifiability (moderate).** Smoothing of meal carbohydrate load, post-meal activity, meal sequencing; in diabetes, agents that reduce excursions (GLP-1 agonists, DPP-4 inhibitors). — *effect:* Dietary restructuring can reduce CGM CV by several percentage points within weeks.; *timeframe:* Days to weeks; *evidence:* RCT. CGM-guided dietary intervention trials reporting variability endpoints.
- **Outcome evidence.** All-cause mortality HR/SD 1.20 (tier C); CVD HR/SD 1.25 (tier C)

### CGM time-in-range

*Also: TIR, time in range 70-180 mg/dL*

- **Mechanism.** Time-in-range is the percentage of CGM readings within a target glucose band (commonly 70-180 mg/dL). It is an integrated measure of overall glycaemic control; more time in range reflects fewer hyper- and hypoglycaemic excursions and correlates inversely with diabetic complications.
- **Measurement.** Percentage of CGM glucose values within the target band over a 10-14 day wear period. — *sample:* wearable; *approx. cost:* $40-90.
- **Signal quality.** moderate — Stable when wear time is adequate; depends on sensor accuracy and on the chosen band, and in non-diabetic adults most readings are already in range, compressing its dynamic range.
- **Major confounders.** sensor MARD and calibration; insufficient wear time; ceiling effect in metabolically healthy people; target-band definition differences; meal and activity patterns during wear
- **Testing cadence.** intermittent 2-week wears
- **Modifiability (high).** Lower-glycaemic diet, post-meal activity, meal sequencing; glucose-lowering therapy in diabetes (notably CGM-guided regimens). — *effect:* CGM-guided lifestyle or therapy adjustments can raise time-in-range by 5-15 percentage points within weeks.; *timeframe:* Days to weeks; *evidence:* RCT. CGM-guided management RCTs reporting time-in-range endpoints.
- **Outcome evidence.** All-cause mortality HR/SD 0.88 (tier C); CVD HR/SD 0.85 (tier C)

### CGM postprandial glucose AUC

*Also: postprandial glucose area under the curve, meal glucose response, postprandial glucose excursion*

- **Mechanism.** Postprandial glucose AUC quantifies the size of glucose excursions after meals, integrating peak height and duration. Repeated large postprandial spikes are hypothesized to drive oxidative stress, endothelial dysfunction and beta-cell strain, and postprandial glucose tracks with vascular risk.
- **Measurement.** Incremental or total area under the CGM glucose curve over a defined post-meal window (typically 2-3 h), averaged across meals. — *sample:* wearable; *approx. cost:* $40-90.
- **Signal quality.** noisy — Highly dependent on meal content, timing, prior activity and the AUC definition (incremental vs total, window length); not standardized, with large within-person variability across meals.
- **Major confounders.** meal carbohydrate load and composition; AUC definition and window length differences; prior meal and physical activity; sensor lag and noise; no standardized reference ranges
- **Testing cadence.** intermittent, meal-challenge or 2-week wears
- **Modifiability (high).** Lower-glycaemic meals, reduced meal carbohydrate, meal sequencing (protein/vegetables first), post-meal walking; acarbose/GLP-1 agents in diabetes. — *effect:* Meal restructuring and post-meal activity can roughly halve postprandial glucose excursions.; *timeframe:* Immediate (per meal) to days; *evidence:* RCT. Meal-sequencing and post-meal-activity crossover trials with CGM endpoints.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD (documented, no pooled HR)

### CGM dawn-phenomenon magnitude

*Also: dawn phenomenon, early-morning glucose rise*

- **Mechanism.** The dawn phenomenon is the early-morning rise in glucose driven by the overnight surge in cortisol, growth hormone and catecholamines increasing hepatic glucose output. Its magnitude on CGM is a proposed marker of hepatic insulin resistance and counter-regulatory tone.
- **Measurement.** Difference between pre-dawn nadir and morning glucose on the overnight CGM trace, averaged over several nights. — *sample:* wearable; *approx. cost:* $40-90.
- **Signal quality.** noisy — No standardized definition (nadir-to-peak window varies), large night-to-night variability, and sensitivity to sleep timing, late meals and sensor noise make the metric exploratory.
- **Major confounders.** no standardized computation; sleep timing and quality; late-evening meals or snacks; night-to-night variability; sensor noise overnight
- **Testing cadence.** research/exploratory only
- **Modifiability (moderate).** Earlier and lower-carbohydrate evening meals, improved sleep, weight loss; in diabetes, basal insulin timing or metformin. — *effect:* Evening-meal timing and basal therapy adjustments can blunt the dawn rise, though evidence is limited.; *timeframe:* Days to weeks; *evidence:* observational. Small CGM studies of the dawn phenomenon and evening-meal/sleep manipulation.
- **Outcome evidence.** none recorded

### Fructosamine

*Also: glycated serum protein, glycated albumin (related)*

- **Mechanism.** Fructosamine measures glycated serum proteins (mainly albumin) and reflects average glycaemia over the preceding 2-3 weeks, a shorter window than HbA1c. It is an integrated glycaemic-exposure marker useful when HbA1c is unreliable (haemoglobinopathies, recent change in control).
- **Measurement.** Colorimetric (nitroblue tetrazolium) assay on serum; no fasting required. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Reproducible assay but influenced by serum protein turnover and concentration; less standardized and less outcome-validated than HbA1c, with a shorter and noisier integration window.
- **Major confounders.** low serum albumin or altered protein turnover (nephrotic syndrome, liver disease); thyroid dysfunction; acute inflammation; shorter integration window than HbA1c; no fasting requirement but assay-method differences
- **Testing cadence.** as indicated when HbA1c is unreliable
- **Modifiability (high).** Same glycaemia-lowering levers as HbA1c: weight loss, low-glycaemic diet, activity, glucose-lowering drugs. — *effect:* Responds to glycaemic change faster than HbA1c (2-3 week window); magnitude tracks glucose control.; *timeframe:* 2-3 weeks; *evidence:* observational. Studies of fructosamine response to glucose-lowering therapy.
- **Outcome evidence.** All-cause mortality HR/SD 1.09 (tier B); CVD HR/SD 1.17 (tier B); Dementia HR/SD 1.10 (tier C)

### C-peptide

*Also: connecting peptide, fasting C-peptide*

- **Mechanism.** C-peptide is co-secreted with insulin in equimolar amounts and, being cleared more slowly and not extracted by the liver, gives a stable index of endogenous insulin secretion. Elevated C-peptide indicates compensatory hyperinsulinaemia and insulin resistance; very low values indicate beta-cell failure.
- **Measurement.** Immunoassay on fasting serum or plasma; more stable than insulin and not affected by exogenous insulin therapy. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** moderate — More stable and standardized than insulin, but right-skewed, affected by renal clearance, and its outcome interpretation is non-monotonic (both extremes are adverse).
- **Major confounders.** renal impairment (raises C-peptide); fasting state; right-skewed distribution; U-shaped interpretation (low = beta-cell failure); obesity strongly raises levels
- **Testing cadence.** baseline + as indicated for metabolic assessment
- **Modifiability (moderate).** Weight loss, exercise, carbohydrate restriction reduce hyperinsulinaemia and C-peptide where beta-cell function is intact. — *effect:* Lifestyle change can lower elevated C-peptide substantially in insulin-resistant individuals.; *timeframe:* Weeks to months; *evidence:* observational. Weight-loss and exercise studies reporting C-peptide change.
- **Outcome evidence.** All-cause mortality HR/SD 1.21 (tier B); CVD HR/SD 1.23 (tier B); Cancer HR/SD 1.13 (tier C)

### Adiponectin

*Also: ADIPOQ, GBP-28, AdipoQ*

- **Mechanism.** Adiponectin is an adipocyte-derived hormone that enhances insulin sensitivity, promotes fatty-acid oxidation and has anti-inflammatory and anti-atherogenic actions. Higher levels mark better metabolic health, yet in older and chronically-ill populations high adiponectin paradoxically predicts higher mortality.
- **Measurement.** Immunoassay (ELISA/RIA) on serum or plasma; total or high-molecular-weight adiponectin; assays not fully harmonized. — *sample:* venous blood; *approx. cost:* $30-80.
- **Signal quality.** moderate — Stable analyte with low within-person variability, but assays differ across platforms and the paradoxical, context-dependent mortality association complicates interpretation.
- **Major confounders.** right-skewed distribution; reverse causation: cachexia, heart failure, renal disease raise adiponectin; natriuretic peptides raise adiponectin; sex differences (higher in women); assay/isoform differences
- **Testing cadence.** baseline + as indicated
- **Modifiability (moderate).** Weight loss, exercise, omega-3 and certain drugs (pioglitazone strongly raises adiponectin); smoking cessation. — *effect:* Pioglitazone can double adiponectin; weight loss and exercise raise it 10-30%.; *timeframe:* Weeks to months; *evidence:* RCT meta-analysis. Thiazolidinedione RCTs and weight-loss/exercise meta-analyses of adiponectin.
- **Outcome evidence.** All-cause mortality HR/SD 1.23 (tier A); CVD HR/SD 1.14 (tier B); Dementia HR/SD 1.29 (tier C); Frailty HR/SD 1.13 (tier C)

### Leptin

*Also: LEP, OB protein*

- **Mechanism.** Leptin is an adipocyte hormone signalling energy stores to the hypothalamus to regulate appetite and energy expenditure. In obesity leptin is high but the brain becomes leptin-resistant; elevated leptin marks adiposity and is linked to inflammation, insulin resistance and vascular risk.
- **Measurement.** Immunoassay (ELISA/RIA) on serum or plasma; assays not fully harmonized; strongly sex- and fat-mass dependent. — *sample:* venous blood; *approx. cost:* $30-70.
- **Signal quality.** moderate — Stable analyte but values are dominated by fat mass and sex, strongly right-skewed, and the outcome association is U-shaped and heavily confounded by adiposity.
- **Major confounders.** fat mass dominates the signal; strong sex difference (higher in women); right-skewed distribution; reverse causation: weight loss/cachexia lowers leptin; acute energy balance and sleep
- **Testing cadence.** baseline + as indicated
- **Modifiability (high).** Weight loss is the dominant lever; leptin falls roughly in proportion to fat-mass loss. — *effect:* Substantial weight loss can lower leptin 40-60%; exercise contributes independently.; *timeframe:* Weeks to months; *evidence:* RCT meta-analysis. Weight-loss RCT meta-analyses reporting leptin change.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.20 (tier C); Cancer HR/SD 1.08 (tier C); Dementia HR/SD 0.75 (tier C); Frailty HR/SD 1.11 (tier C)

### Leptin:adiponectin ratio

*Also: LAR, L/A ratio*

- **Mechanism.** The leptin:adiponectin ratio combines a marker of adiposity/leptin resistance (leptin) with an inverse marker of insulin sensitivity (adiponectin), giving an integrated index of adipose-tissue dysfunction. A high ratio reflects insulin resistance, atherogenic risk and inflammatory adipose tissue.
- **Measurement.** Calculated from separately assayed leptin and adiponectin; inherits both assays' lack of harmonization. — *sample:* venous blood; *approx. cost:* $60-140.
- **Signal quality.** moderate — Combines two skewed, sex-dimorphic analytes; the ratio is a better insulin-resistance correlate than either alone but lacks standardized reference ranges and outcome validation.
- **Major confounders.** strong sex difference; right-skewed distribution; fat mass dominates leptin component; reverse causation via the adiponectin paradox; no harmonized assays or cutpoints
- **Testing cadence.** baseline + as indicated
- **Modifiability (high).** Weight loss and exercise lower leptin and raise adiponectin, moving the ratio favourably; pioglitazone raises adiponectin. — *effect:* Weight loss can lower the ratio substantially as leptin falls and adiponectin rises.; *timeframe:* Weeks to months; *evidence:* observational. Weight-loss and exercise studies reporting leptin:adiponectin ratio change.
- **Outcome evidence.** All-cause mortality HR/SD 1.13 (tier C); CVD HR/SD 1.19 (tier C)


## Lipids / Cardiovascular Risk

### ApoB

*Also: apolipoprotein B, apolipoprotein B-100, total apoB*

- **Mechanism.** ApoB is the structural protein of all atherogenic lipoproteins (LDL, VLDL, IDL, Lp(a)); each particle carries exactly one apoB molecule, so apoB concentration equals the number of atherogenic particles. Particle number drives arterial wall infiltration and atherosclerosis, making apoB a more direct causal proxy for atherogenic burden than cholesterol mass.
- **Measurement.** Immunoturbidimetric or immunonephelometric assay on a standard chemistry analyzer; non-fasting acceptable; standardized internationally. — *sample:* venous blood; *approx. cost:* $15-40.
- **Signal quality.** clean — Low biological and analytic variability; non-fasting measurement valid; well-standardized assays. More reproducible than calculated LDL-C.
- **Major confounders.** acute illness/inflammation; pregnancy; recent very high-fat meal (modest); lipid-lowering medication
- **Testing cadence.** baseline + annual; more frequent during lipid-lowering titration
- **Modifiability (high).** Statins, ezetimibe, and PCSK9 inhibitors lower apoB; PCSK9 inhibitors and high-intensity statins produce the largest reductions. — *effect:* High-intensity statin lowers apoB ~35-45%; adding ezetimibe +15-20%; PCSK9 inhibitors lower apoB 45-55%.; *timeframe:* 4-6 weeks to steady state; *evidence:* RCT meta-analysis. CTT Collaboration statin meta-analyses; FOURIER (Sabatine 2017) and ODYSSEY OUTCOMES (Schwartz 2018) PCSK9 RCTs.
- **Outcome evidence.** All-cause mortality HR/SD 1.11 (tier B); CVD HR/SD 1.43 (tier A); Cancer HR/SD 1.10 (tier C); Dementia HR/SD 1.05 (tier B)

### LDL-C

*Also: LDL cholesterol, low-density lipoprotein cholesterol, calculated LDL*

- **Mechanism.** LDL-C measures the cholesterol mass carried by LDL particles, the dominant atherogenic lipoprotein. Cumulative lifetime LDL-C exposure causally drives atherosclerotic plaque formation; genetic, epidemiologic and trial evidence is concordant.
- **Measurement.** Calculated (Friedewald/Martin-Hopkins equation from total cholesterol, HDL-C, triglycerides) or direct enzymatic assay; fasting traditionally preferred for calculated values. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Calculated LDL-C is biased at high triglycerides or low LDL-C; cholesterol content per particle varies, causing discordance with particle number. Direct assays improve reliability.
- **Major confounders.** non-fasting state (affects calculated value); high triglycerides; acute illness; pregnancy; lipid-lowering drugs
- **Testing cadence.** baseline + annual; more frequent during treatment titration
- **Modifiability (high).** Statins (first-line), ezetimibe, PCSK9 inhibitors, bempedoic acid; diet (reduced saturated fat, fiber, plant sterols). — *effect:* High-intensity statin lowers LDL-C ~50%; ezetimibe ~15-20% additional; PCSK9 inhibitor ~50-60%; diet alone ~5-15%.; *timeframe:* 4-6 weeks to steady state; *evidence:* RCT meta-analysis. Cholesterol Treatment Trialists' (CTT) Collaboration meta-analyses, Lancet.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.38 (tier A); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR)

### LDL particle number (LDL-P)

*Also: LDL-P, LDL particle concentration, NMR LDL particle number*

- **Mechanism.** LDL-P counts the actual number of LDL particles, the unit that physically infiltrates the arterial wall. When LDL-P and LDL-C are discordant (e.g., insulin resistance, high triglycerides), LDL-P better tracks atherosclerotic risk because particle number, not cholesterol mass, drives plaque.
- **Measurement.** Nuclear magnetic resonance (NMR) lipoprotein spectroscopy; specialty lab assay. — *sample:* venous blood; *approx. cost:* $50-120.
- **Signal quality.** moderate — NMR is reproducible but less standardized across platforms than apoB; biologically near-equivalent to apoB (both index particle number).
- **Major confounders.** insulin resistance / metabolic syndrome; high triglycerides; acute illness; lipid-lowering drugs; assay platform differences
- **Testing cadence.** baseline; repeat if discordance with LDL-C suspected
- **Modifiability (high).** Same lipid-lowering therapies as LDL-C (statins, ezetimibe, PCSK9 inhibitors); addressing insulin resistance reduces discordance. — *effect:* Statins lower LDL-P broadly in parallel with LDL-C/apoB (~30-50%); residual discordant LDL-P may persist with insulin resistance.; *timeframe:* 4-6 weeks to steady state; *evidence:* RCT. Statin trials with NMR substudies; CTT-class evidence extrapolated via apoB equivalence.
- **Outcome evidence.** CVD HR/SD 1.32 (tier B)

### LDL particle size

*Also: LDL size, small dense LDL, sdLDL, LDL phenotype B*

- **Mechanism.** Small, dense LDL particles penetrate the arterial wall more readily, are retained longer, and are more susceptible to oxidation than large buoyant LDL. A predominance of small dense LDL (pattern B) marks an atherogenic dyslipidemia typical of insulin resistance.
- **Measurement.** NMR spectroscopy, gradient gel electrophoresis, or ion mobility; small-dense LDL-cholesterol can also be measured by a direct homogeneous assay. — *sample:* venous blood; *approx. cost:* $50-120.
- **Signal quality.** noisy — Size is strongly correlated with triglycerides and particle number; once particle number (apoB/LDL-P) is accounted for, independent predictive value of size is weak and inconsistent across studies.
- **Major confounders.** triglyceride level; insulin resistance; particle number (collinearity); assay method differences
- **Testing cadence.** baseline; limited incremental value over apoB/LDL-P
- **Modifiability (moderate).** Triglyceride-lowering measures: weight loss, carbohydrate reduction, exercise, glycemic control; fibrates and high-dose omega-3 shift LDL toward larger particles. — *effect:* Lifestyle and triglyceride reduction can shift LDL phenotype from pattern B to A; magnitude depends on baseline metabolic state.; *timeframe:* 3-6 months; *evidence:* observational. Intervention studies of low-carbohydrate diets and fibrates on LDL subfractions.
- **Outcome evidence.** CVD HR/SD 1.18 (tier B)

### HDL-C

*Also: HDL cholesterol, high-density lipoprotein cholesterol, good cholesterol*

- **Mechanism.** HDL-C measures cholesterol carried by HDL particles, historically viewed as protective via reverse cholesterol transport. However, Mendelian randomization shows HDL-C is not causally protective; it is largely a marker of metabolic health, and both low and very high HDL-C associate with excess mortality.
- **Measurement.** Direct homogeneous enzymatic assay on standard chemistry analyzer; non-fasting acceptable. — *sample:* venous blood; *approx. cost:* $10-25.
- **Signal quality.** moderate — Analytically reliable, but biologically a poor causal target; HDL-C reflects metabolic, inflammatory and lifestyle status rather than a modifiable causal pathway.
- **Major confounders.** alcohol intake; smoking; inflammation/acute illness; metabolic syndrome; estrogen status; genetic variants (CETP)
- **Testing cadence.** baseline + periodic with standard lipid panel
- **Standardization (U-shaped).** HDL-C has a U-shaped association with mortality and dementia — both low and very high levels raise risk — so no single per-SD slope represents it; the off-scale cells use an extreme-vs-reference contrast.
- **Modifiability (low).** Exercise, smoking cessation, weight loss, moderate alcohol modestly raise HDL-C; however raising HDL-C has not been shown to reduce events. — *effect:* Lifestyle changes raise HDL-C ~5-10%; CETP inhibitors raise it >50% but without clinical benefit.; *timeframe:* 3-6 months for lifestyle change; *evidence:* RCT. AIM-HIGH, HPS2-THRIVE (niacin); dal-OUTCOMES, ACCELERATE (CETP inhibitors) - all neutral for events.
- **Outcome evidence.** All-cause mortality native HR 1.68 [categorical, tier B]; CVD HR/SD 0.78 (tier A); Dementia native HR 1.66 [categorical, tier B]

### HDL functionality (cholesterol efflux capacity)

*Also: cholesterol efflux capacity, CEC, HDL function*

- **Mechanism.** Cholesterol efflux capacity is a functional assay measuring HDL's ability to accept cholesterol from macrophages, the first step of reverse cholesterol transport. It captures HDL function rather than HDL-C mass, and is inversely associated with cardiovascular events independent of HDL-C level.
- **Measurement.** Ex vivo assay: serum/apoB-depleted plasma incubated with cholesterol-labeled macrophages; efflux quantified by radiolabel or fluorescence. Research-grade, not standardized for routine clinical use. — *sample:* venous blood; *approx. cost:* $research assay; not routinely priced (estimate 150-400 if offered).
- **Signal quality.** noisy — No standardized clinical assay; results vary by cell type, acceptor preparation and lab protocol; primarily a research tool.
- **Major confounders.** assay methodology differences; inflammation; diabetes; HDL particle concentration; lack of assay standardization
- **Testing cadence.** research use; not part of routine clinical monitoring
- **Modifiability (low).** Exercise, weight loss, and statins modestly improve efflux capacity; no proven therapy reliably raises it with demonstrated clinical benefit. — *effect:* Lifestyle and statin effects on CEC are small and variable.; *timeframe:* months; *evidence:* observational. Intervention substudies of exercise and statins on HDL function.
- **Outcome evidence.** CVD HR/SD 0.68 (tier B)

### Triglycerides

*Also: TG, serum triglycerides, fasting triglycerides*

- **Mechanism.** Triglycerides index triglyceride-rich lipoproteins (VLDL and remnants). The cholesterol carried in these remnant particles is atherogenic; Mendelian randomization links triglyceride-rich lipoprotein pathways causally to coronary disease, though much of TG's risk is mediated through remnant cholesterol and apoB.
- **Measurement.** Enzymatic colorimetric assay; fasting traditionally preferred, non-fasting increasingly accepted with adjusted thresholds. — *sample:* venous blood; *approx. cost:* $10-25.
- **Signal quality.** noisy — High within-person biological variability (diet, alcohol, fasting state); single measurements are unstable and skewed.
- **Major confounders.** fasting/non-fasting state; recent alcohol; recent meal; uncontrolled diabetes; acute illness; pregnancy; estrogen therapy
- **Testing cadence.** baseline + annual; repeat for confirmation given variability
- **Modifiability (high).** Weight loss, carbohydrate/alcohol reduction, exercise, glycemic control; fibrates, high-dose omega-3 fatty acids (icosapent ethyl); statins lower TG modestly. — *effect:* Lifestyle can lower TG 20-50%; fibrates ~30-50%; icosapent ethyl ~20%.; *timeframe:* 4-12 weeks; *evidence:* RCT. REDUCE-IT (icosapent ethyl, Bhatt 2019); fibrate trial meta-analyses.
- **Outcome evidence.** All-cause mortality HR/SD 1.08 (tier B); CVD HR/SD 1.37 (tier A); Dementia HR/SD 0.96 (tier C)

### Triglyceride:HDL ratio

*Also: TG/HDL ratio, TG:HDL-C ratio, triglyceride to HDL cholesterol ratio*

- **Mechanism.** The TG/HDL-C ratio is a surrogate for insulin resistance and atherogenic dyslipidemia (high remnants, small dense LDL, low HDL). It captures the metabolic syndrome lipid pattern in a single number and predicts cardiovascular events and incident diabetes.
- **Measurement.** Calculated from standard lipid panel (triglycerides divided by HDL-C, both in mg/dL or both in mmol/L - units must match). — *sample:* venous blood; *approx. cost:* $10-25.
- **Signal quality.** noisy — Inherits high variability of triglycerides; ratio is unit-dependent and population-dependent (thresholds differ by ethnicity and mg/dL vs mmol/L).
- **Major confounders.** fasting state; alcohol; recent meal; insulin resistance status; ethnicity-specific thresholds; unit convention (mg/dL vs mmol/L)
- **Testing cadence.** baseline + periodic with lipid panel
- **Modifiability (moderate).** Weight loss, low-carbohydrate diet, exercise, glycemic control improve the ratio mainly by lowering triglycerides and raising HDL-C. — *effect:* Lifestyle interventions can lower the ratio substantially (often 25-50%) by improving insulin sensitivity.; *timeframe:* 8-12 weeks; *evidence:* observational. Low-carbohydrate diet and exercise intervention studies on TG/HDL-C.
- **Outcome evidence.** All-cause mortality HR/SD 1.01 (tier B); CVD HR/SD 1.12 (tier B); Dementia HR/SD 1.12 (tier C)

### Lp(a)

*Also: lipoprotein(a), lipoprotein little a, Lp little a*

- **Mechanism.** Lp(a) is an LDL-like particle with apolipoprotein(a) attached to apoB; it is pro-atherogenic, pro-thrombotic and pro-inflammatory. Levels are ~80-90% genetically determined by the LPA gene. Genetic and epidemiologic evidence establishes Lp(a) as a causal, independent risk factor for ASCVD and aortic stenosis.
- **Measurement.** Immunoassay; ideally reported in nmol/L (molar) with an isoform-insensitive assay; non-fasting acceptable; usually a one-time test given genetic stability. — *sample:* venous blood; *approx. cost:* $20-80.
- **Signal quality.** clean — Highly genetically stable over a lifetime; main analytic concern is isoform-dependent mass assays - molar (nmol/L) assays preferred. Once measured, rarely needs repeating.
- **Major confounders.** assay isoform sensitivity (mass vs molar units); acute inflammation (modest rise); renal disease; estrogen status
- **Testing cadence.** once in a lifetime for most adults (genetically stable); recheck only if new therapy or clinical change
- **Modifiability (low).** Largely non-modifiable by lifestyle or diet. PCSK9 inhibitors lower Lp(a) ~20-25%; niacin ~20-30%; targeted antisense/siRNA agents (pelacarsen, olpasiran) lower it 80-90% but cardiovascular outcome trials are ongoing. Lipoprotein apheresis lowers it acutely in severe cases. — *effect:* Diet/exercise: essentially no effect. PCSK9i: ~20-25%. Investigational RNA therapies: 80-90% (outcomes not yet proven).; *timeframe:* Lifestyle: no meaningful change. Drug therapy: weeks.; *evidence:* RCT. Pelacarsen/olpasiran phase 2 trials (Tsimikas 2020; O'Donoghue 2022); HORIZON outcomes trial ongoing.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.13 (tier A); Dementia HR/SD 0.94 (tier C)

### Non-HDL cholesterol

*Also: non-HDL-C, non-high-density lipoprotein cholesterol*

- **Mechanism.** Non-HDL-C (total cholesterol minus HDL-C) captures the cholesterol carried by ALL atherogenic lipoproteins - LDL, VLDL, IDL, remnants, and Lp(a). It is a more complete measure of atherogenic cholesterol burden than LDL-C, especially when triglycerides are elevated, and requires no fasting or extra assay.
- **Measurement.** Calculated from a standard lipid panel (total cholesterol minus HDL-C); no fasting required; no additional cost beyond the panel. — *sample:* venous blood; *approx. cost:* $0-10 (derived from standard panel).
- **Signal quality.** clean — Derived from two well-standardized, low-variability measurements; valid non-fasting; not biased at high triglycerides as calculated LDL-C is.
- **Major confounders.** acute illness/inflammation; pregnancy; lipid-lowering medication
- **Testing cadence.** baseline + annual with lipid panel
- **Modifiability (high).** Statins, ezetimibe, PCSK9 inhibitors; diet (reduced saturated fat, fiber); triglyceride-lowering lifestyle changes also reduce non-HDL-C. — *effect:* High-intensity statin lowers non-HDL-C ~45-50%; further reductions with add-on agents.; *timeframe:* 4-6 weeks to steady state; *evidence:* RCT meta-analysis. CTT Collaboration statin meta-analyses; non-HDL-C is a guideline secondary treatment target.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.50 (tier A)

### Remnant cholesterol

*Also: remnant-C, remnant lipoprotein cholesterol, TRL-cholesterol*

- **Mechanism.** Remnant cholesterol is the cholesterol carried by triglyceride-rich lipoprotein remnants (VLDL, IDL, chylomicron remnants), calculated as total cholesterol minus LDL-C minus HDL-C. Remnant particles penetrate and are retained in the arterial wall and also drive inflammation; Mendelian randomization implicates remnant cholesterol as causal for ischemic heart disease.
- **Measurement.** Calculated from standard lipid panel (total cholesterol minus HDL-C minus LDL-C); can also be directly measured. Non-fasting samples capture remnant burden well. — *sample:* venous blood; *approx. cost:* $0-15 (derived from standard panel).
- **Signal quality.** moderate — Derived value; inherits triglyceride variability and depends on accuracy of the LDL-C estimate, but non-fasting measurement is a strength.
- **Major confounders.** fasting/non-fasting state; recent alcohol/meal; diabetes; method used to estimate LDL-C; acute illness
- **Testing cadence.** baseline + periodic with lipid panel
- **Modifiability (moderate).** Weight loss, reduced carbohydrate and alcohol intake, exercise, glycemic control; statins, fibrates, high-dose omega-3 and (investigationally) ANGPTL3/APOC3 inhibitors lower remnant cholesterol. — *effect:* Lifestyle can lower remnant-C substantially via triglyceride reduction; statins lower it ~20-30%.; *timeframe:* 8-12 weeks; *evidence:* observational. Lifestyle and lipid-lowering intervention studies; MR supports causal benefit of lowering.
- **Outcome evidence.** All-cause mortality HR/SD 1.17 (tier B); CVD HR/SD 1.67 (tier A); Dementia HR/SD 1.04 (tier B)

### Oxidized LDL

*Also: oxLDL, oxidized low-density lipoprotein, ox-LDL*

- **Mechanism.** Oxidized LDL is LDL modified by oxidative stress; it is taken up by macrophage scavenger receptors to form foam cells and triggers endothelial dysfunction and vascular inflammation. Circulating oxLDL is thought to reflect both atherogenic lipoprotein burden and oxidative/inflammatory activity within plaque.
- **Measurement.** ELISA using monoclonal antibodies (e.g., 4E6 or related) against oxidized epitopes on apoB; specialty lab assay. — *sample:* venous blood; *approx. cost:* $40-150.
- **Signal quality.** noisy — Different antibody assays measure different epitopes and are not interchangeable; circulating oxLDL correlates strongly with LDL-C/apoB, so independent signal is modest. Not standardized for clinical use.
- **Major confounders.** overall LDL/apoB level (collinearity); assay/antibody differences; acute inflammation; smoking; diabetes; lack of standardization
- **Testing cadence.** research use; not part of routine clinical monitoring
- **Modifiability (moderate).** Statins lower circulating oxLDL (largely by lowering LDL particle number); smoking cessation and improved glycemic control reduce oxidative modification. — *effect:* Statins reduce oxLDL roughly in proportion to LDL-C reduction.; *timeframe:* weeks to months; *evidence:* observational. Statin intervention studies measuring oxLDL as a secondary endpoint.
- **Outcome evidence.** CVD HR/SD 1.44 (tier B)

### Coronary artery calcium score (CAC)

*Also: CAC, coronary calcium score, Agatston score, calcium scoring CT*

- **Mechanism.** CAC quantifies calcified atherosclerotic plaque in the coronary arteries via non-contrast CT. It is a direct measure of accumulated subclinical atherosclerotic disease burden, integrating lifetime exposure to all risk factors; it is among the strongest predictors of future cardiac events and reclassifies risk beyond standard risk scores.
- **Measurement.** Non-contrast gated cardiac CT; calcified plaque quantified by the Agatston method; low radiation dose (~1 mSv). — *sample:* imaging; *approx. cost:* $100-400.
- **Signal quality.** clean — Highly reproducible imaging measure with low scan-to-scan variability; directly visualizes disease burden rather than a circulating proxy. Main limitation: detects only calcified plaque, can miss early non-calcified plaque.
- **Major confounders.** age (calcium accumulates with age); chronic kidney disease (vascular calcification); scanner/protocol differences; statin use can increase plaque calcification density
- **Testing cadence.** one-time for risk stratification in intermediate-risk adults; repeat scanning generally not recommended at short intervals
- **Standardization (Non-normal).** Coronary artery calcium has a large point mass at zero and a long skewed tail, and risk rises across graded clinical strata rather than along a smooth slope — it is reported by CAC category, not per SD.
- **Modifiability (fixed).** CAC reflects established calcified plaque and does not regress; the score itself cannot be lowered. Statins paradoxically increase calcified plaque density while stabilizing plaque. The actionable response is aggressive risk-factor and lipid management triggered by an elevated score. — *effect:* CAC score is essentially fixed/non-regressing; clinical value is in guiding intensity of preventive therapy.; *timeframe:* not applicable (score does not decrease); *evidence:* observational. Statin-CAC progression studies; guideline use of CAC to guide therapy intensity.
- **Outcome evidence.** All-cause mortality native HR 2.60 [categorical, tier A]; CVD native HR 6.84 [categorical, tier A]; Cancer native HR 2.00 [categorical, tier C]; Dementia native HR 1.71 [categorical, tier B]

### Carotid intima-media thickness (cIMT)

*Also: cIMT, carotid IMT, CCA-IMT, carotid intima-media thickness*

- **Mechanism.** cIMT measures the combined thickness of the intima and media layers of the carotid artery wall by ultrasound. Increased cIMT reflects subclinical atherosclerosis and arterial aging and is associated with future stroke and myocardial infarction, though it adds limited incremental predictive value over standard risk factors.
- **Measurement.** B-mode carotid ultrasound; standardized measurement of mean common carotid artery wall thickness; operator-dependent. — *sample:* imaging; *approx. cost:* $100-300.
- **Signal quality.** moderate — Operator- and protocol-dependent; measurement definitions vary (near vs far wall, segments included), limiting comparability. Reproducible within standardized protocols.
- **Major confounders.** age; operator/sonographer variability; measurement protocol differences (segment, wall); ultrasound equipment differences
- **Testing cadence.** one-time for risk assessment; serial measurement of progression has not proven clinically useful
- **Modifiability (moderate).** Statins, blood pressure control, and lifestyle changes slow or modestly reverse cIMT progression. — *effect:* Statins can slow progression and produce small regression of cIMT; absolute changes are small.; *timeframe:* 1-2 years for measurable change; *evidence:* RCT meta-analysis. Meta-analyses of statin and antihypertensive trials measuring cIMT progression (surrogate endpoint).
- **Outcome evidence.** All-cause mortality HR/SD 1.22 (tier B); CVD HR/SD 1.27 (tier A); Dementia HR/SD 1.08 (tier B)

### ApoE genotype

*Also: APOE genotype, apolipoprotein E genotype, APOE e4, ApoE4 carrier status*

- **Mechanism.** APOE has three common alleles (e2, e3, e4) encoding apolipoprotein E, which mediates clearance of triglyceride-rich and remnant lipoproteins. The e4 allele raises LDL-C and remnant levels and modestly increases coronary heart disease risk; it is also the strongest common genetic risk factor for late-onset Alzheimer's disease. e2 lowers LDL-C but can predispose to type III hyperlipoproteinemia.
- **Measurement.** Genotyping of the two APOE SNPs (rs429358, rs7412) by PCR or array; a single test, results are lifelong. — *sample:* venous blood; *approx. cost:* $50-200.
- **Signal quality.** clean — Genotype is fixed and measured with essentially no biological or analytic variability; the test is definitive once performed.
- **Major confounders.** none for the genotype itself (fixed); phenotypic expression modified by environment, sex, ancestry
- **Testing cadence.** once in a lifetime (genotype is fixed)
- **Standardization (Genotype).** ApoE is a genotype (e2/e3/e4 allele status), not a continuous measurement — there is no standard deviation to standardize against; HRs are e4-carrier group contrasts.
- **Modifiability (fixed).** Genotype is fixed and cannot be changed. e4 carriers benefit from earlier and more intensive management of modifiable risk factors (LDL-C lowering, blood pressure, exercise, sleep) for both cardiovascular and dementia risk. — *effect:* Not applicable - genotype is immutable; downstream lipid and cognitive risk is partially modifiable.; *timeframe:* not applicable; *evidence:* mechanistic. APOE biology; risk-factor modification trials in e4 carriers.
- **Outcome evidence.** All-cause mortality native HR 1.16 [categorical, tier A]; CVD native HR 1.09 [categorical, tier B]; Dementia native HR 3.20 [categorical, tier A]


## Mitochondrial / Cardiorespiratory Fitness

### VO2max (relative)

*Also: maximal oxygen uptake, cardiorespiratory fitness, peak VO2, CRF, aerobic capacity*

- **Mechanism.** VO2max integrates the capacity of the lungs, heart, vasculature and skeletal-muscle mitochondria to deliver and utilize oxygen. It is the single most robust functional predictor of all-cause and cardiovascular mortality, and declining VO2max tracks the loss of physiologic reserve that defines biological aging.
- **Measurement.** Gold standard: graded maximal exercise test (treadmill or cycle ergometer) with expired-gas analysis (cardiopulmonary exercise testing, CPET). Submaximal estimates and treadmill-time-derived METs are common surrogates in large cohorts. — *sample:* functional test; *approx. cost:* $150-400.
- **Signal quality.** clean — Directly measured CPET VO2max is highly reproducible (test-retest CV ~3-5%). Estimated/submaximal values introduce moderate noise; effort dependence and ceiling at true max are the main measurement issues.
- **Major confounders.** age; sex; body composition / adiposity; submaximal vs maximal effort; altitude; beta-blocker use; test modality (treadmill vs cycle); reverse causation from subclinical disease
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (high).** Structured aerobic / endurance training, especially moderate-to-vigorous continuous exercise and high-intensity interval training. — *effect:* Typically 10-25% increase in VO2max with 3-6 months of training; larger relative gains in deconditioned and older individuals, attenuated in already-fit and very old.; *timeframe:* 8-12 weeks for measurable gains; 6 months for substantial improvement.; *evidence:* RCT meta-analysis. Numerous RCT meta-analyses of aerobic and HIIT training (e.g., Bacon 2013 PLoS One; Milanovic 2015 Sports Med) consistently show ~10-25% VO2max improvement.
- **Outcome evidence.** All-cause mortality HR/SD 0.76 (tier A); CVD HR/SD 0.72 (tier A); Cancer HR/SD 0.64 (tier B); Dementia HR/SD 0.81 (tier B); Frailty (documented, no pooled HR)

### Ventilatory thresholds (VT1, VT2)

*Also: ventilatory threshold, anaerobic threshold, gas exchange threshold, respiratory compensation point, VT1, VT2*

- **Mechanism.** Ventilatory thresholds mark the exercise intensities at which ventilation rises disproportionately as metabolic acidosis develops (VT1 ~aerobic/gas-exchange threshold, VT2 ~respiratory compensation point). They index submaximal mitochondrial oxidative capacity and the workload sustainable without progressive acidosis, and they decline with age and deconditioning.
- **Measurement.** Identified during cardiopulmonary exercise testing from ventilatory equivalents and gas-exchange slopes (V-slope method). — *sample:* functional test; *approx. cost:* $150-400.
- **Signal quality.** moderate — Threshold determination is somewhat operator- and method-dependent with inter-rater variability; physiologically meaningful but harder to standardize than VO2max.
- **Major confounders.** age; sex; training status; test protocol / ramp rate; rater subjectivity; medications affecting ventilation
- **Testing cadence.** baseline + annual (with CPET)
- **Modifiability (high).** Aerobic endurance training, particularly threshold and tempo work. — *effect:* Ventilatory thresholds shift to higher absolute and relative workloads with training, often improving more than VO2max in percentage terms.; *timeframe:* 6-12 weeks.; *evidence:* RCT. Endurance training trials consistently show upward shifts in ventilatory/lactate thresholds (general exercise physiology literature).
- **Outcome evidence.** CVD (documented, no pooled HR)

### Resting metabolic rate

*Also: RMR, basal metabolic rate, BMR, resting energy expenditure*

- **Mechanism.** RMR is the energy expended at rest, dominated by organ and lean-tissue mitochondrial activity. Its relationship to aging is complex: low RMR is linked to longevity in cross-species comparisons, yet in humans higher RMR predicts mortality, an association largely confounded by fat mass and subclinical disease.
- **Measurement.** Indirect calorimetry (ventilated hood) under fasted, rested conditions; predictive equations are common but less accurate. — *sample:* in-office device; *approx. cost:* $75-250.
- **Signal quality.** noisy — Strongly determined by lean and fat mass; once body composition is accounted for, residual RMR signal is small and easily confounded by recent food, sleep, temperature, thyroid status and illness.
- **Major confounders.** fat mass; lean mass; thyroid function; subclinical illness / inflammation; ambient temperature; recent food intake; menstrual phase; medications
- **Testing cadence.** not routinely tracked; baseline if metabolic assessment indicated
- **Modifiability (low).** Resistance training to increase lean mass; RMR otherwise tightly regulated and not a useful intervention target. — *effect:* Small (a few percent) changes via lean-mass gain; adaptive thermogenesis lowers RMR with weight loss.; *timeframe:* Months.; *evidence:* observational. Body-composition and weight-change literature on RMR adaptation.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR)

### Substrate utilization (RER / maximal fat oxidation)

*Also: respiratory exchange ratio, RER, RQ, maximal fat oxidation, MFO, FATmax, metabolic flexibility*

- **Mechanism.** The respiratory exchange ratio and maximal fat oxidation reflect the balance between fat and carbohydrate combustion and the ability to switch fuels (metabolic flexibility). Impaired fat oxidation and blunted metabolic flexibility are features of insulin resistance and mitochondrial dysfunction, mechanistically tied to cardiometabolic aging.
- **Measurement.** Indirect calorimetry at rest (RER) and during graded exercise to derive maximal fat oxidation rate and FATmax intensity. — *sample:* in-office device; *approx. cost:* $100-350.
- **Signal quality.** noisy — FATmax/MFO have well-documented reproducibility problems (test protocol, prior diet, glycogen status, day-to-day variability); RER is highly sensitive to recent carbohydrate intake.
- **Major confounders.** recent diet / glycogen status; fasting duration; exercise protocol; training status; fitness level; caffeine; time of day
- **Testing cadence.** not routinely tracked
- **Modifiability (moderate).** Endurance and FATmax-intensity exercise training; weight loss; low-carbohydrate dietary periodization. — *effect:* Maximal fat oxidation and metabolic flexibility improve meaningfully with training; FATmax-intensity training programs reduce body fat, visceral fat and HOMA-IR.; *timeframe:* 6-12 weeks.; *evidence:* RCT meta-analysis. Systematic review/meta-analysis of FATmax-intensity exercise training for cardiometabolic health in overweight/obesity (Sports Med / Med Sci Sports Exerc literature, 2024-2025).
- **Outcome evidence.** none recorded

### Exercise economy

*Also: running economy, movement economy, oxygen cost of locomotion, metabolic cost of walking*

- **Mechanism.** Exercise economy is the oxygen (energy) cost of locomotion at a fixed submaximal speed. Poorer walking economy means greater relative effort for daily tasks, contributes to fatigue and activity avoidance, and rises with age - a determinant of mobility decline and loss of independence.
- **Measurement.** Indirect calorimetry during steady-state submaximal treadmill walking or running at standardized speeds. — *sample:* functional test; *approx. cost:* $100-300.
- **Signal quality.** moderate — Reasonably reproducible within a session, but absolute values depend on test speed, footwear, gait mechanics and familiarization; cross-study comparability is limited.
- **Major confounders.** test speed; body mass; gait mechanics / orthopedic status; footwear; training status; age
- **Testing cadence.** not routinely tracked
- **Modifiability (moderate).** Endurance training, strength and gait training, weight loss. — *effect:* Walking/running economy improves modestly with training and weight reduction; gains are smaller and slower than VO2max changes.; *timeframe:* Weeks to months.; *evidence:* observational. Exercise physiology literature on economy adaptation; aging mobility studies.
- **Outcome evidence.** Frailty (documented, no pooled HR)

### Lactate threshold and lactate kinetics

*Also: lactate threshold, LT1, LT2, onset of blood lactate accumulation, OBLA, maximal lactate steady state, lactate clearance*

- **Mechanism.** The lactate threshold is the exercise intensity above which blood lactate accumulates faster than it can be cleared, reflecting the balance between glycolytic flux and mitochondrial oxidative capacity. A higher threshold and faster lactate clearance indicate better skeletal-muscle oxidative function, which declines with age and deconditioning.
- **Measurement.** Serial capillary (fingertip/earlobe) blood lactate sampling during a graded exercise test; thresholds identified from the lactate-workload curve. — *sample:* capillary blood; *approx. cost:* $100-300.
- **Signal quality.** moderate — Lactate values are reproducible but threshold definitions vary widely (fixed 2/4 mmol/L vs individualized methods), and results depend on protocol, prior diet and glycogen status.
- **Major confounders.** test protocol / stage duration; prior diet and glycogen; training status; sampling site; hydration; ambient temperature
- **Testing cadence.** not routinely tracked; used in athletic monitoring
- **Modifiability (high).** Endurance training, especially threshold-intensity and high-volume aerobic training. — *effect:* Lactate threshold shifts substantially to higher workloads with training; often the most trainable submaximal aerobic parameter.; *timeframe:* 6-12 weeks.; *evidence:* RCT. Endurance-training literature consistently shows rightward shift of the lactate curve.
- **Outcome evidence.** none recorded

### Grip strength

*Also: handgrip strength, hand grip strength, muscle strength*

- **Mechanism.** Grip strength is a simple, reliable proxy for total-body muscle strength and overall physiologic reserve. Low grip strength reflects sarcopenia, neuromuscular decline and systemic frailty, and is one of the strongest single-measure predictors of mortality and disability in aging.
- **Measurement.** Maximal isometric force with a calibrated hand dynamometer (e.g., Jamar), best of multiple trials, typically dominant hand. — *sample:* in-office device; *approx. cost:* $20-80.
- **Signal quality.** clean — Highly reproducible with standardized protocol and calibrated dynamometer; main variability sources are posture, instructions and effort, all controllable.
- **Major confounders.** age; sex; body size / height; hand dominance; arthritis or hand injury; dynamometer type and protocol; reverse causation from subclinical disease
- **Testing cadence.** baseline + annual
- **Modifiability (high).** Progressive resistance training, including grip- and forearm-specific loading. — *effect:* Resistance training increases grip and overall strength by ~10-30% in older adults within months; among the most trainable functional markers.; *timeframe:* 8-12 weeks.; *evidence:* RCT meta-analysis. Resistance-training RCT meta-analyses in older adults consistently show meaningful strength gains.
- **Outcome evidence.** All-cause mortality HR/SD 1.27 (tier A); CVD HR/SD 1.29 (tier A); Cancer HR/SD 1.16 (tier B); Dementia HR/SD 1.54 (tier B); Frailty HR/SD 1.29 (tier B)

### Lower-body strength / power relative to body mass

*Also: leg press strength, sit-to-stand power, relative muscle power, lower-extremity strength, leg power*

- **Mechanism.** Lower-body strength and power relative to body mass govern the ability to rise from a chair, climb stairs and recover balance. Power (force x velocity) declines faster than strength with age and is a strong determinant of mobility, fall risk and loss of independence.
- **Measurement.** Leg-press 1-RM or isokinetic dynamometry; field surrogates include the 5- or 30-second sit-to-stand and sit-to-stand power estimated from chair-rise time, height and body mass. — *sample:* functional test; *approx. cost:* $0-150.
- **Signal quality.** moderate — Laboratory dynamometry is reproducible; field power estimates from chair-rise are practical but introduce estimation error and depend on technique and motivation.
- **Major confounders.** body mass; joint pain / orthopedic limitation; test technique; motivation / effort; age; sex
- **Testing cadence.** baseline + annual
- **Standardization (Clinical cutpoint).** The available HRs contrast a clinically-defined low-strength group against normal — a threshold split of unknown SD width, neither a quantile nor a per-unit estimate, so it cannot be standardized.
- **Modifiability (high).** Progressive resistance training, especially high-velocity / power-oriented training and explosive sit-to-stand exercise. — *effect:* Lower-body strength improves 25-100%+ and muscle power improves substantially in older adults with 8-12 weeks of resistance/power training; power-oriented training is more effective than slow heavy training for power gains.; *timeframe:* 8-12 weeks.; *evidence:* RCT meta-analysis. Power- and resistance-training RCT meta-analyses in older adults.
- **Outcome evidence.** All-cause mortality native HR 1.57 [categorical, tier B]; CVD native HR 0.86 [categorical, tier B]; Frailty (documented, no pooled HR)

### Muscle mass (appendicular lean mass index)

*Also: appendicular lean mass, ALM, ALMI, skeletal muscle mass index, lean body mass, sarcopenia*

- **Mechanism.** Appendicular lean mass index (limb lean mass / height^2) quantifies skeletal-muscle quantity, the substrate for strength, metabolic glucose disposal and protein reserve. Low muscle mass (sarcopenia) reflects anabolic decline and predicts disability and mortality, though strength predicts outcomes more strongly than mass alone.
- **Measurement.** Dual-energy X-ray absorptiometry (DEXA) appendicular lean mass; bioelectrical impedance is a lower-cost surrogate. — *sample:* imaging; *approx. cost:* $75-250.
- **Signal quality.** moderate — DEXA lean mass is reproducible (CV ~1-2%), but hydration status affects readings and lean mass includes non-muscle tissue; muscle quality/strength is not captured.
- **Major confounders.** body size / height normalization; hydration status; fat mass (obesity); DEXA vs BIA method; age; sex; ethnicity
- **Testing cadence.** baseline + every 1-2 years
- **Standardization (Cutpoint / rank).** Muscle-mass HRs here use a sarcopenia diagnostic cutoff (CVD) and a percentile-rank exposure (cancer); neither carries the population SD needed to convert to a per-SD basis.
- **Modifiability (moderate).** Progressive resistance training, with adequate dietary protein; combined with anabolic stimulus in deficient states. — *effect:* Resistance training increases appendicular lean mass by roughly 0.5-2 kg over 3-6 months in older adults; gains are modest in magnitude relative to strength gains.; *timeframe:* 12-24 weeks.; *evidence:* RCT meta-analysis. Resistance-training RCT meta-analyses (e.g., Peterson 2011 Ageing Res Rev) show significant lean-mass gains in older adults.
- **Outcome evidence.** All-cause mortality HR/SD 0.50 (tier B); CVD native HR 1.43 [categorical, tier B]; Cancer native HR 0.87 [unconvertible, tier C]; Dementia HR/SD 1.33 (tier B); Frailty (documented, no pooled HR)

### Gait speed

*Also: walking speed, usual gait speed, habitual gait speed, 4-meter walk*

- **Mechanism.** Usual gait speed is an integrative 'vital sign' that summarizes the function of the cardiovascular, musculoskeletal, neurological and energetic systems. It declines with age, predicts survival and disability across populations, and is one of the most validated single functional measures in geriatrics.
- **Measurement.** Time to walk a measured course (typically 4 m or 6 m) at usual pace, often the better of two trials; can also be measured over longer distances or with instrumented walkways. — *sample:* functional test; *approx. cost:* $0-50.
- **Signal quality.** clean — Highly reproducible and standardized; main variability is course length, static vs dynamic start and usual-pace instruction, all controllable.
- **Major confounders.** age; height / leg length; orthopedic and neurological conditions; course length and start protocol; footwear; test environment
- **Testing cadence.** baseline + annual
- **Modifiability (moderate).** Multicomponent exercise (aerobic + resistance + balance/gait training); resistance and power training in deconditioned older adults. — *effect:* Exercise interventions improve gait speed by roughly 0.05-0.10+ m/s in older adults; ~0.05 m/s is a clinically meaningful change.; *timeframe:* 12-24 weeks.; *evidence:* RCT meta-analysis. Exercise-intervention RCT meta-analyses in older adults show significant gait-speed improvement.
- **Outcome evidence.** All-cause mortality HR/SD 0.71 (tier A); CVD HR/SD 1.23 (tier B); Cancer HR/SD 1.00 (tier B); Dementia HR/SD 1.59 (tier B); Frailty (documented, no pooled HR)

### Sit-to-stand (chair-rise test)

*Also: chair-rise test, five-times sit-to-stand, 5xSTS, 30-second chair stand, sitting-rising test*

- **Mechanism.** The chair-rise / sit-to-stand test measures lower-body strength, power and balance during a functionally essential task. Poor performance signals loss of mobility reserve, fall risk and frailty, and predicts mortality independent of conventional risk factors.
- **Measurement.** Time for five chair rises (5xSTS), number of rises in 30 seconds, or the floor-based sitting-rising test scored 0-10. — *sample:* functional test; *approx. cost:* $0-25.
- **Signal quality.** moderate — Reproducible with standardized protocol; depends on chair height, arm use, technique and motivation, and different versions (timed, count, floor-based) are not interchangeable.
- **Major confounders.** chair height / arm use; joint pain and orthopedic limitation; test version; motivation / effort; age; balance impairment
- **Testing cadence.** baseline + annual
- **Modifiability (high).** Progressive resistance and power training targeting the lower body; practiced sit-to-stand exercise. — *effect:* Chair-rise time and repetitions improve substantially with 8-12 weeks of resistance/power training in older adults.; *timeframe:* 8-12 weeks.; *evidence:* RCT meta-analysis. Resistance- and power-training RCT meta-analyses in older adults.
- **Outcome evidence.** All-cause mortality HR/SD 1.87 (tier B); CVD HR/SD 1.81 (tier C); Frailty (documented, no pooled HR)

### VO2 kinetics (mean response time)

*Also: oxygen uptake kinetics, VO2 mean response time, VO2 on-kinetics, VO2 recovery kinetics, tau*

- **Mechanism.** VO2 kinetics describe how rapidly oxygen uptake rises at the onset of submaximal exercise (or recovers afterward). Faster kinetics indicate efficient oxygen delivery and mitochondrial oxidative responsiveness; sluggish kinetics reflect impaired microvascular and mitochondrial function and predict poorer functional mobility in older and chronically ill adults.
- **Measurement.** Breath-by-breath gas analysis during constant-load submaximal exercise transitions or recovery, fitted to exponential models to derive mean response time / tau. — *sample:* functional test; *approx. cost:* $150-400.
- **Signal quality.** noisy — Breath-by-breath data are inherently noisy; reliable kinetic parameter estimation requires multiple transitions and signal averaging, and modeling choices affect results. Not standardized for routine use.
- **Major confounders.** work-rate domain (moderate vs heavy); number of transitions averaged; modeling/curve-fitting method; fitness status; prior exercise (priming effect); age
- **Testing cadence.** not routinely tracked; research/clinical-physiology use
- **Modifiability (moderate).** Endurance / aerobic training; high-intensity interval training. — *effect:* VO2 on-kinetics speed up (shorter tau / mean response time) with aerobic training, particularly in older and deconditioned individuals.; *timeframe:* 6-12 weeks.; *evidence:* RCT. Exercise-training studies show acceleration of VO2 kinetics with endurance training.
- **Outcome evidence.** Frailty (documented, no pooled HR)

### Heart rate recovery

*Also: HRR, post-exercise heart rate recovery, heart rate recovery at 1 minute*

- **Mechanism.** Heart rate recovery is the fall in heart rate in the first 1-2 minutes after peak exercise, driven mainly by parasympathetic (vagal) reactivation. Attenuated recovery indicates autonomic dysfunction, an independent marker of cardiovascular risk and mortality.
- **Measurement.** Difference between peak heart rate and heart rate at 1 (or 2) minutes after a graded exercise test, with a standardized cool-down. — *sample:* functional test; *approx. cost:* $100-300.
- **Signal quality.** moderate — Reproducible when cool-down protocol (active vs passive recovery, posture) is standardized; protocol differences affect absolute values and abnormality cutoffs vary across studies.
- **Major confounders.** cool-down protocol (active vs passive); body position during recovery; beta-blockers and rate-limiting drugs; fitness level; peak effort achieved; age
- **Testing cadence.** baseline + annual (with exercise testing)
- **Modifiability (moderate).** Aerobic / endurance exercise training, which enhances vagal tone. — *effect:* Heart rate recovery improves by several bpm with regular aerobic training; trained individuals show faster recovery.; *timeframe:* 8-12 weeks.; *evidence:* RCT. Exercise-training trials demonstrate improved heart rate recovery with aerobic training.
- **Outcome evidence.** All-cause mortality HR/SD 1.11 (tier A); CVD HR/SD 1.16 (tier B); Dementia HR/SD 0.90 (tier B)

### Heart rate variability (HRV)

*Also: HRV, SDNN, RMSSD, vagal tone, autonomic balance*

- **Mechanism.** Heart rate variability quantifies beat-to-beat variation in heart rate, indexing autonomic (mainly parasympathetic) regulation. Higher HRV reflects healthy vagal tone and physiologic adaptability; HRV declines with age and disease and lower values predict mortality and cardiovascular events.
- **Measurement.** Time-domain (SDNN, RMSSD) and frequency-domain measures from ECG or photoplethysmography over short recordings (~5 min), 24-hour Holter, or consumer wearables. — *sample:* wearable; *approx. cost:* $0-300.
- **Signal quality.** noisy — HRV is highly state-dependent (posture, breathing, stress, sleep, caffeine, time of day) with large day-to-day variability; recording length and device differ across studies. Trends are more informative than single values.
- **Major confounders.** age; posture and breathing rate; physical activity and stress before recording; recording length; device / measurement method; medications (beta-blockers); alcohol and caffeine; time of day / circadian phase
- **Testing cadence.** trend monitoring (e.g., daily wearable) or baseline + annual clinical recording
- **Modifiability (moderate).** Aerobic exercise training; secondarily slow-paced breathing / HRV-biofeedback and improved sleep and stress management. — *effect:* Aerobic training produces modest increases in HRV (RMSSD/SDNN); biofeedback can acutely raise HRV but durable resting changes are smaller.; *timeframe:* 8-16 weeks.; *evidence:* RCT meta-analysis. Meta-analyses of exercise training and of HRV-biofeedback show modest HRV improvements.
- **Outcome evidence.** All-cause mortality HR/SD 1.19 (tier B); CVD HR/SD 1.21 (tier B); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)


## Inflammation

### hsCRP

*Also: high-sensitivity C-reactive protein, CRP, hs-CRP*

- **Mechanism.** Acute-phase protein synthesized by the liver in response to IL-6; serves as a downstream integrator of systemic low-grade inflammation ('inflammaging'). Elevated hsCRP marks the chronic inflammatory state that accompanies atherosclerosis, metabolic dysfunction, and biological aging, though Mendelian randomization indicates it is a marker rather than a cause of disease.
- **Measurement.** High-sensitivity immunoturbidimetric or immunonephelometric assay on a standard clinical chemistry analyzer. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Analytically robust and well-standardized, but high within-person biological variability: levels spike acutely with infection, trauma, or recent exercise. Single measurements are noisy; serial measurements or values >10 mg/L should be excluded as acute reactions. Tertile/quartile stability over years is moderate.
- **Major confounders.** acute infection or recent illness; obesity and adiposity; smoking; estrogen / hormone therapy; recent vigorous exercise or trauma; chronic conditions (RA, periodontal disease); statin and NSAID use
- **Testing cadence.** baseline + annual; repeat if >10 mg/L to exclude acute inflammation
- **Modifiability (high).** Statin therapy, weight loss, and regular aerobic exercise; canakinumab/colchicine for residual inflammatory risk. — *effect:* Rosuvastatin lowered hsCRP ~37% in JUPITER; diet-induced weight loss produces a comparable reduction (roughly proportional to kg lost); regular exercise lowers hsCRP by ~20-30%.; *timeframe:* Weeks for statins; 3-6 months for weight loss and exercise.; *evidence:* RCT. Ridker PM et al. JUPITER trial. N Engl J Med. 2008;359:2195-2207. POUNDS LOST trial (weight loss & hsCRP), Am J Clin Nutr 2012.
- **Outcome evidence.** All-cause mortality HR/SD 1.54 (tier A); CVD HR/SD 1.37 (tier A); Cancer HR/SD 1.11 (tier B); Dementia HR/SD 1.37 (tier B); Frailty HR/SD 1.27 (tier C)

### IL-6

*Also: interleukin-6, IL6*

- **Mechanism.** Pleiotropic pro-inflammatory cytokine and the principal driver of hepatic acute-phase protein synthesis (including CRP and fibrinogen). A central mediator of 'inflammaging'; unlike CRP, human genetic (IL6R) evidence supports a causal role in atherosclerotic disease, making the IL-6 pathway a validated therapeutic target.
- **Measurement.** ELISA or high-sensitivity immunoassay (e.g., electrochemiluminescence); proteomic panels (Olink, SomaScan) also quantify it. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** moderate — Biologically meaningful but analytically more variable than CRP: short half-life, diurnal variation, sensitivity to recent exercise/infection, and assay-platform differences. Single measurements are noisy; cross-assay comparability is limited.
- **Major confounders.** acute infection / illness; obesity (adipose tissue secretes IL-6); recent exercise; smoking; diurnal rhythm; assay platform variability; chronic inflammatory disease
- **Testing cadence.** baseline + annual; not yet routine in standard care
- **Modifiability (high).** Direct pathway blockade (IL-6R antibody tocilizumab; IL-1beta antibody canakinumab; colchicine) lowers IL-6 axis activity; lifestyle (exercise, weight loss) produces smaller reductions. — *effect:* Canakinumab (CANTOS) reduced hsCRP ~37% and IL-6 substantially with a 15% cut in major cardiovascular events; tocilizumab strongly suppresses CRP/fibrinogen. Exercise and weight loss lower IL-6 modestly (~10-25%).; *timeframe:* Days to weeks for biologics; months for lifestyle change.; *evidence:* RCT. Ridker PM et al. CANTOS trial (canakinumab). N Engl J Med. 2017;377:1119-1131. IL6R MR consortium, Lancet 2012.
- **Outcome evidence.** All-cause mortality HR/SD 1.48 (tier A); CVD HR/SD 1.25 (tier A); Cancer HR/SD 1.20 (tier C); Dementia HR/SD 1.32 (tier B); Frailty HR/SD 1.14 (tier C)

### TNF-alpha

*Also: tumor necrosis factor alpha, TNF-a, TNFa*

- **Mechanism.** Pro-inflammatory cytokine central to innate immunity, apoptosis signaling, insulin resistance, and cachexia. Implicated in inflammaging and tissue catabolism; its soluble receptors (sTNFR1/sTNFR2) are more stable, prognostically informative analytes than free TNF-alpha itself.
- **Measurement.** High-sensitivity ELISA or multiplex immunoassay; soluble TNF receptors (sTNFR1/2) often preferred for stability. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** noisy — Free circulating TNF-alpha is present at very low concentrations with a short half-life, marked diurnal variation, and poor assay reproducibility near the detection limit. Soluble receptors are far more reliable; most robust prognostic data come from sTNFR measurements rather than TNF-alpha itself.
- **Major confounders.** assay sensitivity limits and platform differences; diurnal variation; acute infection; obesity / adipose secretion; renal function (affects soluble receptor clearance); smoking
- **Testing cadence.** research / specialist use; not routine
- **Modifiability (moderate).** Anti-TNF biologics (etanercept, infliximab, adalimumab) in inflammatory disease; lifestyle (weight loss, exercise) for modest reductions in the general population. — *effect:* Anti-TNF therapy strongly suppresses TNF signaling and, in observational meta-analysis of rheumatoid arthritis cohorts, was associated with reduced cardiovascular events (RR ~0.46); lifestyle change yields modest reductions.; *timeframe:* Days to weeks for biologics; months for lifestyle.; *evidence:* observational. Barnabe C et al. Systematic review and meta-analysis: anti-TNF therapy and cardiovascular events in rheumatoid arthritis. Arthritis Care Res. 2011;63(4):522-529. PMID: 20957658.
- **Outcome evidence.** All-cause mortality HR/SD 1.34 (tier C); CVD HR/SD 1.63 (tier C); Cancer HR/SD 1.32 (tier C); Dementia HR/SD 1.55 (tier C); Frailty (documented, no pooled HR)

### GlycA

*Also: glycoprotein acetyls, NMR glycoprotein acetylation signal, N-acetyl glycoprotein*

- **Mechanism.** Composite NMR signal arising from N-acetyl methyl groups of glycosylated acute-phase proteins (alpha-1-acid glycoprotein, haptoglobin, alpha-1-antitrypsin, alpha-1-antichymotrypsin, transferrin). Integrates both protein abundance and glycosylation state, capturing a stable, time-averaged measure of chronic systemic inflammation that predicts long-term disease risk.
- **Measurement.** Proton (1H) nuclear magnetic resonance spectroscopy of plasma/serum, typically as part of a multi-analyte NMR metabolomics panel. — *sample:* venous blood; *approx. cost:* $30-100.
- **Signal quality.** clean — Lower within-person biological variability than CRP because it averages several glycoproteins with longer half-lives; NMR quantification is highly reproducible. Less sensitive to transient acute spikes, making it a comparatively stable inflammation index. Newer marker with less standardization across labs.
- **Major confounders.** obesity / adiposity; smoking; acute infection (less than CRP); metabolic syndrome and diabetes; NMR platform / lab differences; statin use
- **Testing cadence.** baseline + every few years; not yet routine clinical use
- **Modifiability (moderate).** Weight loss, exercise, and statin therapy lower GlycA, paralleling effects on other acute-phase markers. — *effect:* Statins and lifestyle interventions reduce GlycA modestly (roughly 5-15%); magnitude tracks reductions in adiposity and CRP.; *timeframe:* Months.; *evidence:* observational. Observational and secondary-analysis data (e.g., statin trials with NMR substudies); no dedicated large RCT with GlycA as primary endpoint.
- **Outcome evidence.** All-cause mortality HR/SD 1.24 (tier B); CVD HR/SD 1.22 (tier B); Cancer HR/SD 1.08 (tier C); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)

### Neutrophil:lymphocyte ratio (NLR)

*Also: NLR, neutrophil-to-lymphocyte ratio*

- **Mechanism.** Ratio derived from a standard complete blood count; integrates the innate (neutrophil) and adaptive (lymphocyte) immune compartments. A high ratio reflects active inflammatory/stress responses with relative lymphopenia, a pattern associated with immunosenescence, physiologic stress, and adverse aging trajectories.
- **Measurement.** Calculated from an automated complete blood count with differential (absolute neutrophil count divided by absolute lymphocyte count). — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** noisy — Cheap and ubiquitous (derived from routine CBC), but highly non-specific and labile: shifts acutely with infection, physical/psychological stress, corticosteroids, and time of day. No standardized cut-points; tertile/quartile definitions vary widely across studies, limiting comparability.
- **Major confounders.** acute infection or illness; corticosteroid and other immunomodulatory drugs; physical / psychological stress; diurnal variation; smoking; absence of standardized thresholds; hematologic conditions
- **Testing cadence.** opportunistic -- available on any CBC; interpret with clinical context
- **Modifiability (low).** No NLR-specific intervention; it falls with resolution of underlying inflammation/illness and improves modestly with exercise and weight loss. — *effect:* Not well quantified; NLR normalizes when an acute or chronic inflammatory driver is treated.; *timeframe:* Days to weeks (acute resolution) to months (lifestyle).; *evidence:* observational. No dedicated RCT targeting NLR as an endpoint; inferred from CBC normalization in observational data.
- **Outcome evidence.** All-cause mortality HR/SD 1.05 (tier B); CVD HR/SD 1.04 (tier B); Cancer HR/SD 1.01 (tier B); Dementia HR/SD 1.16 (tier B); Frailty HR/SD 1.33 (tier C)

### Fibrinogen

*Also: plasma fibrinogen, factor I*

- **Mechanism.** Hepatically synthesized acute-phase glycoprotein and the terminal substrate of the coagulation cascade. Sits at the intersection of inflammation and thrombosis: elevated levels increase blood viscosity, platelet aggregation, and clot formation, mechanistically linking systemic inflammation to atherothrombotic events.
- **Measurement.** Clauss clotting assay (functional) or immunoassay (antigen) on a coagulation analyzer. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Reasonably standardized assay but, as an acute-phase reactant, levels rise with infection, smoking, pregnancy, and inflammation. Within-person variability is moderate; longer half-life than CRP gives somewhat more stable readings.
- **Major confounders.** acute infection / inflammation; smoking; pregnancy and estrogen / hormone therapy; obesity; diabetes; assay method (Clauss vs antigen); age
- **Testing cadence.** baseline + periodic; available within coagulation panels
- **Modifiability (moderate).** Smoking cessation, weight loss, and regular physical activity lower fibrinogen; no fibrinogen-specific drug therapy in primary prevention. — *effect:* Smoking cessation and exercise reduce fibrinogen modestly (roughly 0.2-0.5 g/L over months); fibrate drugs lower it somewhat but without proven outcome benefit attributable to fibrinogen.; *timeframe:* Months.; *evidence:* observational. Observational and secondary trial data; no RCT demonstrates outcome benefit from fibrinogen lowering per se.
- **Outcome evidence.** All-cause mortality HR/SD 1.70 (tier A); CVD HR/SD 1.94 (tier A); Cancer HR/SD 1.12 (tier B); Dementia HR/SD 1.30 (tier B); Frailty HR/SD 1.31 (tier C)

### Homocysteine

*Also: total homocysteine, tHcy, hyperhomocysteinemia*

- **Mechanism.** Sulfur-containing amino acid intermediate of methionine metabolism. Elevated levels reflect impaired one-carbon metabolism (folate/B12/B6 status, renal function, MTHFR genotype) and have been linked to endothelial dysfunction, oxidative stress, and neurotoxicity. Often grouped with inflammatory/vascular risk markers though it is not a classic acute-phase reactant.
- **Measurement.** Immunoassay or LC-MS/MS on fasting plasma; sample must be processed promptly as levels rise with delayed serum separation. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** moderate — Analytically reliable but strongly determined by modifiable nutritional status (folate, B12, B6) and renal function; pre-analytic handling matters (levels rise if red cells are not separated quickly). Genetically influenced (MTHFR). Less labile day-to-day than cytokines.
- **Major confounders.** folate / vitamin B12 / B6 status; renal impairment (major determinant); MTHFR C677T genotype; age and male sex; hypothyroidism; delayed sample processing; certain drugs (methotrexate, antiepileptics)
- **Testing cadence.** baseline; recheck after B-vitamin repletion if elevated
- **Standardization (Cutpoint (frailty)).** Homocysteine standardizes normally for mortality, CVD and dementia; only the frailty cell is off-scale, because that study reported a clinical >=15 micromol/L cutpoint contrast rather than a continuous estimate.
- **Modifiability (high).** Folic acid plus vitamin B12 (and B6) supplementation reliably lowers plasma homocysteine. — *effect:* B-vitamin supplementation lowers homocysteine by roughly 25-30%; however this lowering did NOT translate into reduced cardiovascular events or mortality in large RCTs.; *timeframe:* Weeks for the biochemical change.; *evidence:* RCT meta-analysis. Marti-Carvajal AJ et al. Cochrane Database Syst Rev. 2017;8:CD006612. JAMA Intern Med 2010 meta-analysis of 8 RCTs.
- **Outcome evidence.** All-cause mortality HR/SD 1.26 (tier B); CVD HR/SD 1.14 (tier B); Cancer (documented, no pooled HR); Dementia HR/SD 1.12 (tier B); Frailty native HR 1.25 [categorical, tier C]

### Ferritin (as inflammation marker)

*Also: serum ferritin, hyperferritinemia*

- **Mechanism.** Primary intracellular iron-storage protein; circulating ferritin reflects body iron stores but is ALSO a positive acute-phase reactant that rises with inflammation independent of iron. This dual identity makes elevated ferritin a mixed signal of iron overload, hepatocellular injury, and/or systemic inflammation.
- **Measurement.** Immunoassay (chemiluminescent or immunoturbidimetric) on a standard analyzer. — *sample:* venous blood; *approx. cost:* $15-40.
- **Signal quality.** noisy — Interpretation is genuinely ambiguous: an elevated value may indicate iron overload, an acute-phase response, liver disease, malignancy, or metabolic syndrome. Without concurrent CRP, transferrin saturation, and clinical context, a single ferritin is difficult to interpret as an inflammation marker specifically.
- **Major confounders.** iron stores / iron overload (primary determinant); acute-phase response / inflammation; liver disease and alcohol; malignancy; metabolic syndrome and fatty liver; recent transfusion; age and sex
- **Testing cadence.** baseline; interpret alongside CRP and transferrin saturation
- **Modifiability (low).** Phlebotomy / venesection for genuine iron overload; treatment of the underlying inflammatory or hepatic condition when ferritin elevation is reactive. — *effect:* Phlebotomy reliably lowers ferritin when iron-driven; reactive (inflammatory) hyperferritinemia falls only when the inflammatory cause resolves.; *timeframe:* Weeks to months.; *evidence:* mechanistic. No outcome RCT shows mortality benefit from lowering ferritin per se in the general population; iron-reduction trials (e.g., VA cooperative trial in PAD) were largely null.
- **Outcome evidence.** All-cause mortality HR/SD 1.05 (tier B); CVD HR/SD 1.05 (tier B); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)

### suPAR

*Also: soluble urokinase plasminogen activator receptor, suPAR*

- **Mechanism.** Circulating cleaved form of the membrane-bound urokinase plasminogen activator receptor (uPAR), shed mainly from immune and endothelial cells. Reflects chronic immune activation; mechanistically implicated in atherosclerosis and as a pathogenic mediator of kidney injury. Marketed as a stable, broad 'biomarker of immune activation and aging'.
- **Measurement.** ELISA (suPARnostic) or proteomic immunoassay; quantified in plasma or serum. — *sample:* venous blood; *approx. cost:* $30-100.
- **Signal quality.** clean — Comparatively stable analyte with lower within-person and diurnal variability than CRP or IL-6, and less reactive to acute infection -- a key reason it is promoted as a marker of chronic inflammation and biological aging. Main limitations are cost and limited assay standardization / availability.
- **Major confounders.** renal function (suPAR rises as GFR falls); smoking; age; obesity; HIV and chronic infection; chronic inflammatory disease
- **Testing cadence.** baseline; specialist / research use, not yet routine
- **Modifiability (low).** No suPAR-specific therapy; levels fall with smoking cessation and treatment of underlying chronic inflammation/infection. Experimental anti-uPAR strategies exist for chronic kidney disease. — *effect:* Smoking cessation lowers suPAR modestly; no intervention has demonstrated outcome benefit via suPAR lowering.; *timeframe:* Months.; *evidence:* mechanistic. No RCT targets suPAR as a modifiable endpoint in the general population; modifiability inferred from observational determinants.
- **Outcome evidence.** All-cause mortality HR/SD 1.43 (tier B); CVD HR/SD 1.35 (tier B); Cancer HR/SD 1.55 (tier C); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)


## Biological Age Clocks

### Horvath clock

*Also: Horvath 2013 multi-tissue clock, DNAmAge, pan-tissue epigenetic clock, first-generation epigenetic clock*

- **Mechanism.** Horvath's 2013 multi-tissue clock estimates age from methylation at 353 CpG sites and was trained purely to predict chronological age across tissues. Epigenetic age acceleration (the residual of clock age on calendar age) reflects deviation from expected methylation aging, but because the clock was not trained on any health outcome it captures mortality and disease risk only weakly. It is a first-generation clock.
- **Measurement.** Genome-wide DNA methylation array (Illumina 450K/EPIC) on extracted DNA; age computed via the published 353-CpG elastic-net predictor, then residualized on chronological age. — *sample:* venous blood; *approx. cost:* $200-500.
- **Signal quality.** moderate — Array methylation is technically reproducible but the clock has appreciable technical noise at the individual-CpG level; test-retest reliability of age acceleration is only moderate, and batch effects and blood cell composition strongly influence estimates.
- **Major confounders.** blood cell type composition; array batch effects; smoking; DNA extraction method; tissue source
- **Testing cadence.** research use; if tracked, every 1-2 years (within-person change is slow relative to measurement noise)
- **Modifiability (low).** First-generation clock age acceleration is only weakly responsive to intervention; modest reductions reported with caloric restriction (CALERIE) and lifestyle programs, but Horvath specifically shows smaller and less consistent change than DunedinPACE. — *effect:* Small: most RCTs show no significant change in Horvath age acceleration; a few report shifts of well under 1 year.; *timeframe:* Months to years; change is slow and hard to distinguish from measurement noise.; *evidence:* RCT. CALERIE caloric restriction trial (Waziry 2023, Nat Aging); lifestyle and diet RCTs with epigenetic clock substudies.
- **Outcome evidence.** All-cause mortality HR/SD 1.14 (tier A); CVD HR/SD 1.20 (tier B); Cancer HR/SD 1.14 (tier C); Dementia (documented, no pooled HR); Frailty HR/SD 1.00 (tier B)

### Hannum clock

*Also: Hannum 2013 clock, blood epigenetic clock, first-generation epigenetic clock*

- **Mechanism.** The Hannum 2013 clock estimates age from 71 blood methylation sites and, like Horvath's, was trained only to predict chronological age. Its age acceleration captures mortality risk somewhat better than Horvath's but still far less than outcome-trained second/third-generation clocks. It is a first-generation clock.
- **Measurement.** Genome-wide DNA methylation array (450K/EPIC) on blood DNA; age via the 71-CpG predictor, residualized on chronological age. — *sample:* venous blood; *approx. cost:* $200-500.
- **Signal quality.** moderate — Reproducible array chemistry but moderate test-retest reliability of age acceleration; sensitive to blood cell composition and batch effects.
- **Major confounders.** blood cell type composition; array batch effects; smoking; DNA extraction method
- **Testing cadence.** research use; if tracked, every 1-2 years
- **Modifiability (low).** First-generation clock age acceleration shows little consistent response to intervention; small reductions reported in some caloric-restriction and lifestyle trials. — *effect:* Small and inconsistent across RCTs; typically not statistically significant.; *timeframe:* Months to years.; *evidence:* RCT. CALERIE caloric restriction trial; lifestyle/diet RCTs with epigenetic clock substudies.
- **Outcome evidence.** All-cause mortality HR/SD 1.15 (tier A); CVD HR/SD 1.23 (tier B); Cancer HR/SD 1.11 (tier C); Dementia (documented, no pooled HR); Frailty HR/SD 1.00 (tier B)

### PhenoAge

*Also: Levine PhenoAge, DNAm PhenoAge, second-generation epigenetic clock, phenotypic age clock*

- **Mechanism.** PhenoAge is a second-generation epigenetic clock: methylation (513 CpGs) trained to predict a composite 'phenotypic age' derived from clinical chemistry and mortality (albumin, creatinine, glucose, CRP, white-cell measures, etc.). Because it was trained on a mortality-linked phenotype rather than chronological age alone, its age acceleration predicts mortality and disease substantially better than first-generation clocks.
- **Measurement.** Genome-wide DNA methylation array (450K/EPIC) on blood DNA; age via the 513-CpG PhenoAge predictor, residualized on chronological age. (A separate blood-chemistry PhenoAge can be computed without methylation.) — *sample:* venous blood; *approx. cost:* $200-500.
- **Signal quality.** moderate — More outcome-relevant than first-generation clocks but still has measurable technical noise; test-retest reliability of age acceleration is moderate. Sensitive to blood cell composition and batch.
- **Major confounders.** blood cell type composition; array batch effects; acute inflammation/illness; smoking; DNA extraction method
- **Testing cadence.** research/consumer use; if tracked, every 1-2 years. Commercially available via TruDiagnostic TruAge and myDNAge.
- **Modifiability (moderate).** PhenoAge acceleration is modestly responsive to lifestyle intervention; caloric restriction (CALERIE) and some diet/exercise programs slow it, and methylation-PhenoAge change tracks improvements in the underlying clinical-chemistry phenotype. — *effect:* Modest: RCTs report reductions on the order of ~1 year of age acceleration, though not always statistically robust.; *timeframe:* Months to a few years.; *evidence:* RCT. Waziry R et al. Effect of caloric restriction on epigenetic aging (CALERIE). Nat Aging 2023;3:248-257; diet/lifestyle RCTs with PhenoAge substudies.
- **Outcome evidence.** All-cause mortality HR/SD 1.14 (tier A); CVD HR/SD 1.17 (tier B); Cancer HR/SD 1.11 (tier C); Dementia HR/SD 1.05 (tier B); Frailty HR/SD 1.00 (tier B)

### GrimAge / GrimAge2

*Also: DNAm GrimAge, GrimAge2, GrimAgeV2, third-generation epigenetic clock, mortality-trained epigenetic clock*

- **Mechanism.** GrimAge is a third-generation epigenetic clock built in two stages: methylation surrogates of seven plasma proteins and smoking pack-years are first estimated, then combined in a model trained directly on time-to-death. This mortality-trained design makes GrimAge (and the refined GrimAge2) the strongest epigenetic predictor of mortality and age-related disease, consistently outperforming first- and second-generation clocks.
- **Measurement.** Genome-wide DNA methylation array (450K/EPIC) on blood DNA; GrimAge computed from DNAm surrogate-protein and smoking components, then residualized on chronological age. — *sample:* venous blood; *approx. cost:* $200-500.
- **Signal quality.** moderate — Best-validated of the epigenetic clocks for outcome prediction; test-retest reliability of GrimAge age acceleration is good relative to other clocks. Still affected by blood cell composition, batch effects, and (by design) heavily influenced by smoking surrogates.
- **Major confounders.** smoking (built into the clock); blood cell type composition; array batch effects; DNA extraction method
- **Testing cadence.** research/consumer use; if tracked, every 1-2 years. Commercially available via TruDiagnostic TruAge (GrimAge2) and other providers.
- **Modifiability (moderate).** GrimAge acceleration responds to smoking cessation (its biggest lever, by design) and modestly to lifestyle and weight-loss interventions; some supplement and multi-component trials report small reductions. — *effect:* Smoking cessation produces meaningful reductions over years; lifestyle interventions yield small effects (under ~1 year of age acceleration) in RCTs.; *timeframe:* Months to years; smoking-driven component changes slowly after cessation.; *evidence:* RCT. CALERIE (Waziry 2023, Nat Aging); lifestyle/diet RCTs with GrimAge substudies; observational smoking-cessation data.
- **Outcome evidence.** All-cause mortality HR/SD 1.50 (tier A); CVD HR/SD 1.55 (tier A); Cancer HR/SD 1.37 (tier B); Dementia HR/SD 1.30 (tier B); Frailty HR/SD 1.43 (tier B)

### DunedinPACE

*Also: DunedinPoAm, Pace of Aging Computed from the Epigenome, rate-of-aging clock, third-generation epigenetic clock*

- **Mechanism.** DunedinPACE is a third-generation epigenetic clock trained not on chronological age but on the longitudinal RATE of decline across 19 organ-system biomarkers measured repeatedly in the Dunedin birth cohort. It estimates the current pace of biological aging (1.0 = one year of biological aging per calendar year). Its rate-based design makes it sensitive to recent and modifiable changes, and it carries among the strongest mortality and disease evidence of any clock.
- **Measurement.** Genome-wide DNA methylation array (EPIC) on blood DNA; pace computed from the published DunedinPACE algorithm (a single-time-point estimate of the rate of aging). — *sample:* venous blood; *approx. cost:* $200-500.
- **Signal quality.** moderate — Highest test-retest reliability of the major epigenetic clocks (designed for reliability), but still subject to array technical noise, blood cell composition, and batch effects. Being a rate measure, it is more responsive to recent exposures than age-trained clocks.
- **Major confounders.** blood cell type composition; array batch effects; smoking; obesity; DNA extraction method
- **Testing cadence.** research/consumer use; can be tracked more frequently than age-trained clocks (e.g. every 6-12 months) because it captures current pace. Commercially available via TruDiagnostic TruAge.
- **Modifiability (moderate).** As a rate measure, DunedinPACE is the most intervention-responsive clock: caloric restriction (CALERIE) significantly slowed DunedinPACE, and diet, exercise and weight-loss programs show measurable slowing. — *effect:* CALERIE: ~2-3% slowing of DunedinPACE with sustained ~12% caloric restriction; lifestyle programs report similar small-to-moderate slowing.; *timeframe:* Months to ~2 years.; *evidence:* RCT. Waziry R et al. Effect of long-term caloric restriction on DNA methylation measures of biological aging (CALERIE). Nat Aging 2023;3:248-257. doi:10.1038/s43587-022-00357-y.
- **Outcome evidence.** All-cause mortality HR/SD 1.65 (tier A); CVD HR/SD 1.39 (tier A); Cancer (documented, no pooled HR); Dementia HR/SD 1.27 (tier B); Frailty HR/SD 1.00 (tier B)

### Telomere length

*Also: leukocyte telomere length, LTL, TL*

- **Mechanism.** Telomeres are repetitive DNA caps that shorten with each cell division and with oxidative/inflammatory stress; leukocyte telomere length is a marker of replicative history and cumulative cellular stress. Shorter telomeres are associated with higher mortality and age-related disease, but telomere length is a weak, noisy individual-level predictor and is upstream of many confounded exposures.
- **Measurement.** Most commonly qPCR (telomere-to-single-copy-gene ratio); also Southern blot (terminal restriction fragment) and flow-FISH. Methods are not interchangeable and qPCR has high measurement error. — *sample:* venous blood; *approx. cost:* $100-300.
- **Signal quality.** noisy — qPCR telomere measurement has substantial within-method coefficient of variation and poor agreement across labs and platforms; a large fraction of measured between-person variance is technical. This noise attenuates HRs and is the central epistemic caveat for this biomarker.
- **Major confounders.** measurement platform / lab (major); blood cell composition; age; sex; smoking; obesity; chronic inflammation; reverse causation from prevalent disease
- **Testing cadence.** research/consumer use; not recommended for repeat clinical tracking given measurement noise. Commercially offered by some longevity testing companies.
- **Modifiability (low).** Telomere length is only modestly and inconsistently modifiable; some lifestyle programs (exercise, diet, stress reduction) and weight loss report attenuated shortening or small lengthening, but results are mixed and prone to measurement-noise artifacts. — *effect:* Small and uncertain; many RCTs show no significant change. Apparent 'lengthening' is often within measurement error.; *timeframe:* Years; telomere attrition is slow and hard to detect against assay noise.; *evidence:* RCT. Lifestyle-intervention RCTs with telomere substudies (e.g. Ornish 2013 Lancet Oncol pilot); generally low-quality or underpowered evidence.
- **Outcome evidence.** All-cause mortality HR/SD 1.09 (tier A); CVD HR/SD 1.22 (tier A); Cancer HR/SD 1.03 (tier B); Dementia HR/SD 1.06 (tier B); Frailty HR/SD 1.00 (tier C)

### GlycanAge

*Also: glycan age, IgG glycan clock, IgG N-glycosylation age, glycomic age*

- **Mechanism.** GlycanAge is a biological-age estimate derived from the N-glycans attached to immunoglobulin G (IgG). IgG glycosylation acts as a functional switch between pro- and anti-inflammatory antibody activity; an age-related shift toward pro-inflammatory (agalactosylated) glycoforms ('inflammaging') underlies the clock. GlycanAge therefore indexes chronic low-grade inflammation rather than a cell-division or methylation process.
- **Measurement.** Chromatographic or mass-spectrometric profiling of IgG N-glycans (e.g. UHPLC after glycan release); age computed from the relative abundance of glycan structures via the GlycanAge algorithm. — *sample:* capillary blood; *approx. cost:* $300-600.
- **Signal quality.** moderate — IgG glycan profiling is analytically reproducible, but the GlycanAge index has a large age-prediction error (~9-10 years) and is strongly influenced by sex, hormones and acute inflammation; it explains roughly 60% of chronological-age variance.
- **Major confounders.** acute infection / inflammation; sex and sex hormones (estrogen); BMI / adiposity; smoking; autoimmune disease
- **Testing cadence.** consumer/wellness use; tracked every 6-12 months. Commercially available as the GlycanAge test.
- **Standardization (Sparse outcomes).** GlycanAge has good cross-sectional disease associations but very limited prospective hazard-ratio evidence for hard outcomes; most cells are unverified or null because no convertible time-to-event HR with a CI was found.
- **Modifiability (moderate).** GlycanAge is among the more intervention-responsive biological-age measures: substantial weight loss markedly shifts IgG glycosylation toward a younger profile, and estrogen/hormone therapy and lifestyle change also move the index. — *effect:* Large in some studies: extensive weight loss (bariatric surgery) produced multi-year reductions in GlycanAge; hormone therapy and lifestyle programs show smaller shifts.; *timeframe:* Months.; *evidence:* observational. Greto VL et al. Extensive weight loss reduces glycan age by altering IgG N-glycosylation. Int J Obes 2021;45:1521-1531. PMID:33941843; estradiol/GlycanAge studies (PMC7732334).
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD native HR 2.22 [categorical, tier C]; Frailty (documented, no pooled HR)

### Proteomic age clocks

*Also: ProtAge, proteomic aging clock, plasma proteomic clock, protein age gap, organ-specific proteomic clocks*

- **Mechanism.** Proteomic age clocks estimate biological age from hundreds to thousands of plasma proteins (e.g. Olink or SomaScan panels) using machine-learning models trained to predict chronological age. The 'proteomic age gap' (predicted minus chronological age) reflects deviation of the circulating proteome from age norms and captures multi-organ aging; organ-specific variants estimate the age of individual organ systems (brain, heart, kidney, immune).
- **Measurement.** High-throughput plasma proteomics (Olink Explore or SomaScan aptamer panels), then a published proteomic-age model; the age gap is the residual on chronological age. — *sample:* venous blood; *approx. cost:* $300-1000.
- **Signal quality.** moderate — Proteomic platforms are reproducible but the clocks are new, model- and platform-specific, and not yet standardized across studies; cross-platform and cross-cohort comparability is unestablished. Age prediction is accurate (r~0.94) but outcome validation is shallower than for epigenetic clocks.
- **Major confounders.** proteomics platform (Olink vs SomaScan); acute illness / inflammation; kidney function; BMI; model-specific protein selection; fasting status
- **Testing cadence.** research use; not yet a consumer product. If tracked, every 1-2 years.
- **Modifiability (moderate).** Proteomic age gap is plausibly modifiable -- it tracks weight, fitness and metabolic health -- but dedicated intervention trials using a proteomic clock as the primary endpoint are largely lacking. — *effect:* Not well quantified; cross-sectional data link healthier lifestyle to a younger proteomic age, but RCT change estimates are not established.; *timeframe:* Unknown; presumably months to years.; *evidence:* observational. Argentieri MA et al. Nat Med 2024 (lifestyle correlates of ProtAgeGap); no dedicated RCT with a proteomic clock as primary endpoint.
- **Outcome evidence.** All-cause mortality HR/SD 1.60 (tier B); CVD HR/SD 1.50 (tier B); Cancer HR/SD 1.20 (tier C); Dementia HR/SD 1.60 (tier B); Frailty (documented, no pooled HR)


## Hormonal Axes

### Cortisol (morning serum)

*Also: AM cortisol, 8am serum cortisol, basal cortisol*

- **Mechanism.** Cortisol is the principal glucocorticoid output of the HPA axis. Chronically elevated cortisol drives visceral adiposity, insulin resistance, sarcopenia, immunosuppression, hippocampal atrophy and accelerated vascular aging; very low cortisol can reflect adrenal insufficiency. The aging HPA axis tends toward flatter, higher-set cortisol exposure.
- **Measurement.** Morning (07:00-09:00) venous draw, immunoassay or LC-MS/MS; strongly time-of-day dependent so sampling must be standardized. — *sample:* venous blood; *approx. cost:* $15-50.
- **Signal quality.** noisy — Large diurnal and pulsatile variation, acute stress reactivity and venipuncture stress make a single morning value an imprecise index of chronic exposure.
- **Major confounders.** time of sampling; acute psychological/physical stress; exogenous glucocorticoids; oral estrogen/pregnancy (raises CBG); shift work; depression; critical illness
- **Testing cadence.** baseline; repeat only if HPA pathology suspected
- **Modifiability (moderate).** Stress-reduction (mindfulness/MBSR, CBT), improved sleep, regular moderate exercise; treatment of overt Cushing's or adrenal insufficiency where present. — *effect:* Mindfulness/relaxation programs lower cortisol output modestly (~10-25% in elevated groups); aerobic training normalizes diurnal slope.; *timeframe:* 8-12 weeks for measurable change; *evidence:* RCT meta-analysis. Koncz A et al. and Pascoe MC et al. meta-analyses of mindfulness/relaxation interventions on cortisol.
- **Outcome evidence.** All-cause mortality HR/SD 1.31 (tier B); CVD HR/SD 1.18 (tier B); Cancer (documented, no pooled HR); Dementia HR/SD 1.20 (tier C); Frailty HR/SD 1.27 (tier C)

### Cortisol diurnal slope (salivary)

*Also: diurnal cortisol slope, salivary cortisol rhythm, flattened slope*

- **Mechanism.** The diurnal cortisol slope quantifies the normal decline from a morning peak to an evening trough. A flatter (less negative) slope reflects HPA-axis dysregulation and chronic stress load and is associated with inflammation, metabolic dysfunction and accelerated aging.
- **Measurement.** Multiple same-day saliva samples (wake, +30-45 min, evening) by passive drool or swab; slope computed by regression across the day. — *sample:* saliva; *approx. cost:* $60-150.
- **Signal quality.** moderate — Captures rhythm rather than a single level, but depends on adherence to sampling times and is variable day-to-day; multi-day sampling improves reliability.
- **Major confounders.** sampling-time adherence; wake time/shift work; smoking; acute stress; exogenous glucocorticoids; oral contraceptives
- **Testing cadence.** baseline; repeat over 2-3 days for reliability
- **Standardization (Non-standard).** Slope is a rate-of-change metric (units of cortisol per hour); a flatter slope is the risk direction and per-SD estimates vary widely by sampling protocol.
- **Modifiability (moderate).** Sleep regularization, light exposure timing, stress reduction, aerobic exercise; treating depression where present. — *effect:* Behavioral programs can steepen a blunted slope modestly; effect sizes small and variable.; *timeframe:* 8-12 weeks; *evidence:* observational. Adam EK et al. systematic review of diurnal cortisol slopes and health outcomes, Psychoneuroendocrinology 2017.
- **Outcome evidence.** All-cause mortality HR/SD 1.31 (tier B); CVD HR/SD 1.30 (tier C); Cancer native HR 2.08 [categorical, tier C]; Frailty HR/SD 1.23 (tier C)

### Hair cortisol

*Also: hair cortisol concentration, HCC*

- **Mechanism.** Hair cortisol concentration integrates systemic cortisol exposure over the preceding 1-3 months, providing a retrospective index of chronic HPA-axis activation that avoids diurnal and acute-stress noise.
- **Measurement.** Proximal 1-3 cm hair segment, cortisol extracted and quantified by immunoassay or LC-MS/MS. — *sample:* hair; *approx. cost:* $60-150.
- **Signal quality.** moderate — Reflects chronic exposure well, but affected by hair treatments, washing frequency, hair color and growth rate; assays not fully standardized.
- **Major confounders.** hair dye/bleaching; washing frequency; hair growth rate; topical/systemic corticosteroids; season
- **Testing cadence.** baseline; repeat at 3-6 months if monitoring chronic stress
- **Standardization (Cross-sectional OR).** The hair-cortisol cardiovascular cell is a cross-sectional odds ratio for prevalent disease across cortisol categories, not an incidence hazard ratio, so it is not placed on the per-SD scale.
- **Modifiability (moderate).** Sustained stress reduction, sleep, exercise; effective only over multi-month horizons given the integrative window. — *effect:* Modest reductions with sustained behavioral change; data limited.; *timeframe:* 3+ months (reflects cumulative exposure); *evidence:* observational. Stalder T et al. Stress-related and basic determinants of hair cortisol in humans: a meta-analysis, Psychoneuroendocrinology 2017.
- **Outcome evidence.** All-cause mortality HR/SD 1.20 (tier C); CVD native HR 2.70 [categorical, tier C]; Frailty HR/SD 1.30 (tier C)

### DHEA-S

*Also: dehydroepiandrosterone sulfate, DHEAS*

- **Mechanism.** DHEA-S is the most abundant circulating steroid and the sulfated reservoir of DHEA, an adrenal androgen precursor. It declines steeply with age ('adrenopause'); low levels are a marker of biological aging and have been linked to reduced anabolic tone, immune function and vascular health.
- **Measurement.** Immunoassay or LC-MS/MS on venous serum; stable through the day so timing is not critical. — *sample:* venous blood; *approx. cost:* $30-70.
- **Signal quality.** clean — Long half-life, minimal diurnal variation and high concentration make DHEA-S a stable, reproducible analyte.
- **Major confounders.** age (strong); exogenous DHEA/steroids; critical illness; anorexia; smoking
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (low).** DHEA supplementation raises DHEA-S, but RCTs show little benefit on hard outcomes, body composition or cognition in healthy older adults. — *effect:* Supplementation normalizes levels; clinically meaningful outcome benefit not demonstrated.; *timeframe:* Weeks to raise levels; outcome benefit unproven; *evidence:* RCT meta-analysis. Nair KS et al. DHEA in elderly women and men, NEJM 2006;355:1647-1659; Cochrane/meta-analyses of DHEA supplementation.
- **Outcome evidence.** All-cause mortality HR/SD 1.21 (tier B); CVD HR/SD 1.15 (tier C); Cancer HR/SD 1.05 (tier C); Dementia HR/SD 1.10 (tier C); Frailty HR/SD 1.30 (tier C)

### Cortisol:DHEA-S ratio

*Also: cortisol/DHEAS ratio, catabolic-anabolic balance index*

- **Mechanism.** The cortisol:DHEA-S ratio indexes the balance between catabolic glucocorticoid tone and anabolic adrenal-androgen tone. A rising ratio with age reflects shift toward catabolism and has been linked to immunosenescence, sarcopenia and poorer stress resilience.
- **Measurement.** Computed from concurrently measured serum cortisol and DHEA-S; both from a single morning venous draw. — *sample:* venous blood; *approx. cost:* $45-100.
- **Signal quality.** moderate — Inherits the diurnal noise of cortisol; the ratio has no standardized reference range and units are arbitrary.
- **Major confounders.** time of sampling; acute stress; exogenous steroids; age; critical illness
- **Testing cadence.** baseline; research-oriented use
- **Standardization (Composite ratio).** A derived ratio of two hormones with no standard SD; per-SD estimates depend entirely on how each component is scaled.
- **Modifiability (low).** Lower the ratio indirectly via stress reduction, sleep and exercise (lowering cortisol); DHEA supplementation raises the denominator but without proven outcome benefit. — *effect:* Modest, indirect; no validated intervention targets the ratio itself.; *timeframe:* 8-12 weeks for component changes; *evidence:* observational. Inferred from cortisol and DHEA-S intervention literature.
- **Outcome evidence.** All-cause mortality native HR 1.25 [categorical, tier C]; CVD native HR 1.20 [categorical, tier C]; Frailty native HR 1.35 [categorical, tier C]

### TSH

*Also: thyroid-stimulating hormone, thyrotropin*

- **Mechanism.** TSH is the pituitary signal driving thyroid hormone output and is the most sensitive index of thyroid status. Both subclinical hyperthyroidism (low TSH) and subclinical hypothyroidism (high TSH) perturb cardiac, metabolic and skeletal physiology; thyroid status shifts with age.
- **Measurement.** Third-generation immunoassay on venous serum; mild diurnal variation. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** clean — Well-standardized, sensitive assay; log-normal distribution and modest diurnal/seasonal variation are the main caveats.
- **Major confounders.** non-thyroidal illness; biotin supplements (assay interference); amiodarone/lithium; pregnancy; time of day; iodine status
- **Testing cadence.** baseline + every 1-2 years; sooner if abnormal
- **Standardization (U-shaped).** Mortality and CVD rise at both low (hyperthyroid) and high (hypothyroid) TSH; the relationship is U/J-shaped and not summarizable as one per-SD slope.
- **Modifiability (high).** Levothyroxine for hypothyroidism; antithyroid drugs/radioiodine for hyperthyroidism. RCTs (TRUST) show levothyroxine for mild subclinical hypothyroidism does not improve symptoms or outcomes in older adults. — *effect:* TSH fully normalizable with treatment; outcome benefit clear only for overt disease and TSH >10.; *timeframe:* 6-8 weeks to re-equilibrate after dose change; *evidence:* RCT. Stott DJ et al. TRUST trial, NEJM 2017;376:2534-2544; ATA guidelines.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD native HR 1.58 [categorical, tier A]; Cancer HR/SD 1.10 (tier C); Dementia (documented, no pooled HR); Frailty native HR 1.20 [categorical, tier C]

### Free T4

*Also: free thyroxine, FT4*

- **Mechanism.** Free T4 is the unbound, biologically available fraction of the main thyroid hormone. Even within the reference range, higher free T4 reflects greater metabolic and cardiac drive and is associated with atrial fibrillation, bone loss and accelerated aging phenotypes.
- **Measurement.** Immunoassay (or equilibrium dialysis LC-MS/MS) on venous serum. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** moderate — Assay estimates of free hormone are method-dependent and affected by binding-protein abnormalities and non-thyroidal illness.
- **Major confounders.** non-thyroidal illness; binding-protein changes (pregnancy, estrogen); heparin; amiodarone; assay method
- **Testing cadence.** baseline + with TSH every 1-2 years
- **Modifiability (low).** Free T4 is treated only in the context of overt thyroid disease (antithyroid drugs lower it; levothyroxine raises it). Within-range variation is not a clinical treatment target. — *effect:* Fully modifiable by thyroid-directed therapy when disease is present; not modified in euthyroid people.; *timeframe:* Weeks after dose change; *evidence:* RCT. ATA thyroid guidelines; TRUST trial context.
- **Outcome evidence.** All-cause mortality HR/SD 1.16 (tier B); CVD HR/SD 1.24 (tier B); Dementia HR/SD 1.04 (tier C); Frailty HR/SD 1.18 (tier C)

### Free T3

*Also: free triiodothyronine, FT3*

- **Mechanism.** Free T3 is the most metabolically active thyroid hormone, generated largely by peripheral deiodination of T4. Low free T3 is a hallmark of non-thyroidal illness syndrome and a strong marker of illness severity and frailty in older and hospitalized populations.
- **Measurement.** Immunoassay on venous serum; free-hormone estimates are method-dependent. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** noisy — Strongly depressed by acute and chronic illness independent of thyroid disease, making it a confounded index of intrinsic thyroid status.
- **Major confounders.** non-thyroidal illness (major); caloric restriction/fasting; age; drugs (glucocorticoids, beta-blockers, amiodarone); assay method
- **Testing cadence.** baseline; interpret with TSH and clinical context
- **Standardization (U-shaped).** Low free T3 marks non-thyroidal illness ('low-T3 syndrome') and predicts mortality, but very high free T3 is also harmful; the association is non-monotonic and confounded by reverse causation.
- **Modifiability (low).** No intervention targets free T3 in non-thyroidal illness; T3 replacement in low-T3 syndrome has not improved outcomes in RCTs. Free T3 normalizes when underlying illness resolves. — *effect:* Not a treatment target; reflects underlying health.; *timeframe:* Tracks recovery from illness; *evidence:* RCT. Trials of T3 supplementation in critical illness and heart failure (no consistent benefit).
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD native HR 2.40 [categorical, tier C]; Frailty native HR 1.50 [categorical, tier C]

### Reverse T3

*Also: rT3, reverse triiodothyronine*

- **Mechanism.** Reverse T3 is the metabolically inactive product of T4 deiodination. Its production rises during illness, starvation and stress as the body downregulates active T3, so elevated rT3 is a sensitive marker of non-thyroidal illness and catabolic state.
- **Measurement.** Immunoassay or LC-MS/MS on venous serum; not routinely available and assay variability is high. — *sample:* venous blood; *approx. cost:* $60-150.
- **Signal quality.** noisy — Assays poorly standardized; the value is dominated by acute and chronic illness, making it an indirect, confounded biomarker.
- **Major confounders.** non-thyroidal illness (major); fasting/caloric restriction; glucocorticoids/amiodarone; critical illness; assay variability
- **Testing cadence.** not routinely indicated
- **Modifiability (fixed).** No intervention targets reverse T3; it normalizes when underlying illness or caloric deprivation resolves. — *effect:* Not a treatment target.; *timeframe:* Tracks recovery; *evidence:* mechanistic. Non-thyroidal illness syndrome literature.
- **Outcome evidence.** All-cause mortality HR/SD 1.13 (tier C); Frailty HR/SD 1.17 (tier C)

### Thyroid antibodies (TPO/Tg)

*Also: TPOAb, anti-thyroid peroxidase antibody, thyroglobulin antibody, TgAb*

- **Mechanism.** Antibodies to thyroid peroxidase (TPO) and thyroglobulin (Tg) mark autoimmune thyroiditis (Hashimoto's) and predict progression to overt hypothyroidism. They index thyroid autoimmunity rather than thyroid hormone level.
- **Measurement.** Immunoassay on venous serum; reported as a titre with a positivity cutoff. — *sample:* venous blood; *approx. cost:* $25-70.
- **Signal quality.** moderate — Reliable for detecting autoimmunity, but titres fluctuate and a positive result is common in older women without clinical disease.
- **Major confounders.** age; female sex; other autoimmune disease; iodine status; assay/cutoff differences
- **Testing cadence.** once to characterize cause of thyroid dysfunction; not routinely repeated
- **Standardization (Categorical).** Thyroid antibodies are interpreted as positive vs negative (autoimmune thyroiditis); there is no meaningful per-SD scale.
- **Modifiability (low).** Antibody titres are not a treatment target; selenium supplementation modestly lowers TPOAb in some trials without proven clinical benefit. Levothyroxine treats the resulting hypothyroidism. — *effect:* Selenium lowers TPOAb titres modestly; no demonstrated outcome benefit.; *timeframe:* 3-6 months for titre change; *evidence:* RCT meta-analysis. Meta-analyses of selenium supplementation in autoimmune thyroiditis.
- **Outcome evidence.** All-cause mortality native HR 1.02 [categorical, tier C]; CVD native HR 1.15 [categorical, tier C]; Cancer (documented, no pooled HR)

### Total testosterone

*Also: serum testosterone, total T*

- **Mechanism.** Total testosterone is the principal circulating androgen, mostly bound to SHBG and albumin. It supports muscle mass, bone density, erythropoiesis and metabolic health; levels decline gradually with age in men, and low testosterone is both a cause and a consequence of poor health.
- **Measurement.** LC-MS/MS (preferred) or immunoassay on a morning venous draw; diurnal and day-to-day variation require standardized timing. — *sample:* venous blood; *approx. cost:* $25-80.
- **Signal quality.** moderate — Diurnal variation, acute illness suppression and assay differences (immunoassay vs LC-MS/MS) reduce reliability of single measurements.
- **Major confounders.** time of day; acute illness; obesity (lowers via SHBG); opioids/glucocorticoids; sleep deprivation; assay method
- **Testing cadence.** morning sample, confirm low values on a second draw
- **Standardization (U-shaped).** In men, mortality is lowest at mid-normal testosterone and rises at both low and very high concentrations; per-SD HRs apply best within the low-to-normal range.
- **Modifiability (moderate).** Testosterone replacement therapy in confirmed hypogonadism; weight loss and treating sleep apnea raise endogenous levels in obese men. — *effect:* TRT restores levels to mid-normal; ~10-15% weight loss can raise total testosterone substantially in obese men.; *timeframe:* Weeks (TRT); 3-6 months (weight loss); *evidence:* RCT. TRAVERSE trial (Lincoff 2023, NEJM); TTrials (Snyder 2016, NEJM); weight-loss trials.
- **Outcome evidence.** All-cause mortality HR/SD 1.14 (tier A); CVD native HR 1.25 [categorical, tier A]; Cancer HR/SD 1.00 (tier B); Dementia HR/SD 1.10 (tier C); Frailty HR/SD 1.30 (tier B)

### Free testosterone

*Also: free T, calculated free testosterone, bioavailable testosterone*

- **Mechanism.** Free testosterone is the unbound, immediately bioavailable androgen fraction, usually calculated from total testosterone, SHBG and albumin. It is considered a more physiologically relevant index of androgen action, especially when SHBG is abnormal (obesity, aging).
- **Measurement.** Calculated from total testosterone, SHBG and albumin (Vermeulen equation), or equilibrium dialysis (reference method); morning sample. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** moderate — Calculated values depend on SHBG assay accuracy and the chosen equation; equilibrium dialysis is accurate but rarely used.
- **Major confounders.** time of day; SHBG assay accuracy; obesity; acute illness; calculation method
- **Testing cadence.** morning sample; confirm low values
- **Standardization (Category contrast).** The free-testosterone frailty cell is a low-vs-higher category contrast in older men, not a continuous per-SD estimate.
- **Modifiability (moderate).** Testosterone replacement in hypogonadism; weight loss raises free testosterone partly by lowering SHBG-bound fraction. — *effect:* TRT restores free testosterone; weight loss produces modest gains.; *timeframe:* Weeks (TRT); months (weight loss); *evidence:* RCT. TTrials (Snyder 2016, NEJM); TRAVERSE (Lincoff 2023, NEJM).
- **Outcome evidence.** All-cause mortality HR/SD 1.13 (tier B); CVD HR/SD 1.15 (tier C); Dementia HR/SD 1.15 (tier C); Frailty native HR 1.34 [categorical, tier B]

### SHBG

*Also: sex hormone-binding globulin*

- **Mechanism.** SHBG is a hepatic glycoprotein that binds sex steroids and regulates their bioavailability. Low SHBG is a sensitive marker of insulin resistance, hepatic fat and metabolic syndrome; very high SHBG accompanies aging, low body weight and frailty.
- **Measurement.** Immunoassay on venous serum; stable through the day. — *sample:* venous blood; *approx. cost:* $25-60.
- **Signal quality.** clean — Low diurnal variation and a robust assay; interpretation is the main complexity given its U-shaped associations.
- **Major confounders.** insulin resistance/obesity (lowers); oral estrogen/pregnancy (raises); thyroid status; liver disease; age
- **Testing cadence.** baseline; with sex-steroid panels
- **Standardization (U-shaped).** SHBG shows a U-shaped relationship with mortality - low SHBG marks insulin resistance, very high SHBG marks frailty/illness - so a single per-SD slope is only approximate.
- **Modifiability (moderate).** Weight loss and improved insulin sensitivity raise low SHBG; reducing hepatic fat is the main lever. — *effect:* Significant weight loss can raise SHBG 20-50% in insulin-resistant individuals.; *timeframe:* 3-6 months; *evidence:* RCT. Weight-loss and lifestyle-intervention trials reporting SHBG change.
- **Outcome evidence.** All-cause mortality HR/SD 1.15 (tier B); CVD HR/SD 0.90 (tier B); Cancer native HR 1.20 [categorical, tier C]; Frailty HR/SD 1.25 (tier C)

### Estradiol

*Also: E2, 17-beta-estradiol*

- **Mechanism.** Estradiol is the dominant estrogen, central to reproductive physiology and with vascular, skeletal, metabolic and neurocognitive actions. Its associations with longevity outcomes depend heavily on sex, menopausal status and the source (endogenous vs exogenous).
- **Measurement.** LC-MS/MS (preferred for low postmenopausal/male levels) or immunoassay on venous serum. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** noisy — Immunoassays are unreliable at the low concentrations seen in men and postmenopausal women; menstrual-cycle phase causes large variation in premenopausal women.
- **Major confounders.** sex and menopausal status; menstrual-cycle phase; exogenous hormones; adiposity (aromatization); assay sensitivity
- **Testing cadence.** context-dependent; cycle-timed in premenopausal women
- **Modifiability (high).** Menopausal hormone therapy raises estradiol; aromatase inhibitors lower it. For longevity outcomes, the route, dose and timing dominate the risk-benefit balance. — *effect:* Hormone therapy fully restores premenopausal-range estradiol.; *timeframe:* Days to weeks; *evidence:* RCT. WHI and ELITE/KEEPS hormone-therapy trials.
- **Outcome evidence.** All-cause mortality HR/SD 1.15 (tier C); CVD (documented, no pooled HR); Cancer HR/SD 1.16 (tier A)

### Progesterone

*Also: P4*

- **Mechanism.** Progesterone is the principal progestogen, central to the luteal phase and pregnancy, with secondary actions on the brain (neurosteroid), breast and vasculature. Endogenous levels are highly cycle-dependent in premenopausal women and very low otherwise.
- **Measurement.** Immunoassay or LC-MS/MS on venous serum; must be timed to cycle phase to be interpretable. — *sample:* venous blood; *approx. cost:* $25-70.
- **Signal quality.** noisy — Extreme within-cycle variation and pulsatile secretion make a single measurement weakly informative for chronic-disease risk.
- **Major confounders.** menstrual-cycle phase (major); pregnancy; exogenous progestogens; menopausal status
- **Testing cadence.** cycle-timed (luteal phase) when assessing ovulation
- **Modifiability (low).** Exogenous progesterone/progestins are used in hormone therapy and contraception; endogenous progesterone is not a longevity treatment target. — *effect:* Exogenous dosing fully controls levels.; *timeframe:* Days; *evidence:* RCT. Hormone-therapy trials (WHI combined arm).
- **Outcome evidence.** Cancer HR/SD 1.00 (tier C)

### LH

*Also: luteinizing hormone*

- **Mechanism.** LH is a pituitary gonadotropin driving gonadal steroidogenesis. Rising LH signals declining gonadal feedback (menopause, primary hypogonadism, aging testes); elevated gonadotropins have been hypothesized to contribute to neurodegeneration.
- **Measurement.** Immunoassay on venous serum; pulsatile secretion means single values are noisy. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** noisy — Pulsatile secretion (90-120 min pulses) and cycle-phase dependence limit the information in a single measurement.
- **Major confounders.** pulsatile secretion; menstrual-cycle phase; menopausal status; GnRH analogues; pituitary disease
- **Testing cadence.** context-dependent; not routine for risk stratification
- **Modifiability (low).** LH is suppressed by sex-steroid replacement or GnRH analogues; it is not an independent longevity treatment target. — *effect:* Fully suppressible pharmacologically.; *timeframe:* Weeks; *evidence:* mechanistic. Endocrine physiology of the HPG axis.
- **Outcome evidence.** All-cause mortality HR/SD 1.10 (tier C); Dementia (documented, no pooled HR)

### FSH

*Also: follicle-stimulating hormone*

- **Mechanism.** FSH is a pituitary gonadotropin driving gametogenesis; it rises sharply with ovarian aging and the menopausal transition. Beyond the gonad, FSH has been linked experimentally to bone resorption, fat accumulation and possibly atherosclerosis.
- **Measurement.** Immunoassay on venous serum; less pulsatile than LH so single values are somewhat more stable. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** moderate — More stable than LH but still cycle- and menopausal-status dependent; rises steeply across the menopause transition.
- **Major confounders.** menopausal status (major); menstrual-cycle phase; GnRH analogues; ovarian reserve; pituitary disease
- **Testing cadence.** context-dependent; used to confirm menopausal status
- **Modifiability (low).** FSH is suppressed by estrogen/hormone therapy; experimental FSH-blocking antibodies are preclinical. Not currently a clinical treatment target. — *effect:* Hormone therapy lowers FSH substantially.; *timeframe:* Weeks; *evidence:* mechanistic. Menopause hormone-therapy pharmacology; preclinical anti-FSH studies.
- **Outcome evidence.** All-cause mortality HR/SD 1.05 (tier C); CVD HR/SD 1.12 (tier C); Frailty HR/SD 1.15 (tier C)

### IGF-1

*Also: insulin-like growth factor 1, somatomedin C*

- **Mechanism.** IGF-1 mediates most growth-hormone actions and drives the conserved GH/IGF-1/insulin signalling pathway central to longevity biology. Low IGF-1 signalling extends lifespan in model organisms but in humans low IGF-1 also marks frailty; high IGF-1 promotes cell proliferation and cancer risk.
- **Measurement.** Immunoassay or LC-MS/MS on venous serum; relatively stable through the day, strongly age-dependent. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** moderate — Stable analyte but strongly age- and nutrition-dependent; assays are not fully harmonized, so age-specific reference ranges are essential.
- **Major confounders.** age (strong); nutritional status/protein intake; liver disease; GH disorders; assay standardization
- **Testing cadence.** baseline + every 1-2 years; interpret against age norms
- **Standardization (U-shaped).** IGF-1 has a U-shaped mortality relationship - both low and high levels raise risk - and opposite directions for cancer (higher = risk) versus frailty (lower = risk).
- **Modifiability (moderate).** IGF-1 is lowered by protein/calorie restriction and raised by higher protein intake, GH therapy and resistance exercise. There is no consensus optimal target given U-shaped risk. — *effect:* Protein restriction lowers IGF-1 ~10-25%; GH therapy or high protein raises it.; *timeframe:* Weeks to months; *evidence:* RCT. Dietary protein/IGF-1 trials; Levine ME et al. low-protein intake and IGF-1, Cell Metab 2014.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD native HR 1.15 [categorical, tier B]; Cancer HR/SD 1.09 (tier A); Dementia HR/SD 1.12 (tier C); Frailty HR/SD 1.22 (tier B)

### IGFBP-3

*Also: insulin-like growth factor binding protein 3*

- **Mechanism.** IGFBP-3 is the main carrier of circulating IGF-1, regulating its half-life and bioavailability, and also has IGF-independent pro-apoptotic actions. It is interpreted largely in relation to IGF-1 (the IGF-1/IGFBP-3 ratio indexes free IGF-1).
- **Measurement.** Immunoassay on venous serum; stable analyte. — *sample:* venous blood; *approx. cost:* $30-80.
- **Signal quality.** moderate — Stable and reproducible, but biologically meaningful chiefly in the context of IGF-1; independent interpretation is limited.
- **Major confounders.** age; nutritional status; GH disorders; liver disease; co-variation with IGF-1
- **Testing cadence.** baseline; with IGF-1
- **Modifiability (low).** IGFBP-3 rises with GH therapy and varies with nutrition; it is not an independent treatment target. — *effect:* Modest changes with GH/nutritional manipulation.; *timeframe:* Weeks to months; *evidence:* observational. GH-therapy and nutritional studies reporting IGFBP-3.
- **Outcome evidence.** All-cause mortality HR/SD 1.10 (tier C); Cancer HR/SD 0.95 (tier C); Frailty HR/SD 1.15 (tier C)

### Growth hormone

*Also: GH, somatotropin, human growth hormone*

- **Mechanism.** Growth hormone, secreted in pulses by the pituitary, drives IGF-1 production and has direct lipolytic and anti-insulin actions. Both GH excess (acromegaly) and severe deficiency shorten life, while reduced GH/IGF-1 signalling extends lifespan in model organisms - a central paradox of longevity biology.
- **Measurement.** Random serum GH is rarely useful; assessment requires stimulation or suppression tests, or measuring IGF-1 as a stable proxy. — *sample:* venous blood; *approx. cost:* $30-200.
- **Signal quality.** noisy — Pulsatile secretion (peaks mostly in slow-wave sleep) and a ~15-minute half-life mean a single random GH value carries little information; dynamic testing is needed.
- **Major confounders.** pulsatile secretion; sleep/exercise/fasting (acute peaks); age; adiposity (blunts GH); assay standardization
- **Testing cadence.** not measured as a single value; dynamic testing only when GH disorder suspected
- **Modifiability (low).** GH excess is treated with surgery, somatostatin analogues or GH-receptor antagonists; GH deficiency with GH replacement. GH supplementation in healthy older adults is not recommended. — *effect:* Disease-state GH is fully treatable; in healthy aging, GH supplementation gives small body-composition change with notable adverse effects.; *timeframe:* Months; *evidence:* RCT meta-analysis. Liu H et al. Ann Intern Med 2007;146:104-115; acromegaly treatment guidelines.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD (documented, no pooled HR); Cancer (documented, no pooled HR); Frailty (documented, no pooled HR)

### Melatonin

*Also: N-acetyl-5-methoxytryptamine*

- **Mechanism.** Melatonin, secreted nocturnally by the pineal gland, is the principal circadian signal and a direct antioxidant. Nocturnal melatonin output declines with age and with light-at-night exposure; low or mistimed melatonin is implicated in circadian disruption, metabolic dysfunction and possibly hormone-sensitive cancer.
- **Measurement.** Overnight urinary 6-sulfatoxymelatonin (aMT6s, an integrated measure) or timed salivary/plasma melatonin; a single daytime value is uninformative. — *sample:* urine; *approx. cost:* $40-150.
- **Signal quality.** noisy — Strong circadian variation means the result depends entirely on sampling timing; light exposure, age and beta-blockers strongly affect levels.
- **Major confounders.** sampling timing; light-at-night/shift work; age; beta-blockers/NSAIDs; exogenous melatonin
- **Testing cadence.** timed overnight collection when assessed at all
- **Modifiability (moderate).** Light hygiene (bright light by day, darkness at night), consistent sleep timing and limiting light-at-night support endogenous melatonin; exogenous melatonin supplements circadian signalling but does not restore the endogenous rhythm. — *effect:* Behavioral measures modestly preserve nocturnal output; exogenous melatonin reliably shifts circadian timing.; *timeframe:* Days to weeks; *evidence:* RCT. Circadian light-therapy and melatonin RCTs.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.12 (tier C); Cancer HR/SD 1.06 (tier B)


## Body Composition & Anthropometrics

### DEXA total body fat %

*Also: DXA body fat percentage, percent body fat, %BF, total adiposity*

- **Mechanism.** Total body fat percentage measured by dual-energy X-ray absorptiometry quantifies whole-body adiposity directly, separating fat mass from lean and bone. Excess adiposity drives insulin resistance, chronic low-grade inflammation, adipokine dysregulation and ectopic lipid deposition, all of which accelerate cardiometabolic disease and biological aging; however very low fat mass in older adults often signals occult illness, wasting and frailty.
- **Measurement.** Whole-body dual-energy X-ray absorptiometry scan; three-compartment estimation of fat, lean and bone mass with low radiation dose. — *sample:* imaging; *approx. cost:* $75-200.
- **Signal quality.** clean — DXA is a precise, low-variability reference method for body fat with excellent reproducibility; minor between-machine and software-version differences exist but within-person tracking is reliable.
- **Major confounders.** reverse causation (illness-related fat loss in older adults); hydration status of lean compartment; scanner/software differences; age- and sex-dependent reference ranges; fat distribution not captured by total %
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (high).** Caloric restriction / dietary change plus aerobic and resistance exercise; GLP-1 receptor agonists (semaglutide, tirzepatide) and bariatric surgery for substantial reductions. — *effect:* Lifestyle programs reduce total body fat ~3-8% of body weight; GLP-1/GIP agonists reduce fat mass 10-20%+; bariatric surgery 20-30%+.; *timeframe:* 3-12 months for meaningful change; *evidence:* RCT meta-analysis. Look AHEAD lifestyle trial; STEP (semaglutide) and SURMOUNT (tirzepatide) RCT programs; bariatric surgery cohorts.
- **Outcome evidence.** All-cause mortality HR/SD 1.00 (tier B); CVD HR/SD 1.21 (tier B); Cancer HR/SD 1.10 (tier B); Dementia HR/SD 1.00 (tier C); Frailty HR/SD 1.00 (tier C)

### Visceral adipose tissue (VAT)

*Also: VAT, visceral fat, intra-abdominal fat, visceral fat area, visceral fat mass*

- **Mechanism.** Visceral adipose tissue is metabolically active fat surrounding the abdominal organs. It drains via the portal vein, delivering free fatty acids and pro-inflammatory adipokines (TNF-alpha, IL-6, resistin) directly to the liver, driving hepatic insulin resistance, dyslipidemia, systemic inflammation and ectopic lipid deposition. VAT is the most pathogenic fat depot and predicts cardiometabolic disease beyond total adiposity or BMI.
- **Measurement.** Direct quantification by DXA visceral fat algorithm, abdominal CT or MRI; estimated by some bioimpedance devices. CT/MRI are reference standards. — *sample:* imaging; *approx. cost:* $75-500.
- **Signal quality.** clean — CT/MRI and DXA visceral algorithms quantify VAT precisely with low measurement error; bioimpedance estimates are noisier. Strong, depot-specific pathophysiologic signal.
- **Major confounders.** sex (men accumulate more VAT); age; ethnicity (Asian populations higher VAT at given BMI); menopausal status; reverse causation in advanced illness; measurement modality differences
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (high).** Caloric restriction, low-refined-carbohydrate diet, and especially aerobic exercise preferentially mobilize visceral fat; GLP-1/GIP agonists and bariatric surgery produce large VAT reductions. — *effect:* Aerobic exercise without weight loss can reduce VAT ~6-10%; combined lifestyle 15-30%; GLP-1 agonists and bariatric surgery reduce VAT substantially more than subcutaneous fat.; *timeframe:* 3-6 months; *evidence:* RCT meta-analysis. Meta-analyses of exercise and visceral fat (Vissers D et al. PLoS One 2013); STEP/SURMOUNT GLP-1 imaging substudies.
- **Outcome evidence.** All-cause mortality HR/SD 1.21 (tier B); CVD HR/SD 1.32 (tier A); Cancer HR/SD 1.20 (tier C); Dementia HR/SD 1.15 (tier C); Frailty HR/SD 1.18 (tier C)

### Bone mineral density

*Also: BMD, DXA bone density, femoral neck BMD, T-score*

- **Mechanism.** Bone mineral density reflects skeletal mineral mass per unit area and is the principal determinant of bone strength and fracture resistance. Low BMD (osteopenia/osteoporosis) is both a direct cause of fragility fractures and a marker of systemic aging, low physical activity, undernutrition, inflammation and hormonal decline. Fractures, particularly hip fractures, are followed by sharply elevated mortality and loss of independence.
- **Measurement.** Dual-energy X-ray absorptiometry of the lumbar spine and proximal femur; reported as areal density (g/cm2) and as T-score / Z-score relative to a reference population. — *sample:* imaging; *approx. cost:* $75-150.
- **Signal quality.** clean — DXA BMD is a precise, well-standardized measurement with low short-term variability and validated reference databases; degenerative spine changes can spuriously raise lumbar BMD in older adults.
- **Major confounders.** spinal osteoarthritis/calcification (spuriously raises lumbar BMD); age and sex reference effects; body size (areal BMD partly size-dependent); glucocorticoid use; reverse causation (illness-related bone loss)
- **Testing cadence.** baseline + every 2 years (more frequent on treatment)
- **Modifiability (high).** Resistance and weight-bearing exercise, adequate calcium and vitamin D, and pharmacotherapy (bisphosphonates, denosumab, anabolic agents teriparatide/romosozumab) for osteoporosis. — *effect:* Antiresorptive drugs raise BMD ~3-8% over 1-3 years and cut fracture risk 30-70%; exercise produces smaller (~1-3%) BMD gains but improves fall risk.; *timeframe:* 1-3 years for BMD change; fracture-risk reduction within 1 year on potent agents; *evidence:* RCT meta-analysis. FIT (alendronate), FREEDOM (denosumab), FRAME (romosozumab) RCTs; Cochrane reviews of exercise and BMD.
- **Outcome evidence.** All-cause mortality HR/SD 1.20 (tier B); CVD HR/SD 1.24 (tier B); Cancer (documented, no pooled HR); Dementia HR/SD 1.21 (tier C); Frailty HR/SD 1.40 (tier A)

### Waist circumference

*Also: WC, abdominal circumference, central adiposity*

- **Mechanism.** Waist circumference is a simple anthropometric proxy for central (abdominal) adiposity, capturing both visceral and subcutaneous abdominal fat. Central fat is the metabolically harmful depot, driving insulin resistance, dyslipidemia, hypertension and systemic inflammation. Waist circumference predicts cardiometabolic disease and mortality more strongly and more monotonically than BMI, especially within normal-BMI ranges.
- **Measurement.** Tape measure at the midpoint between the lowest rib and iliac crest (or at the umbilicus per protocol); requires standardized technique. — *sample:* in-office device; *approx. cost:* $0-5.
- **Signal quality.** moderate — Inexpensive and predictive, but measurement site and observer technique vary, and it does not separate visceral from subcutaneous fat; standardized protocols improve reliability.
- **Major confounders.** measurement site/protocol differences; observer technique; respiratory phase / posture; sex and ethnicity (different risk thresholds); reverse causation in advanced illness
- **Testing cadence.** baseline + annual
- **Modifiability (high).** Caloric restriction, reduced refined-carbohydrate intake, aerobic and resistance exercise; GLP-1/GIP agonists and bariatric surgery for large reductions. — *effect:* Lifestyle programs reduce waist circumference ~3-7 cm; GLP-1/GIP agonists 7-15+ cm; bariatric surgery substantially more.; *timeframe:* 3-6 months; *evidence:* RCT meta-analysis. Look AHEAD and Diabetes Prevention Program lifestyle trials; STEP/SURMOUNT GLP-1/GIP RCTs.
- **Outcome evidence.** All-cause mortality HR/SD 1.20 (tier A); CVD HR/SD 1.16 (tier A); Cancer HR/SD 1.14 (tier B); Dementia HR/SD 1.10 (tier C); Frailty HR/SD 1.20 (tier C)

### Waist:hip ratio

*Also: WHR, waist-to-hip ratio*

- **Mechanism.** Waist:hip ratio contrasts central (abdominal) fat against gluteofemoral fat. A high ratio reflects an android, visceral-predominant fat pattern, whereas gluteofemoral (hip) fat is metabolically protective. WHR therefore captures fat distribution rather than total adiposity and predicts cardiometabolic disease and mortality strongly, including in people with normal BMI. Mendelian randomization indicates a causal effect of central fat distribution on cardiometabolic disease.
- **Measurement.** Ratio of tape-measured waist circumference to hip circumference; standardized landmarks required for both measurements. — *sample:* in-office device; *approx. cost:* $0-5.
- **Signal quality.** moderate — Captures fat distribution well and is BMI-independent, but depends on two measurements each with site/observer variability; hip circumference is influenced by both fat and gluteal muscle/bone.
- **Major confounders.** measurement technique for both waist and hip; sex (different distributions and thresholds); ethnicity; menopausal/hormonal status; muscularity of hip/gluteal region
- **Testing cadence.** baseline + annual
- **Modifiability (moderate).** Caloric restriction and aerobic/resistance exercise modestly improve WHR by mobilizing central fat; GLP-1/GIP agonists and bariatric surgery produce larger shifts. WHR is harder to change than waist circumference because hip dimension also shifts. — *effect:* Lifestyle interventions improve WHR modestly (~0.01-0.03 absolute); pharmacologic and surgical weight loss produce larger changes.; *timeframe:* 6-12 months; *evidence:* RCT meta-analysis. Meta-analyses of exercise/diet and body fat distribution; GLP-1/GIP weight-loss RCT programs.
- **Outcome evidence.** All-cause mortality HR/SD 1.15 (tier A); CVD HR/SD 1.27 (tier A); Cancer HR/SD 1.10 (tier C); Dementia HR/SD 1.12 (tier C); Frailty HR/SD 1.18 (tier C)

### BMI

*Also: body mass index, Quetelet index*

- **Mechanism.** Body mass index (weight / height squared) is a simple proxy for overall body size and adiposity. It is cheap, universal and predictive of cardiometabolic disease at the population level, but it does not distinguish fat from lean mass or capture fat distribution. Its relationship with all-cause mortality is U-shaped: both obesity and low BMI raise risk, the latter heavily confounded by smoking, illness and frailty.
- **Measurement.** Calculated from measured weight and height; weight in kg divided by height in metres squared. — *sample:* in-office device; *approx. cost:* $0-5.
- **Signal quality.** moderate — Highly reproducible to measure, but a poor discriminator of body composition: it conflates muscle and fat and misses central adiposity, so it misclassifies risk in muscular, older and sarcopenic individuals.
- **Major confounders.** muscle vs fat mass (not distinguished); fat distribution not captured; smoking (drives low-BMI mortality); reverse causation (illness-related weight loss); age (shifts the optimal range upward); ethnicity (different risk thresholds)
- **Testing cadence.** baseline + annual
- **Standardization (U-shaped).** BMI has a U-shaped 'obesity paradox' relation with all-cause mortality and frailty, and an age-dependent, sign-reversing relation with dementia, so those cells use no per-SD slope; the CVD and cancer cells reflect the monotone upper arm of the curve and are converted from per-5-unit HRs.
- **Modifiability (high).** Caloric restriction and physical activity; GLP-1/GIP receptor agonists and bariatric surgery for substantial reductions in those with obesity. The goal is body-composition improvement, not BMI minimization. — *effect:* Lifestyle programs lower BMI ~1-3 kg/m2; GLP-1/GIP agonists 4-8+ kg/m2; bariatric surgery more.; *timeframe:* 3-12 months; *evidence:* RCT meta-analysis. Diabetes Prevention Program and Look AHEAD lifestyle trials; STEP/SURMOUNT GLP-1/GIP RCTs; bariatric surgery cohorts.
- **Outcome evidence.** All-cause mortality native HR 1.00 [categorical, tier A]; CVD HR/SD 1.19 (tier A); Cancer HR/SD 1.09 (tier A); Dementia native HR 1.00 [categorical, tier B]; Frailty native HR 1.00 [categorical, tier B]

### Bioimpedance phase angle

*Also: phase angle, BIA phase angle, PhA*

- **Mechanism.** Phase angle is derived from bioelectrical impedance analysis as the arctangent of reactance over resistance. It reflects the integrity and mass of cell membranes and the ratio of intracellular to extracellular water. A higher phase angle indicates healthier, better-hydrated cells with intact membranes and more body cell mass; a low phase angle marks cellular deterioration, malnutrition, inflammation and is an emerging prognostic marker of frailty, sarcopenia and mortality.
- **Measurement.** Single- or multi-frequency bioelectrical impedance analysis; phase angle computed from measured resistance and reactance at (typically) 50 kHz. — *sample:* in-office device; *approx. cost:* $20-80.
- **Signal quality.** moderate — A raw, device-measured BIA parameter (no compartment model assumptions), but values depend on device, electrode placement, frequency, hydration, recent food/exercise and ambient conditions; reference ranges are age- and sex-specific.
- **Major confounders.** hydration status; device and frequency differences; electrode placement; recent food intake / exercise; age and sex (strong reference effects); body temperature
- **Testing cadence.** baseline + every 6-12 months
- **Modifiability (moderate).** Resistance exercise and adequate protein/energy intake raise phase angle by increasing body cell mass and cell-membrane integrity; treatment of underlying inflammation and correction of malnutrition also improve it. — *effect:* Resistance-training and nutritional interventions raise phase angle modestly (~0.2-0.5 degrees) over months in older or clinical populations.; *timeframe:* 3-6 months; *evidence:* RCT. RCTs of resistance training and nutritional supplementation reporting improvements in BIA phase angle in older and clinical populations.
- **Outcome evidence.** All-cause mortality HR/SD 1.45 (tier B); CVD HR/SD 1.35 (tier C); Cancer HR/SD 1.50 (tier B); Dementia (documented, no pooled HR); Frailty HR/SD 1.55 (tier B)


## Cardiovascular & Autonomic Function

### Resting blood pressure

*Also: office blood pressure, clinic blood pressure, systolic blood pressure, SBP*

- **Mechanism.** Chronically elevated blood pressure imposes pulsatile and steady-state load on the arterial wall, heart and microvasculature, driving atherosclerosis, left-ventricular hypertrophy, small-vessel cerebral disease and renal injury. It is one of the largest single modifiable contributors to global cardiovascular and all-cause mortality, and mid-life hypertension is a well-established risk factor for later dementia.
- **Measurement.** Seated, rested, validated oscillometric or auscultatory cuff measurement; multiple readings averaged. Systolic BP is the most prognostic single component in middle-aged and older adults. — *sample:* in-office device; *approx. cost:* $0-30.
- **Signal quality.** moderate — High within-person visit-to-visit variability, white-coat and masked effects, cuff/technique error and regression dilution attenuate single-occasion measurements; averaging multiple readings improves reliability.
- **Major confounders.** white-coat / masked hypertension; measurement technique and cuff size; time of day and recent activity/caffeine; antihypertensive medication; acute stress / anxiety; regression dilution from single readings
- **Testing cadence.** baseline + at least annual; more frequent if elevated or treated
- **Standardization (U-shaped (frailty)).** Most BP cells convert cleanly per-SD, but the frailty cell reflects a U-shaped, age-dependent relation (low BP marks frailty in the oldest-old) and is kept as a categorical hypertension contrast.
- **Modifiability (high).** Sodium reduction, weight loss, DASH-style diet, aerobic exercise, alcohol moderation, and antihypertensive medication (the most reliably effective lever). — *effect:* Lifestyle change lowers SBP ~4-11 mmHg; single antihypertensive agents lower SBP ~8-10 mmHg; combination therapy can lower SBP 20+ mmHg.; *timeframe:* Days to weeks for drugs; weeks to months for lifestyle change.; *evidence:* RCT meta-analysis. Blood Pressure Lowering Treatment Trialists' Collaboration meta-analyses; DASH and lifestyle-intervention RCTs.
- **Outcome evidence.** All-cause mortality HR/SD 1.19 (tier A); CVD HR/SD 1.87 (tier A); Cancer HR/SD 1.13 (tier C); Dementia HR/SD 1.49 (tier B); Frailty native HR 1.14 [categorical, tier C]

### Ambulatory blood pressure

*Also: ambulatory blood pressure monitoring, ABPM, 24-hour blood pressure, nighttime blood pressure, nocturnal blood pressure*

- **Mechanism.** Ambulatory monitoring captures blood pressure across daily activity and sleep, removing white-coat effect and capturing nocturnal BP and dipping pattern. Nighttime and 24-hour averages are more strongly and more linearly related to cardiovascular events and mortality than office BP, and reveal masked hypertension and non-dipping that office measurement misses.
- **Measurement.** Oscillometric cuff worn 24 hours, measuring at fixed intervals (e.g., every 15-30 min awake, 30-60 min asleep); daytime, nighttime and 24-hour averages and dipping status are derived. — *sample:* wearable; *approx. cost:* $50-250.
- **Signal quality.** clean — Multiple readings over 24 h sharply reduce measurement error and white-coat/masked misclassification relative to office BP; main limitations are sleep disturbance from the cuff and occasional failed readings.
- **Major confounders.** sleep quality and disruption by the device; activity level during recording; shift work / irregular sleep schedule; antihypertensive medication and dosing time; arm movement artefact
- **Testing cadence.** as indicated for diagnosis/treatment monitoring; not routinely repeated
- **Standardization (Dipping pattern).** Average ambulatory SBP cells convert cleanly per-SD, but the dementia signal rests on a categorical nocturnal-dipping pattern that cannot be placed on a per-SD scale.
- **Modifiability (high).** Same levers as office BP (sodium reduction, weight loss, exercise, antihypertensives); evening dosing of antihypertensives can specifically lower nighttime BP and improve dipping in some patients. — *effect:* Antihypertensive therapy lowers 24-hour SBP ~10-15 mmHg; lifestyle change ~5-10 mmHg over 24-hour average.; *timeframe:* Weeks for drug effects; weeks to months for lifestyle change.; *evidence:* RCT meta-analysis. Ambulatory-BP substudies of antihypertensive RCTs; BPLTTC meta-analyses.
- **Outcome evidence.** All-cause mortality HR/SD 1.29 (tier A); CVD HR/SD 1.30 (tier A); Cancer (documented, no pooled HR); Dementia native HR 1.31 [categorical, tier C]; Frailty (documented, no pooled HR)

### Pulse wave velocity (arterial stiffness)

*Also: carotid-femoral pulse wave velocity, cfPWV, aortic stiffness, arterial stiffness*

- **Mechanism.** Carotid-femoral pulse wave velocity is the reference measure of aortic stiffness. A stiffer aorta transmits pulsatile pressure into the cerebral and renal microcirculation, raises central pulse pressure and afterload, and reflects cumulative arterial aging. It is a strong, partly independent predictor of cardiovascular events, and aortic stiffness contributes to small-vessel brain damage and cognitive decline.
- **Measurement.** Transit time of the pressure wave between carotid and femoral sites divided by path distance, using applanation tonometry or oscillometric/Doppler devices; reported in m/s. — *sample:* imaging; *approx. cost:* $100-300.
- **Signal quality.** moderate — Reproducible with trained operators and standardized distance estimation, but device type, path-length method, heart rate and blood pressure at the time of measurement all influence the value; intra-individual variability is moderate.
- **Major confounders.** blood pressure at time of measurement; heart rate; device and path-length estimation method; age (strong); operator technique
- **Testing cadence.** baseline + every few years; not routinely repeated
- **Standardization (Mixed metrics).** CVD evidence is native per-SD and mortality is per-unit-convertible, but dementia/frailty cells rest on heterogeneous categorical high-vs-low contrasts from small cohorts.
- **Modifiability (moderate).** Aerobic exercise, blood-pressure lowering, weight loss, sodium restriction; aerobic training has the most consistent direct destiffening effect. — *effect:* Aerobic exercise lowers cfPWV ~0.5-1.0 m/s in middle-aged/older adults; antihypertensive therapy reduces PWV partly via BP lowering. Stiffness rises ~0.1 m/s per year with aging, so gains are partial.; *timeframe:* 8-12 weeks for exercise effects.; *evidence:* RCT meta-analysis. Exercise-training RCT meta-analyses on arterial stiffness (e.g., aerobic-training PWV meta-analyses).
- **Outcome evidence.** All-cause mortality HR/SD 1.52 (tier A); CVD HR/SD 1.30 (tier A); Cancer (documented, no pooled HR); Dementia native HR 1.50 [categorical, tier C]; Frailty native HR 1.45 [categorical, tier C]

### Endothelial function (FMD / peripheral arterial tonometry)

*Also: flow-mediated dilation, FMD, brachial artery reactivity, peripheral arterial tonometry, reactive hyperemia index, RHI, EndoPAT*

- **Mechanism.** Endothelial function reflects the capacity of the vascular endothelium to release nitric oxide and dilate in response to increased shear stress. Endothelial dysfunction is an early, potentially reversible step in atherosclerosis, preceding structural disease, and integrates the net effect of cardiovascular risk factors on the vascular wall.
- **Measurement.** Brachial-artery flow-mediated dilation (ultrasound, % diameter change after cuff-occlusion hyperemia) or peripheral arterial tonometry (reactive hyperemia index, e.g., EndoPAT). Both require standardized protocols. — *sample:* functional test; *approx. cost:* $100-400.
- **Signal quality.** noisy — FMD is highly operator-dependent with substantial measurement variability; results depend on baseline diameter, cuff position, temperature, fasting and prior exercise/diet. Tonometry indices are more reproducible but measure a related, not identical, construct.
- **Major confounders.** operator skill and image analysis; baseline arterial diameter; room temperature, fasting state, recent food/exercise/caffeine; menstrual cycle phase; time of day; vasoactive medications
- **Testing cadence.** research / specialized use; not part of routine screening
- **Standardization (Noisy / patient pop.).** FMD has high measurement noise and most outcome evidence comes from patient populations; CVD/mortality cells are per-unit-convertible but the per-SD exponent is large, and dementia evidence is categorical and weak.
- **Modifiability (moderate).** Aerobic exercise, smoking cessation, weight loss, Mediterranean-style diet; statins and antihypertensives also improve endothelial function. — *effect:* Aerobic exercise improves brachial FMD by ~1-3 percentage points in adults with impaired function; effects are larger in those with greater baseline dysfunction.; *timeframe:* 4-12 weeks; some dietary effects within days.; *evidence:* RCT meta-analysis. Exercise-training RCT meta-analyses on flow-mediated dilation.
- **Outcome evidence.** All-cause mortality HR/SD 0.59 (tier B); CVD HR/SD 0.66 (tier B); Cancer (documented, no pooled HR); Dementia native HR 0.79 [categorical, tier C]; Frailty (documented, no pooled HR)

### Resting heart rate

*Also: resting pulse, RHR, heart rate at rest*

- **Mechanism.** Resting heart rate reflects the balance of sympathetic and parasympathetic autonomic tone, cardiorespiratory fitness and overall cardiovascular health. Elevated resting heart rate is associated with greater myocardial oxygen demand, shorter diastolic perfusion time, higher arterial shear stress, and is a marker of reduced vagal tone; it consistently predicts cardiovascular and all-cause mortality.
- **Measurement.** Pulse or ECG/heart-rate measurement after seated rest; widely available from clinical visits and consumer wearables, the latter enabling repeated overnight/resting measurement. — *sample:* wearable; *approx. cost:* $0-30.
- **Signal quality.** moderate — Easy and cheap to measure but sensitive to recent activity, caffeine, stress, posture, time of day, fever and rate-affecting drugs; averaged resting or overnight values from wearables are more stable than a single clinic reading.
- **Major confounders.** physical activity / fitness; caffeine, nicotine, stress, anxiety; beta-blockers and other rate-affecting drugs; fever / acute illness; thyroid status; posture and time of measurement
- **Testing cadence.** baseline + annual; continuous/frequent if wearable-tracked
- **Standardization (Marker not causal).** Most cells convert cleanly per-SD, but resting heart rate is largely a marker of fitness/autonomic tone (MR and rate-lowering trials give limited causal support); the dementia cell uses categorical heart-rate bands.
- **Modifiability (high).** Regular aerobic (endurance) exercise is the primary lever; weight loss, reduced stimulant intake, stress reduction and treatment of conditions raising heart rate also help. Beta-blockers lower heart rate pharmacologically. — *effect:* Endurance training lowers resting heart rate by ~5-12 bpm depending on baseline fitness and training volume.; *timeframe:* Several weeks to a few months of consistent training.; *evidence:* RCT meta-analysis. Endurance-training RCT meta-analyses on resting heart rate (e.g., aerobic-exercise resting-HR meta-analyses).
- **Outcome evidence.** All-cause mortality HR/SD 1.10 (tier A); CVD HR/SD 1.13 (tier A); Cancer HR/SD 1.16 (tier B); Dementia native HR 1.55 [categorical, tier B]; Frailty HR/SD 1.07 (tier C)

### Echocardiographic measures (LV mass / diastolic function)

*Also: left ventricular mass, LV mass index, LV hypertrophy, diastolic function, diastolic dysfunction, E/e' ratio, left atrial volume*

- **Mechanism.** Echocardiographic structural and functional measures - left ventricular mass, geometry and diastolic function (relaxation and filling pressures) - reflect the cumulative effect of pressure load, ischemia, metabolic stress and aging on the heart. Increased LV mass and impaired diastolic function are subclinical markers of cardiac end-organ damage and precede heart failure (especially HFpEF), atrial fibrillation, stroke and death.
- **Measurement.** Transthoracic echocardiography: LV mass from linear or 3D measurements (indexed to body size); diastolic function from mitral inflow, tissue-Doppler e' and E/e', and left atrial volume index. Cardiac MRI is the more precise reference for mass. — *sample:* imaging; *approx. cost:* $200-1000.
- **Signal quality.** moderate — Echocardiographic LV mass and diastolic indices are reproducible with experienced operators and standardized protocols, but depend on image quality, geometric assumptions, loading conditions and inter-observer variability; diastolic grading is partly categorical.
- **Major confounders.** blood pressure / loading conditions at time of scan; body size and indexation method; image quality and acoustic windows; inter-observer variability; age (e' declines with age); atrial fibrillation (limits diastolic assessment)
- **Testing cadence.** as clinically indicated; not part of routine population screening
- **Standardization (Mixed; thin late-life).** LV mass cells are native per-SD, but diastolic-function grading is partly categorical and dementia/frailty evidence rests on present-vs-absent contrasts from heterogeneous cohorts.
- **Modifiability (moderate).** Blood-pressure control is the dominant lever for LV mass regression; weight loss, sodium reduction, treatment of sleep apnea and aerobic exercise also help. Diastolic function is less readily reversible. — *effect:* Effective antihypertensive therapy regresses LV mass by ~10-15% over ~6-12 months (ARB/ACE-inhibitor and CCB regimens most effective); diastolic improvement is more modest.; *timeframe:* 6-12 months for meaningful LV mass regression.; *evidence:* RCT meta-analysis. Antihypertensive RCT meta-analyses of LV mass regression (e.g., Klingbeil 2003 meta-analysis of drug-class effects on LV mass).
- **Outcome evidence.** All-cause mortality HR/SD 1.27 (tier B); CVD HR/SD 1.40 (tier B); Cancer (documented, no pooled HR); Dementia native HR 1.35 [categorical, tier C]; Frailty (documented, no pooled HR)


## Renal & Hepatic Function

### eGFR (cystatin C-based)

*Also: estimated glomerular filtration rate (cystatin C), eGFRcys, CKD-EPI cystatin C eGFR*

- **Mechanism.** Glomerular filtration rate estimated from serum cystatin C indexes kidney clearance function. Unlike creatinine, cystatin C is largely independent of muscle mass, so cystatin C-based eGFR detects kidney impairment earlier and more accurately predicts mortality and cardiovascular events; low eGFR reflects nephron loss, vascular aging, and systemic risk burden.
- **Measurement.** Serum cystatin C immunoturbidimetric/nephelometric assay (IFCC-standardized) entered into the CKD-EPI cystatin C equation; reported as mL/min/1.73m2. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** clean — Cystatin C is less affected by muscle mass and diet than creatinine; assays are internationally standardized. Modestly influenced by thyroid status, corticosteroids, and inflammation.
- **Major confounders.** thyroid dysfunction; corticosteroid / immunosuppressant use; systemic inflammation and adiposity; smoking; acute illness
- **Testing cadence.** baseline + annual; more frequent in CKD or with nephrotoxic exposures
- **Modifiability (moderate).** Blood-pressure control (RAS blockade), glycemic control, SGLT2 inhibitors, avoidance of nephrotoxins, and treatment of the underlying disease slow eGFR decline; SGLT2 inhibitors and finerenone reduce eGFR slope decline in trials. — *effect:* SGLT2 inhibitors slow chronic eGFR decline by ~1-2 mL/min/1.73m2 per year and reduce kidney-failure and mortality endpoints; established eGFR loss is largely not reversible.; *timeframe:* Months to years (slope effect); an acute dip then stabilization within weeks of starting SGLT2i/RAS blockade; *evidence:* RCT meta-analysis. Nuffield Department of Population Health SGLT2 inhibitor Meta-Analysis (SMART-C). Lancet 2022;400:1788-1801. doi:10.1016/S0140-6736(22)02074-8
- **Outcome evidence.** All-cause mortality HR/SD 1.45 (tier A); CVD HR/SD 1.39 (tier A); Cancer HR/SD 1.10 (tier C); Dementia HR/SD 1.16 (tier B); Frailty HR/SD 1.30 (tier C)

### Cystatin C

*Also: serum cystatin C, CST3 protein*

- **Mechanism.** Cystatin C is a low-molecular-weight protein produced by all nucleated cells at a near-constant rate, freely filtered by the glomerulus and catabolized in the tubule, so its serum concentration is an inverse marker of glomerular filtration. It also reflects inflammation, adiposity, and vascular aging, making it a strong, partly non-renal predictor of mortality and aging outcomes.
- **Measurement.** Immunoturbidimetric or immunonephelometric assay on a standard chemistry analyzer; IFCC-standardized reference material. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** clean — Low analytic variability and internationally standardized; independent of muscle mass. Elevated by corticosteroids, hyperthyroidism, inflammation and high adiposity, which adds non-renal signal.
- **Major confounders.** corticosteroid / immunosuppressant use; thyroid dysfunction; systemic inflammation; obesity; smoking
- **Testing cadence.** baseline + annual
- **Modifiability (low).** No direct cystatin C-targeting therapy; the marker improves with control of kidney injury (blood pressure, glycemic control, SGLT2 inhibitors), weight loss, and reduction of corticosteroid exposure and inflammation. — *effect:* Cystatin C tracks the underlying GFR and inflammatory state; meaningful change requires improving kidney function or body composition. Modest reductions seen with weight loss and SGLT2 inhibitors.; *timeframe:* Months; *evidence:* observational. Weight-loss and SGLT2 inhibitor cohort/trial data showing modest cystatin C reduction; no direct cystatin C-lowering therapy exists.
- **Outcome evidence.** All-cause mortality HR/SD 1.30 (tier A); CVD HR/SD 1.32 (tier A); Cancer HR/SD 1.12 (tier C); Dementia HR/SD 1.20 (tier C); Frailty HR/SD 1.28 (tier C)

### Albuminuria (urine albumin:creatinine ratio)

*Also: UACR, urine albumin-to-creatinine ratio, microalbuminuria, ACR*

- **Mechanism.** Urinary albumin excretion reflects glomerular barrier dysfunction and is a sensitive marker of generalized endothelial and microvascular damage, not just kidney disease. Even high-normal albuminuria predicts cardiovascular events, kidney failure, and death, making UACR a window onto systemic vascular health.
- **Measurement.** Albumin and creatinine measured on a single spot urine sample (preferably first-morning); ratio reported as mg/g or mg/mmol. Confirm with repeat testing given day-to-day variability. — *sample:* urine; *approx. cost:* $10-30.
- **Signal quality.** moderate — High within-person biological variability; affected by exercise, fever, posture, hydration and menstrual contamination. Repeat measurement and first-morning samples improve reliability; the marker is strongly right-skewed and analyzed on a log scale.
- **Major confounders.** recent vigorous exercise; fever / acute illness / urinary tract infection; uncontrolled hypertension or hyperglycemia; menstrual blood contamination; hydration status and upright posture
- **Testing cadence.** baseline + annual; more frequent in diabetes, hypertension, or CKD
- **Modifiability (high).** RAS blockade (ACE inhibitors / ARBs), SGLT2 inhibitors, the nonsteroidal MRA finerenone, GLP-1 receptor agonists, blood-pressure and glycemic control, weight loss and sodium restriction reduce albuminuria. — *effect:* RAS blockade lowers UACR ~30-40%; SGLT2 inhibitors and finerenone add further reductions; substantial albuminuria reduction is achievable and tracks lower kidney and CV risk.; *timeframe:* Weeks to a few months; *evidence:* RCT meta-analysis. FIDELIO-DKD/FIGARO-DKD (finerenone) and CREDENCE/DAPA-CKD (SGLT2 inhibitor) trials; RAS blockade albuminuria meta-analyses.
- **Outcome evidence.** All-cause mortality HR/SD 1.50 (tier A); CVD HR/SD 1.59 (tier A); Cancer HR/SD 1.10 (tier C); Dementia HR/SD 1.22 (tier B); Frailty HR/SD 1.25 (tier C)

### ALT

*Also: alanine aminotransferase, SGPT, serum glutamic-pyruvic transaminase*

- **Mechanism.** Alanine aminotransferase is a hepatocyte enzyme released into blood on liver-cell injury; elevated ALT indicates hepatocellular damage, most commonly from metabolic dysfunction-associated steatotic liver disease (MASLD), and tracks insulin resistance and visceral adiposity. Paradoxically, very low ALT in older adults reflects depleted muscle mass and frailty.
- **Measurement.** Enzymatic spectrophotometric assay on a standard chemistry analyzer; part of routine liver function panels. — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** moderate — Widely available and cheap, but reference ranges are population- and sex-dependent and many labs use upper limits that are too high; values vary with adiposity, muscle mass, vigorous exercise and time of day.
- **Major confounders.** muscle mass / sarcopenia (low ALT); vigorous exercise and muscle injury; obesity and insulin resistance; alcohol intake; hepatotoxic medications; sex and age
- **Testing cadence.** baseline + annual; more frequent with known liver disease or hepatotoxic drugs
- **Modifiability (moderate).** Weight loss, exercise, reduced alcohol and treatment of MASLD lower elevated ALT; for low ALT in older adults, resistance training and adequate protein to rebuild muscle mass are the relevant interventions. — *effect:* 7-10% weight loss can normalize elevated ALT in MASLD; resistance training raises muscle mass and can modestly raise low ALT.; *timeframe:* 3-6 months; *evidence:* RCT meta-analysis. Lifestyle/weight-loss RCTs in MASLD (e.g. Vilar-Gomez 2015); resistance-training trials for sarcopenia.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.05 (tier B); Cancer HR/SD 1.08 (tier C); Dementia HR/SD 0.92 (tier C); Frailty HR/SD 0.80 (tier B)

### AST

*Also: aspartate aminotransferase, SGOT, serum glutamic-oxaloacetic transaminase*

- **Mechanism.** Aspartate aminotransferase is released on injury to hepatocytes but also skeletal/cardiac muscle and red cells, making it less liver-specific than ALT. The AST/ALT ratio (De Ritis ratio) rises with advanced fibrosis, alcoholic liver disease, and sarcopenia, and is itself a mortality marker.
- **Measurement.** Enzymatic spectrophotometric assay on a standard chemistry analyzer; part of routine liver function panels. — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** moderate — Less liver-specific than ALT; elevated by muscle injury, hemolysis, vigorous exercise and cardiac events. Reference ranges are population-dependent.
- **Major confounders.** skeletal / cardiac muscle injury; vigorous exercise; hemolysis of the sample; alcohol intake; hepatotoxic medications; muscle mass
- **Testing cadence.** baseline + annual; more frequent with liver disease or hepatotoxic drugs
- **Modifiability (moderate).** Weight loss, exercise, alcohol reduction and treatment of underlying liver disease lower elevated AST; for low AST in older adults, resistance training and adequate protein are relevant. — *effect:* Lifestyle change normalizes mildly elevated AST in MASLD; muscle-building raises low AST modestly.; *timeframe:* 3-6 months; *evidence:* RCT meta-analysis. Lifestyle/weight-loss RCTs in MASLD; resistance-training trials for sarcopenia.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.06 (tier C); Cancer HR/SD 1.10 (tier C); Dementia HR/SD 0.93 (tier C); Frailty HR/SD 0.85 (tier C)

### GGT

*Also: gamma-glutamyl transferase, gamma-GT, GGTP*

- **Mechanism.** Gamma-glutamyl transferase is a cell-surface enzyme central to glutathione metabolism and is induced by oxidative stress, alcohol, and hepatic fat. Circulating GGT integrates liver injury, oxidative stress, and metabolic dysfunction, making it one of the stronger liver-enzyme predictors of mortality and cardiovascular disease.
- **Measurement.** Enzymatic spectrophotometric assay on a standard chemistry analyzer; part of routine liver panels. — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** moderate — Strongly right-skewed and analyzed on a log scale; sensitive but non-specific - induced by alcohol, enzyme-inducing drugs, obesity and biliary disease. Reproducible analytically.
- **Major confounders.** alcohol intake; enzyme-inducing medications (e.g. anticonvulsants); obesity and insulin resistance; biliary obstruction; smoking
- **Testing cadence.** baseline + annual
- **Modifiability (high).** Alcohol reduction or abstinence, weight loss, exercise, treatment of MASLD and withdrawal of enzyme-inducing drugs lower GGT. — *effect:* Alcohol cessation can roughly halve an elevated GGT within weeks; weight loss and exercise produce further reductions.; *timeframe:* 2-8 weeks for alcohol-related GGT; 3-6 months for metabolic causes; *evidence:* observational. Alcohol-cessation and lifestyle intervention cohort studies showing rapid GGT decline.
- **Outcome evidence.** All-cause mortality HR/SD 1.22 (tier A); CVD HR/SD 1.20 (tier A); Cancer HR/SD 1.18 (tier B); Dementia HR/SD 1.12 (tier C); Frailty HR/SD 1.15 (tier C)

### Fibrosis-4 index (FIB-4)

*Also: FIB-4 index, Fibrosis-4 score*

- **Mechanism.** FIB-4 is a composite index of liver fibrosis calculated from age, AST, ALT, and platelet count. It non-invasively estimates the degree of hepatic fibrosis - the structural change that drives liver-related mortality in MASLD and other chronic liver disease - and higher FIB-4 also flags advanced metabolic disease and predicts all-cause and cardiovascular death.
- **Measurement.** Calculated: (age x AST) / (platelet count x sqrt(ALT)); requires only routine blood tests. Validated cut-points ~1.3 (rule out) and ~2.67 (rule in advanced fibrosis); age-adjusted cut-points used over 65. — *sample:* venous blood; *approx. cost:* $0-15 (derived from routine labs).
- **Signal quality.** moderate — Derived index; inherits AST/ALT and platelet variability and loses specificity at older ages (the unadjusted low cut-point over-flags people over 65). Strongly right-skewed; best used with validated cut-points.
- **Major confounders.** age (inflates the score; needs age-specific cut-points); AST/ALT variability (muscle injury, exercise); thrombocytopenia from non-liver causes; acute illness; alcohol intake
- **Testing cadence.** baseline + every 1-3 years in metabolic-risk patients; more frequent if elevated
- **Modifiability (moderate).** Weight loss, exercise, alcohol reduction, treatment of MASLD/MASH (including resmetirom and GLP-1 receptor agonists) and control of metabolic risk can stabilize or reduce FIB-4 by improving the underlying fibrosis and its AST/ALT inputs. — *effect:* Significant weight loss and MASH-directed therapy can lower FIB-4 and reverse early fibrosis; established advanced fibrosis (cirrhosis) is largely fixed.; *timeframe:* 6-24 months; *evidence:* RCT meta-analysis. MASH RCTs (e.g. MAESTRO-NASH resmetirom, GLP-1 RA trials) and lifestyle/weight-loss studies showing fibrosis-marker improvement.
- **Outcome evidence.** All-cause mortality HR/SD 1.45 (tier B); CVD HR/SD 1.28 (tier B); Cancer HR/SD 1.35 (tier C); Dementia HR/SD 1.18 (tier C); Frailty HR/SD 1.30 (tier C)

### Uric acid

*Also: serum urate, serum uric acid, SUA*

- **Mechanism.** Uric acid is the end-product of purine metabolism, cleared mainly by the kidney. High levels cause gout and reflect renal impairment, insulin resistance, and metabolic syndrome, and may promote oxidative stress and endothelial dysfunction; yet urate is also a major plasma antioxidant, so very low levels may be detrimental, producing a U-shaped relationship with mortality.
- **Measurement.** Enzymatic (uricase) colorimetric assay on a standard chemistry analyzer; part of routine metabolic panels. — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** moderate — Analytically reliable, but strongly influenced by renal clearance, diuretics, diet, alcohol and sex; the U-shaped outcome relationship complicates interpretation of a single value.
- **Major confounders.** kidney function / eGFR; diuretic and urate-lowering medication; alcohol and purine-rich diet; sex and menopausal status; obesity and insulin resistance
- **Testing cadence.** baseline + annual; more frequent during urate-lowering therapy
- **Modifiability (high).** Urate-lowering drugs (allopurinol, febuxostat), weight loss, reduced alcohol and purine intake, and avoidance of urate-raising diuretics lower serum uric acid; SGLT2 inhibitors and losartan modestly lower urate. — *effect:* Allopurinol/febuxostat reliably lower urate to target (<6 mg/dL); lifestyle change lowers it modestly. Whether lowering urate improves non-gout outcomes is unproven (neutral cardiovascular trials).; *timeframe:* 2-8 weeks for pharmacologic lowering; *evidence:* RCT meta-analysis. Urate-lowering therapy RCTs (e.g. CARES, FAST for safety; ALL-HEART for cardiovascular outcomes - neutral).
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.13 (tier A); Cancer HR/SD 1.05 (tier C); Dementia HR/SD 0.94 (tier B); Frailty (documented, no pooled HR)


## Cognitive & Neurological

### Processing speed

*Also: cognitive processing speed, psychomotor speed, Digit Symbol Substitution Test, Symbol Digit Modalities Test, perceptual speed, reaction time (cognitive)*

- **Mechanism.** Processing speed is the rate at which simple cognitive operations are executed; it is the cognitive domain most sensitive to aging, cerebral small-vessel disease and white-matter integrity. Slowing reflects diffuse neural and vascular compromise and is an early, sensitive marker of incipient neurodegeneration and global physiologic decline, predicting dementia, disability and death.
- **Measurement.** Timed paper-and-pencil or computerized tasks: Digit Symbol Substitution Test (DSST), Symbol Digit Modalities Test (SDMT), Trail Making Test A, simple/choice reaction time, and NIH Toolbox Pattern Comparison Processing Speed Test. Cambridge Brain Sciences includes speeded reasoning tasks. Scores are typically standardized to within-cohort z-scores. — *sample:* functional test; *approx. cost:* $0-75.
- **Signal quality.** moderate — Speeded tasks are reliable (test-retest r ~0.7-0.85) but sensitive to motor speed, vision, fatigue, practice effects and education; cross-instrument comparability is imperfect, so cells are standardized within study.
- **Major confounders.** age; education / cognitive reserve; depression; vision and motor impairment; medication (sedatives, anticholinergics); practice effects; subclinical cerebrovascular disease (reverse causation)
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (moderate).** Aerobic and combined exercise training, speed-of-processing cognitive training (e.g. ACTIVE-trial Useful Field of View training), and management of vascular risk factors (blood pressure, glycemia). — *effect:* ACTIVE speed-of-processing training produced large immediate gains (effect size ~1.5 SD on trained tasks) with partial 10-year persistence; aerobic exercise yields small-to-moderate gains (~0.2-0.3 SD). Far transfer to function is modest.; *timeframe:* Weeks for trained-task gains; 6-12 months for exercise/vascular effects.; *evidence:* RCT. Rebok GW, et al. Ten-year effects of the ACTIVE cognitive training trial. J Am Geriatr Soc. 2014;62(1):16-24. doi:10.1111/jgs.12607. PMID:24417410
- **Outcome evidence.** All-cause mortality HR/SD 1.27 (tier A); CVD HR/SD 1.20 (tier B); Cancer HR/SD 1.05 (tier C); Dementia HR/SD 1.65 (tier A); Frailty HR/SD 1.35 (tier B)

### Working memory

*Also: working memory capacity, n-back performance, digit span, list sorting, NIH Toolbox List Sorting, Cambridge Brain Sciences Monkey Ladder*

- **Mechanism.** Working memory is the capacity to hold and manipulate information over short intervals; it depends on frontoparietal networks and is sensitive to aging, dopaminergic tone and early Alzheimer pathology. Decline contributes to loss of independence and predicts dementia, though it is a weaker mortality predictor than processing speed.
- **Measurement.** Digit span backward, spatial span, n-back tasks, NIH Toolbox List Sorting Working Memory Test, and Cambridge Brain Sciences spatial-span/Monkey-Ladder tasks. Scored as within-cohort standardized z-scores. — *sample:* functional test; *approx. cost:* $0-75.
- **Signal quality.** moderate — Span tasks are moderately reliable (test-retest r ~0.7-0.8); n-back reliability is lower. Sensitive to attention, motivation, anxiety and practice effects; cross-instrument standardization is required.
- **Major confounders.** age; education / cognitive reserve; depression and anxiety; attention/effort; practice effects; subclinical neurodegeneration (reverse causation)
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (low).** Working-memory and process-based cognitive training, aerobic and combined exercise, and vascular-risk management. — *effect:* Working-memory training produces reliable near-transfer gains (~0.5-0.7 SD on trained tasks) but limited far transfer to untrained abilities; exercise yields small gains (~0.1-0.2 SD).; *timeframe:* Weeks for trained-task gains; 3-6 months for exercise effects.; *evidence:* RCT meta-analysis. Melby-Lervag M, Hulme C. Is working memory training effective? A meta-analytic review. Dev Psychol. 2013;49(2):270-291. doi:10.1037/a0028228. PMID:22612437 (near transfer reliable, far transfer minimal).
- **Outcome evidence.** All-cause mortality HR/SD 1.15 (tier B); CVD HR/SD 1.12 (tier C); Cancer (documented, no pooled HR); Dementia HR/SD 1.55 (tier B); Frailty HR/SD 1.25 (tier C)

### Executive function

*Also: executive control, cognitive flexibility, Trail Making Test B, Stroop, verbal fluency, NIH Toolbox Flanker / DCCS, Cambridge Brain Sciences reasoning*

- **Mechanism.** Executive function comprises inhibition, set-shifting and goal management, dependent on prefrontal-subcortical circuits. It governs self-care, medication adherence and health decision-making, and its decline reflects frontal-lobe and small-vessel pathology; it strongly predicts dementia, disability and mortality.
- **Measurement.** Trail Making Test B (and B-A), Stroop interference, verbal/category fluency, Wisconsin Card Sorting, and NIH Toolbox Flanker Inhibitory Control and Dimensional Change Card Sort; Cambridge Brain Sciences includes planning/reasoning tasks. Scored as within-cohort standardized composite z-scores. — *sample:* functional test; *approx. cost:* $0-75.
- **Signal quality.** moderate — Composite executive scores are moderately reliable (r ~0.7-0.85); individual tasks are less reliable and 'task-impure' (load on processing speed and working memory). Composite standardization mitigates instrument heterogeneity.
- **Major confounders.** age; education / cognitive reserve; processing-speed contamination of timed tasks; depression; practice effects; subclinical cerebrovascular and neurodegenerative disease (reverse causation)
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (moderate).** Aerobic and combined exercise training, multidomain lifestyle interventions (e.g. FINGER trial), and vascular-risk-factor management. — *effect:* Multidomain intervention (FINGER) improved executive-function composite by ~0.15-0.25 SD versus control over 2 years; aerobic exercise yields small gains (~0.1-0.2 SD).; *timeframe:* 6-24 months for measurable change.; *evidence:* RCT. Ngandu T, et al. A 2 year multidomain intervention of diet, exercise, cognitive training, and vascular risk monitoring versus control (FINGER): a randomised controlled trial. Lancet. 2015;385(9984):2255-2263. doi:10.1016/S0140-6736(15)60461-5. PMID:25771249
- **Outcome evidence.** All-cause mortality HR/SD 1.22 (tier A); CVD HR/SD 1.18 (tier B); Cancer HR/SD 1.04 (tier C); Dementia HR/SD 1.60 (tier A); Frailty HR/SD 1.40 (tier B)

### pTau-217

*Also: plasma phosphorylated tau 217, p-tau217, phospho-tau 217*

- **Mechanism.** Plasma phosphorylated tau at threonine-217 rises early in the Alzheimer pathophysiological cascade, tracking cerebral amyloid-beta deposition before tau-tangle accumulation and symptoms. It is the most accurate blood marker of Alzheimer pathology and predicts progression from preclinical/MCI states to dementia.
- **Measurement.** Immunoassay (e.g. Lilly/MSD, ALZpath, Janssen) or mass spectrometry on plasma; results are concentration-dependent and assay-specific, with %p-tau217 ratio variants. Increasingly available as a clinical lab test. — *sample:* venous blood; *approx. cost:* $200-700.
- **Signal quality.** clean — Excellent analytical performance and large effect sizes for Alzheimer pathology (AUC ~0.9-0.95 versus amyloid PET); main caveats are assay-to-assay non-comparability and renal-function/BMI effects on absolute levels.
- **Major confounders.** chronic kidney disease (raises levels); BMI (dilutional, lowers levels); assay platform differences; age; APOE genotype; acute illness
- **Testing cadence.** baseline; repeat at 1-2 years if monitoring
- **Modifiability (low).** Anti-amyloid monoclonal antibodies (lecanemab, donanemab) lower plasma pTau-217 as a treatment-response biomarker; no lifestyle intervention has an established effect on pTau-217. — *effect:* Anti-amyloid therapy reduces plasma pTau-217 by roughly 20-40% over 12-18 months, paralleling amyloid-PET clearance; lifestyle effects unproven.; *timeframe:* 6-18 months with anti-amyloid therapy.; *evidence:* RCT. Donanemab (TRAILBLAZER-ALZ 2; Sims JR, et al. JAMA. 2023;330(6):512-527) and lecanemab (Clarity AD; van Dyck CH, et al. N Engl J Med. 2023;388:9-21) trials report plasma pTau-217 reductions with treatment.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD (documented, no pooled HR); Cancer (documented, no pooled HR); Dementia HR/SD 1.76 (tier B); Frailty (documented, no pooled HR)

### Neurofilament light chain (NfL)

*Also: neurofilament light, NEFL, plasma NfL, serum NfL*

- **Mechanism.** Neurofilament light chain is a structural axonal protein released into blood upon neuro-axonal injury. Plasma NfL is a non-specific but highly sensitive marker of neurodegeneration that rises with age and across diverse neurological diseases, and elevated levels predict cognitive decline, dementia and mortality.
- **Measurement.** Ultrasensitive single-molecule immunoassay (Simoa) on serum or plasma; widely used in research and increasingly clinically available. Levels rise steeply with age, so age-adjusted reference ranges are needed. — *sample:* venous blood; *approx. cost:* $100-400.
- **Signal quality.** moderate — Excellent analytical sensitivity and reproducibility, but biologically non-specific: elevated by any neuro-axonal injury, strongly age- and renal-function-dependent. Disease specificity is low, so it is best read as a general neurodegeneration index.
- **Major confounders.** age (strong); chronic kidney disease (raises levels); BMI / blood volume (dilutional); any concurrent neurological injury or comorbidity; cardiovascular disease
- **Testing cadence.** baseline; repeat at 1-2 years if monitoring
- **Modifiability (low).** No established intervention reliably lowers plasma NfL in non-disease states; disease-modifying therapies (e.g. in MS, SMA) reduce NfL, and vascular-risk control may slow NfL rise. Lifestyle effects are unproven. — *effect:* Disease-modifying drugs in active neurological disease lower NfL substantially (30-70%); no demonstrated lifestyle-driven reduction in healthy aging.; *timeframe:* Months with effective disease-modifying therapy; no established lifestyle timeframe.; *evidence:* RCT. NfL reductions with disease-modifying therapy are established in multiple sclerosis and spinal muscular atrophy trials; no RCT shows lifestyle modification lowers NfL in healthy older adults.
- **Outcome evidence.** All-cause mortality HR/SD 1.40 (tier B); CVD HR/SD 1.30 (tier C); Cancer (documented, no pooled HR); Dementia HR/SD 1.50 (tier B); Frailty HR/SD 1.35 (tier C)

### GFAP

*Also: glial fibrillary acidic protein, plasma GFAP, serum GFAP*

- **Mechanism.** Glial fibrillary acidic protein is an astrocytic intermediate-filament protein released into blood during reactive astrogliosis. Plasma GFAP rises early in Alzheimer disease (tracking amyloid pathology and astrocyte activation) and predicts cognitive decline and incident dementia; it is more Alzheimer-associated than NfL.
- **Measurement.** Ultrasensitive single-molecule immunoassay (Simoa) on serum or plasma; research-grade and increasingly clinically available. Levels rise with age and are affected by renal function and BMI. — *sample:* venous blood; *approx. cost:* $100-400.
- **Signal quality.** moderate — Good analytical performance; elevated by astrogliosis from Alzheimer pathology, traumatic brain injury and other CNS insults, so it is moderately but not fully Alzheimer-specific. Age, kidney function and BMI confound absolute levels.
- **Major confounders.** age; chronic kidney disease; BMI (dilutional); traumatic brain injury and other acute CNS insults; assay platform differences
- **Testing cadence.** baseline; repeat at 1-2 years if monitoring
- **Modifiability (low).** Anti-amyloid monoclonal antibodies lower plasma GFAP modestly as a downstream response biomarker; no lifestyle intervention has an established effect on GFAP. — *effect:* Anti-amyloid therapy produces modest GFAP reductions over 12-18 months; lifestyle effects unproven.; *timeframe:* 12-18 months with anti-amyloid therapy.; *evidence:* RCT. Anti-amyloid trials (donanemab, lecanemab) report plasma GFAP reductions as exploratory biomarkers; no RCT shows lifestyle modification lowers GFAP in healthy aging.
- **Outcome evidence.** All-cause mortality HR/SD 1.35 (tier C); CVD HR/SD 1.20 (tier C); Cancer (documented, no pooled HR); Dementia HR/SD 1.55 (tier B); Frailty (documented, no pooled HR)


## Sleep & Recovery

### Total sleep time

*Also: sleep duration, habitual sleep duration, TST*

- **Mechanism.** Sleep duration is the total time spent asleep per 24h. Sleep supports glymphatic clearance, metabolic and immune regulation, blood-pressure dipping, and neural restoration. Both insufficient and excessively long sleep are associated with adverse cardiometabolic, cognitive and mortality outcomes; long sleep is partly a marker of underlying ill-health, fragmentation or comorbidity rather than a cause.
- **Measurement.** Self-report (questionnaire/diary) in most large cohorts; objectively by actigraphy, wearable accelerometry or polysomnography. Self-report and objective estimates correlate only moderately. — *sample:* wearable; *approx. cost:* $0-150.
- **Signal quality.** moderate — Self-reported duration is prone to recall and rounding bias and over-estimates objective sleep; objective measures vary by device and night. Single-night and single-question estimates are noisy; multi-night averages are more stable.
- **Major confounders.** depression and psychiatric illness; shift work; underlying chronic disease (reverse causation, especially for long sleep); socioeconomic status; obesity and sleep-disordered breathing; age; medications (sedatives)
- **Testing cadence.** baseline + annual; continuous if wearable-tracked
- **Standardization (U-shaped).** Total sleep time has a robust U-shaped relation with mortality and CVD — both short (<6h) and long (>9h) sleep raise risk — so no single per-SD slope is meaningful; the off-scale cells use short-vs-reference and long-vs-reference contrasts.
- **Modifiability (moderate).** Sleep hygiene, cognitive behavioural therapy for insomnia (CBT-I), consistent sleep scheduling, and treatment of underlying disorders (apnea, depression). Extending habitually short sleep toward 7-8h is achievable; shortening pathologically long sleep depends on treating the underlying cause. — *effect:* CBT-I and sleep-extension interventions can increase total sleep time by ~20-60 minutes; effect on hard outcomes is unproven.; *timeframe:* Weeks to a few months for behavioural change.; *evidence:* RCT. Meta-analyses of CBT-I (e.g., van Straten 2018, Sleep Med Rev) and sleep-extension RCTs; outcome trials on mortality are lacking.
- **Outcome evidence.** All-cause mortality native HR 1.12 [categorical, tier A]; CVD native HR 1.11 [categorical, tier A]; Cancer native HR 1.01 [unconvertible, tier B]; Dementia native HR 1.06 [unconvertible, tier B]; Frailty native HR 1.21 [categorical, tier C]

### Sleep efficiency

*Also: SE, time asleep / time in bed, sleep continuity*

- **Mechanism.** Sleep efficiency is the percentage of time in bed actually spent asleep (total sleep time / time in bed x 100). Low efficiency reflects fragmented or difficult-to-initiate/maintain sleep; fragmentation impairs glymphatic clearance, sympathetic down-regulation and metabolic recovery, and predicts mortality and cognitive decline independent of total sleep time.
- **Measurement.** Derived from polysomnography, actigraphy or wearable accelerometry; can be approximated from sleep diaries. Objective measurement is preferred because perceived efficiency diverges from measured. — *sample:* wearable; *approx. cost:* $0-300.
- **Signal quality.** moderate — Objective efficiency is reasonably reproducible across nights but device- and algorithm-dependent; diary-based efficiency is biased. Affected by time-in-bed behaviour, not only sleep quality.
- **Major confounders.** age (efficiency declines with age); depression and anxiety; sleep-disordered breathing; pain and nocturia; time-in-bed behaviour; device/algorithm differences; medications
- **Testing cadence.** baseline + annual; continuous if wearable-tracked
- **Standardization (Clinical cutpoint).** Sleep-efficiency HRs are reported as low-vs-adequate category contrasts (e.g. <80% efficiency, or fragmentation quartiles), not a continuous per-SD exposure.
- **Modifiability (high).** CBT-I (especially sleep-restriction and stimulus-control components, which raise efficiency by design), treatment of pain/nocturia, and treatment of sleep-disordered breathing. — *effect:* CBT-I reliably raises sleep efficiency by ~5-12 percentage points; among the most modifiable sleep parameters.; *timeframe:* 4-8 weeks of CBT-I.; *evidence:* RCT meta-analysis. Trauer JM et al. Cognitive behavioral therapy for chronic insomnia: a systematic review and meta-analysis. Ann Intern Med 2015;163:191-204.
- **Outcome evidence.** All-cause mortality native HR 1.84 [categorical, tier C]; CVD native HR 1.36 [categorical, tier C]; Dementia native HR 1.32 [categorical, tier C]; Frailty native HR 1.40 [categorical, tier C]

### Deep (slow-wave) sleep duration

*Also: slow-wave sleep, N3 sleep, deep sleep, SWS*

- **Mechanism.** Slow-wave sleep (N3) is the deepest non-REM stage, dominated by delta-frequency cortical oscillations. It drives glymphatic clearance of metabolic waste including amyloid-beta, growth-hormone secretion, memory consolidation, parasympathetic dominance and blood-pressure dipping. Slow-wave sleep declines markedly with age; its loss is implicated in neurodegeneration and cardiometabolic risk.
- **Measurement.** Gold standard is polysomnography with EEG staging; consumer wearables estimate 'deep sleep' from accelerometry and heart-rate but agree only modestly with PSG-staged N3. — *sample:* wearable; *approx. cost:* $0-1500.
- **Signal quality.** noisy — PSG N3 staging is reliable but single-night; consumer-wearable 'deep sleep' has poor stage-level agreement with PSG and large between-device variability. Strongly age-dependent.
- **Major confounders.** age (dominant confounder); sleep-disordered breathing; alcohol and sedatives; first-night effect / single-night PSG; device staging error (wearables); depression; prior sleep debt
- **Testing cadence.** baseline + periodic; continuous trend if wearable-tracked
- **Standardization (No usable SD).** Slow-wave-sleep HRs are reported per percentage point or per year of decline with no population SD available, and one cell is a quartile contrast — none can be standardized to per SD.
- **Modifiability (low).** Slow-wave sleep is largely age-determined and not robustly modifiable. Regular aerobic exercise, avoidance of alcohol and certain sedatives near bedtime, and treatment of sleep-disordered breathing modestly preserve or increase slow-wave sleep; acoustic/closed-loop slow-oscillation stimulation can enhance delta power experimentally. — *effect:* Lifestyle measures produce small increases in slow-wave sleep (minutes); no intervention reliably restores youthful slow-wave sleep amounts.; *timeframe:* Weeks to months for exercise effects; acute for acoustic stimulation.; *evidence:* RCT. Exercise-and-sleep RCT meta-analyses; closed-loop acoustic stimulation trials (e.g., Ngo 2013, Neuron). No hard-outcome RCTs.
- **Outcome evidence.** All-cause mortality native HR 1.27 [unconvertible, tier C]; CVD native HR 1.27 [categorical, tier C]; Dementia native HR 1.27 [unconvertible, tier C]

### REM sleep duration

*Also: rapid eye movement sleep, REM%, paradoxical sleep*

- **Mechanism.** REM sleep is characterized by cortical activation, vivid dreaming, muscle atonia and irregular autonomic activity. It supports emotional memory processing, synaptic regulation and thermoregulatory and autonomic functions. Reduced REM has been linked in cohort data to higher mortality and dementia risk; loss of REM atonia (REM behaviour disorder) is a strong prodrome of synucleinopathy.
- **Measurement.** Gold standard is polysomnography with EEG/EOG/EMG staging; consumer wearables estimate REM from heart-rate variability and movement with moderate-to-poor stage agreement vs PSG. — *sample:* wearable; *approx. cost:* $0-1500.
- **Signal quality.** noisy — PSG REM staging is reliable but single-night and sensitive to first-night effects, sleep timing and prior REM deprivation; consumer-wearable REM estimates agree only modestly with PSG.
- **Major confounders.** age; antidepressants and other REM-suppressing drugs; alcohol; sleep-disordered breathing (REM-related desaturation); first-night effect / single-night PSG; depression; circadian timing of recording
- **Testing cadence.** baseline + periodic; continuous trend if wearable-tracked
- **Standardization (No usable SD).** REM-sleep HRs are reported per percentage point of REM with no population SD available to convert them to a per-SD basis.
- **Modifiability (low).** REM duration is largely physiologically and age-determined. Discontinuing or minimizing REM-suppressing medications (many antidepressants) and alcohol near bedtime, treating sleep-disordered breathing, and maintaining adequate total sleep time (REM is concentrated in the last third of the night) modestly preserve REM. — *effect:* Avoiding REM-suppressant exposures and protecting late-night sleep produce modest increases in REM; no intervention reliably and durably augments REM beyond physiologic levels.; *timeframe:* Days to weeks after removing a REM-suppressant.; *evidence:* mechanistic. Pharmacology of REM-suppressing agents; observational sleep-architecture studies. No hard-outcome RCTs targeting REM duration.
- **Outcome evidence.** All-cause mortality native HR 1.13 [unconvertible, tier B]; CVD native HR 1.11 [unconvertible, tier C]; Dementia native HR 1.09 [unconvertible, tier C]

### Sleep apnea severity (ODI / AHI)

*Also: obstructive sleep apnea, OSA, apnea-hypopnea index, AHI, oxygen desaturation index, ODI, home sleep apnea testing*

- **Mechanism.** Obstructive sleep apnea causes repetitive upper-airway collapse with intermittent hypoxia, arousals and intrathoracic pressure swings. Severity is indexed by the apnea-hypopnea index (AHI) or oxygen desaturation index (ODI) — events per hour of sleep. Intermittent hypoxia and sympathetic surges drive hypertension, atrial fibrillation, endothelial dysfunction, insulin resistance, neuroinflammation and oxidative stress, linking severe OSA to cardiovascular, cognitive and mortality risk.
- **Measurement.** Attended in-lab polysomnography (gold standard) or validated home sleep apnea testing (HSAT) recording airflow, respiratory effort and pulse oximetry; ODI is derived from oximetry alone. Severity strata: AHI <5 normal, 5-15 mild, 15-30 moderate, >=30 severe. — *sample:* in-office device; *approx. cost:* $150-3000.
- **Signal quality.** moderate — AHI is reasonably reproducible but varies by hypopnea scoring rule, body position, sleep stage and night-to-night; HSAT can underestimate AHI by dividing events over recording time rather than sleep time. ODI is robust but captures only desaturation events.
- **Major confounders.** obesity (strong shared cause); age and male sex; hypopnea scoring criteria; supine sleep proportion; alcohol and sedatives; single-night variability; REM-related vs non-REM events
- **Testing cadence.** diagnostic test then re-test after intervention or major weight change
- **Standardization (Clinical strata).** Sleep apnea severity is graded into clinical strata (normal/mild/moderate/severe by AHI cut-points) and is highly right-skewed; hard-outcome HRs are reported as severe-vs-none contrasts rather than per-SD slopes.
- **Modifiability (high).** Continuous positive airway pressure (CPAP) is first-line and reliably normalizes AHI/ODI; weight loss (including bariatric surgery and GLP-1 agonists such as tirzepatide), mandibular advancement devices, positional therapy, and hypoglossal-nerve stimulation are alternatives. — *effect:* CPAP reduces AHI to <5 in most adherent users; ~10% weight loss can roughly halve AHI; tirzepatide reduced AHI by ~25-30 events/h in the SURMOUNT-OSA trial.; *timeframe:* Immediate for CPAP; weeks to months for weight-loss approaches.; *evidence:* RCT meta-analysis. CPAP RCTs and meta-analyses (SAVE 2016; meta-analyses on blood pressure and sleepiness); SURMOUNT-OSA tirzepatide RCT (Malhotra 2024, NEJM).
- **Outcome evidence.** All-cause mortality native HR 1.39 [categorical, tier A]; CVD native HR 1.79 [categorical, tier A]; Cancer native HR 1.40 [categorical, tier C]; Dementia native HR 1.26 [categorical, tier B]; Frailty native HR 1.50 [categorical, tier C]

### Wearable-derived sleep regularity

*Also: sleep consistency, sleep regularity index, SRI, sleep timing variability, social jetlag*

- **Mechanism.** Sleep regularity captures the day-to-day consistency of sleep timing and duration. The Sleep Regularity Index (SRI) scores the probability of being in the same state (asleep/awake) at the same clock time on consecutive days (0-100). Irregular sleep misaligns the circadian system from behaviour, disrupting metabolic, autonomic and inflammatory rhythms; emerging cohort data suggest irregularity predicts mortality and cardiometabolic risk independent of sleep duration.
- **Measurement.** Computed from multi-day accelerometry or wearable sleep tracking (typically 7+ days). Common metrics: Sleep Regularity Index, standard deviation of sleep onset/midpoint, and social jetlag (weekday-weekend midpoint difference). — *sample:* wearable; *approx. cost:* $0-300.
- **Signal quality.** moderate — Requires multiple consecutive nights, making it more stable than single-night metrics, but depends on device sleep-detection accuracy and the chosen regularity metric; definitions are not yet standardized across studies.
- **Major confounders.** shift work and irregular work schedules; age; depression and psychiatric illness; caregiving and social obligations; device/algorithm and metric definition differences; underlying illness; weekend/weekday measurement window
- **Testing cadence.** rolling multi-day window; continuous if wearable-tracked
- **Standardization (Index extremes).** Sleep-regularity HRs contrast extreme percentiles of the Sleep Regularity Index (most vs least regular), an index-based contrast rather than a per-SD continuum.
- **Modifiability (moderate).** Consistent sleep-wake scheduling (fixed bedtime and wake time, including weekends), light exposure timing, limiting late-night screen/social activity, and managing shift-work schedules. Wearable feedback and CBT-I components can support regular timing. — *effect:* Behavioural scheduling can substantially raise the Sleep Regularity Index and cut social jetlag within weeks; durability depends on work and social constraints.; *timeframe:* Days to weeks for timing changes; sustained adherence is the limiting factor.; *evidence:* observational. Circadian and sleep-timing intervention studies; no hard-outcome RCTs targeting sleep regularity specifically.
- **Outcome evidence.** All-cause mortality native HR 1.53 [categorical, tier B]; CVD native HR 1.40 [categorical, tier B]; Cancer native HR 1.28 [categorical, tier C]; Dementia native HR 1.53 [categorical, tier C]


## Micronutrient & Metabolic Cofactors

### Vitamin D (25-OH)

*Also: 25-hydroxyvitamin D, 25(OH)D, calcidiol, serum vitamin D*

- **Mechanism.** 25-hydroxyvitamin D is the circulating storage form and the standard index of vitamin D status. Beyond calcium/bone homeostasis, the vitamin D receptor is expressed in immune, vascular and muscle tissue; deficiency is associated with impaired immune regulation, sarcopenia and adverse cardiometabolic profiles. Observational associations are strong but largely non-causal: Mendelian randomization and large supplementation RCTs show little effect on hard outcomes except in frank deficiency.
- **Measurement.** Immunoassay or LC-MS/MS (reference method) on serum; standardized to the Vitamin D Standardization Program. — *sample:* venous blood; *approx. cost:* $30-65.
- **Signal quality.** moderate — Marked seasonal and assay-platform variability; immunoassays differ from LC-MS/MS. Status is heavily confounded by adiposity, sun exposure, skin pigmentation and outdoor activity, so a low value is partly a marker of poor general health.
- **Major confounders.** season and latitude; adiposity (volumetric dilution); skin pigmentation and sun exposure; outdoor physical activity; supplement use; assay platform differences; reverse causation (illness reduces outdoor activity)
- **Testing cadence.** baseline; recheck after supplementation or seasonally if deficient
- **Standardization (U-shaped).** 25(OH)D shows a non-linear, reverse-J/U-shaped relation with mortality (risk concentrated at deficiency; little or no gradient at high-normal); single linear per-SD HRs understate the deficiency-end effect.
- **Modifiability (moderate).** Oral vitamin D3 supplementation; sunlight exposure. Correcting frank deficiency raises 25(OH)D reliably. — *effect:* 800-2000 IU/day typically raises 25(OH)D by ~20-40 nmol/L; status is highly modifiable, but hard-outcome benefit is largely confined to correcting genuine deficiency.; *timeframe:* 8-12 weeks to a new steady state; *evidence:* RCT meta-analysis. VITAL (Manson 2019, NEJM), ViDA, D-Health and the VITAL/D-Health mortality meta-analyses: 25(OH)D is highly modifiable but supplementation did NOT reduce all-cause mortality, CVD events or cancer incidence in vitamin-D-replete populations.
- **Outcome evidence.** All-cause mortality native HR 1.25 [categorical, tier A]; CVD HR/SD 1.00 (tier A); Cancer HR/SD 1.00 (tier A); Dementia native HR 1.19 [categorical, tier B]; Frailty native HR 1.30 [categorical, tier B]

### Vitamin B12

*Also: cobalamin, serum B12, holotranscobalamin*

- **Mechanism.** Vitamin B12 (cobalamin) is a cofactor for methionine synthase and methylmalonyl-CoA mutase; deficiency causes megaloblastic anemia, peripheral neuropathy, cognitive impairment and elevated homocysteine and methylmalonic acid. Elevated serum B12 is rarely a true nutritional state and usually flags hepatic release, myeloproliferative/neoplastic disease or renal impairment, which explains the high-end mortality signal.
- **Measurement.** Competitive-binding immunoassay for total serum B12; holotranscobalamin and functional markers (MMA, homocysteine) used to confirm deficiency. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** noisy — Total serum B12 is a poor index of tissue status: it includes inactive haptocorrin-bound B12, so normal values can mask functional deficiency. Confirmation with MMA/holoTC is often needed; high values are frequently epiphenomena of disease.
- **Major confounders.** liver disease (elevates B12); malignancy and myeloproliferative disorders (elevate B12); renal impairment; supplement and fortified-food intake; metformin and PPI use (lower B12); pregnancy; haptocorrin-bound inactive fraction
- **Testing cadence.** baseline; recheck if symptomatic or on metformin/PPI long-term
- **Standardization (U-shaped).** Serum B12 has a U-shaped mortality association: both deficiency and high/elevated B12 (often reflecting liver disease, malignancy or renal impairment) carry excess risk, so no linear per-SD HR is defensible.
- **Modifiability (moderate).** Oral or intramuscular B12 (cobalamin) for deficiency; address malabsorption (pernicious anemia, metformin/PPI effect). — *effect:* Replacement reliably normalizes serum B12 and corrects deficiency-related anemia/neuropathy; high baseline B12 is not a treatment target (it reflects underlying disease).; *timeframe:* Weeks (hematologic response); months for neurologic recovery; *evidence:* RCT meta-analysis. B-Vitamin Treatment Trialists' Collaboration (homocysteine-lowering trials); B12 deficiency correction is effective for the deficiency syndrome, but supplementation in replete people has not reduced cardiovascular events, cancer or dementia.
- **Outcome evidence.** All-cause mortality native HR 1.00 [categorical, tier B]; CVD (documented, no pooled HR); Cancer native HR 3.00 [categorical, tier B]; Dementia HR/SD 1.00 (tier B); Frailty (documented, no pooled HR)

### Folate

*Also: serum folate, RBC folate, vitamin B9, 5-methyltetrahydrofolate*

- **Mechanism.** Folate is the methyl donor for one-carbon metabolism, required for DNA synthesis/repair, methionine regeneration and homocysteine remethylation. Deficiency causes megaloblastic anemia, elevated homocysteine and (periconceptionally) neural tube defects. Population folate status improved markedly after grain fortification, compressing the deficiency range in fortified countries.
- **Measurement.** Immunoassay for serum folate (recent intake) or red-cell folate (longer-term stores); LC-MS/MS for 5-MTHF. — *sample:* venous blood; *approx. cost:* $20-45.
- **Signal quality.** moderate — Serum folate fluctuates with recent diet; RBC folate is a steadier index of stores. Interpretation is intertwined with B12 status (masking of B12 deficiency) and varies between fortified and non-fortified populations.
- **Major confounders.** recent dietary intake (serum folate); grain fortification status of the population; concurrent B12 status; alcohol use; MTHFR genotype; supplement use; hemolysis (affects RBC folate assay)
- **Testing cadence.** baseline; recheck if anemic or planning pregnancy
- **Standardization (Trial contrast).** The folate cancer cell is a supplementation-vs-placebo randomized contrast (folic acid did not significantly change cancer incidence), not a per-SD estimate of circulating folate.
- **Modifiability (high).** Dietary folate (leafy greens, legumes), folic acid supplements, and population-level grain fortification. — *effect:* Fortification and supplementation reliably raise serum/RBC folate and lower homocysteine ~20-25%; status is highly modifiable.; *timeframe:* 4-8 weeks for serum folate; ~3-4 months for RBC folate; *evidence:* RCT meta-analysis. B-Vitamin Treatment Trialists' Collaboration (Clarke 2010) and Vollset 2013 cancer meta-analysis: folate is highly modifiable but supplementation did not reduce cardiovascular events or cancer; periconceptional benefit (neural tube defects) is the clear exception.
- **Outcome evidence.** All-cause mortality HR/SD 1.06 (tier B); CVD HR/SD 1.00 (tier A); Cancer native HR 1.06 [categorical, tier A]; Dementia HR/SD 1.00 (tier B); Frailty (documented, no pooled HR)

### Methylmalonic acid (MMA)

*Also: MMA, serum methylmalonic acid, functional B12 marker*

- **Mechanism.** Methylmalonic acid accumulates when methylmalonyl-CoA mutase activity is impaired by functional vitamin B12 deficiency; it is the most specific functional marker of cellular B12 status. Elevated MMA also rises with renal impairment and ageing, and emerging work links circulating MMA to mitochondrial dysfunction and a tumor-progression-promoting metabolic state.
- **Measurement.** LC-MS/MS or GC-MS quantification of serum/plasma methylmalonic acid. — *sample:* venous blood; *approx. cost:* $60-150.
- **Signal quality.** moderate — MMA is a specific functional marker of B12 status but is strongly raised by reduced renal clearance and rises with age independent of B12, so an elevated value must be interpreted against eGFR.
- **Major confounders.** renal function / eGFR (major); age; true functional B12 deficiency; gut microbiome-derived propionate; intra-individual variability
- **Testing cadence.** as needed to confirm functional B12 deficiency
- **Modifiability (moderate).** Vitamin B12 replacement lowers MMA when elevation is due to functional B12 deficiency; renally driven elevation is not nutritionally modifiable. — *effect:* B12 repletion reliably lowers MMA when it is the cause; the portion driven by reduced renal clearance or ageing is not modifiable.; *timeframe:* Weeks after B12 repletion; *evidence:* observational. B12 supplementation studies show MMA falls with repletion; no outcome trial has tested whether lowering MMA itself reduces mortality.
- **Outcome evidence.** All-cause mortality HR/SD 1.20 (tier B); CVD HR/SD 1.16 (tier C); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)

### Omega-3 index

*Also: O3I, RBC EPA+DHA, erythrocyte omega-3 index*

- **Mechanism.** The omega-3 index is the percentage of EPA + DHA in red-cell membrane fatty acids, an integrated index of long-chain omega-3 status over ~3-4 months. Higher membrane omega-3 content is associated with lower triglycerides, reduced platelet aggregation, anti-arrhythmic membrane effects and lower systemic inflammation; it is a candidate modifiable substrate for cardiovascular and possibly brain ageing.
- **Measurement.** Gas chromatography of red-cell membrane fatty acids, expressed as % EPA+DHA of total fatty acids. — *sample:* venous blood; *approx. cost:* $60-120.
- **Signal quality.** moderate — RBC omega-3 index is more stable than plasma fatty acids (reflects ~120-day membrane turnover) but assay standardization across labs is imperfect and the index varies with diet, fish/supplement intake and genetics.
- **Major confounders.** fish and seafood intake; fish-oil supplement use; overall diet quality; genetics of fatty-acid desaturation (FADS locus); assay/lab standardization
- **Testing cadence.** baseline; recheck ~3-4 months after dietary change
- **Modifiability (high).** Increased oily-fish intake or EPA/DHA supplementation reliably raises the red-cell omega-3 index. — *effect:* ~1-2 g/day EPA+DHA typically raises the index by 3-5 percentage points (e.g., from ~4% to ~8%); dose-dependent and durable.; *timeframe:* ~3-4 months to a new steady state (RBC membrane turnover); *evidence:* RCT meta-analysis. Dose-response supplementation trials of EPA/DHA and the omega-3 index; outcome trials (VITAL, STRENGTH, REDUCE-IT) show the index is highly modifiable but routine-dose supplementation has mostly NOT reduced hard CVD or mortality endpoints outside high-triglyceride/high-dose contexts.
- **Outcome evidence.** All-cause mortality HR/SD 0.86 (tier A); CVD HR/SD 0.90 (tier A); Cancer HR/SD 0.94 (tier B); Dementia HR/SD 0.90 (tier B); Frailty HR/SD 0.92 (tier C)

### Magnesium (RBC)

*Also: red blood cell magnesium, erythrocyte magnesium, serum magnesium*

- **Mechanism.** Magnesium is a cofactor for >300 enzymes, including those of ATP metabolism, and stabilizes membranes and cardiac electrical activity. Most serum-magnesium epidemiology actually uses serum (not RBC) magnesium; low magnesium status is linked to insulin resistance, hypertension, arrhythmia and vascular calcification. RBC magnesium is a somewhat better index of intracellular stores but is far less studied in outcome cohorts.
- **Measurement.** Colorimetric or atomic-absorption assay; serum magnesium is routine, RBC magnesium requires red-cell isolation. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** noisy — Serum magnesium reflects only ~1% of body magnesium and is tightly homeostatically buffered; RBC magnesium better tracks stores but is poorly standardized and rarely used in outcome studies, so most evidence rests on serum magnesium.
- **Major confounders.** renal function; diuretic and PPI use; diet quality and alcohol; diabetes / insulin resistance; homeostatic buffering of serum magnesium; hemolysis (affects RBC assay)
- **Testing cadence.** baseline; recheck if on diuretics/PPIs or symptomatic
- **Modifiability (moderate).** Magnesium-rich diet (whole grains, nuts, leafy greens, legumes) and oral magnesium supplements; address depleting medications (diuretics, PPIs). — *effect:* Supplementation modestly raises serum and RBC magnesium and lowers blood pressure ~2-4 mmHg; deficiency correction is reliable.; *timeframe:* Weeks to a few months; *evidence:* RCT meta-analysis. RCT meta-analyses of magnesium supplementation and blood pressure / insulin sensitivity show modest physiologic effects; no large trial demonstrates reduced mortality, CVD events or dementia.
- **Outcome evidence.** All-cause mortality HR/SD 1.00 (tier B); CVD HR/SD 0.94 (tier A); Cancer (documented, no pooled HR); Dementia HR/SD 0.96 (tier C); Frailty (documented, no pooled HR)

### Zinc

*Also: serum zinc, plasma zinc*

- **Mechanism.** Zinc is a structural and catalytic cofactor for hundreds of enzymes and transcription factors, essential for immune function, antioxidant defense (Cu/Zn-SOD) and protein synthesis. Low zinc status is linked to impaired immunity, poor wound healing and, in older adults, to sarcopenia and infection risk. Serum zinc is a weak status index because it is homeostatically buffered and falls as an acute-phase response.
- **Measurement.** Atomic-absorption spectroscopy or ICP-MS on serum/plasma; trace-metal-free collection required. — *sample:* venous blood; *approx. cost:* $25-60.
- **Signal quality.** noisy — Serum zinc reflects a small, tightly buffered pool and drops during inflammation/infection (redistribution), so a low value often signals illness rather than dietary deficiency. Diurnal variation and fasting state also affect it.
- **Major confounders.** acute-phase response / inflammation (lowers zinc); fasting state and time of day; hypoalbuminemia; trace-metal contamination of collection tubes; diet and supplement use
- **Testing cadence.** baseline; interpret alongside CRP
- **Modifiability (moderate).** Zinc-rich diet (meat, shellfish, legumes, nuts) and oral zinc supplements; address malabsorption and excess copper/iron competition. — *effect:* Supplementation reliably corrects deficiency; serum zinc is buffered so changes are modest, and excess zinc can cause copper deficiency.; *timeframe:* Weeks to months; *evidence:* observational. Zinc supplementation trials show benefit for deficiency syndromes (immune function, growth) but no large trial demonstrates reduced mortality, CVD or dementia; over-supplementation risks copper depletion.
- **Outcome evidence.** All-cause mortality HR/SD 1.11 (tier B); CVD HR/SD 1.09 (tier C); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR); Frailty HR/SD 1.11 (tier C)

### Copper

*Also: serum copper, plasma copper*

- **Mechanism.** Copper is a cofactor for cytochrome c oxidase, Cu/Zn-SOD, lysyl oxidase and ceruloplasmin. ~90% of serum copper is carried on ceruloplasmin, an acute-phase reactant, so serum copper rises with inflammation, infection, malignancy and estrogen exposure. Elevated serum copper is therefore largely a marker of inflammatory or neoplastic disease rather than copper intake.
- **Measurement.** Atomic-absorption spectroscopy or ICP-MS on serum/plasma. — *sample:* venous blood; *approx. cost:* $25-60.
- **Signal quality.** noisy — Serum copper is dominated by ceruloplasmin-bound copper and behaves as an acute-phase reactant; a high value usually reflects inflammation, malignancy, pregnancy or estrogen use, not dietary copper. Poor index of true copper status.
- **Major confounders.** inflammation / acute-phase response (raises copper); malignancy; estrogen / oral contraceptives / pregnancy; liver disease; ceruloplasmin level; smoking
- **Testing cadence.** baseline; interpret alongside CRP and ceruloplasmin
- **Standardization (U-shaped).** Serum copper shows a non-monotonic mortality relation: both deficiency and (more strongly) elevated copper carry risk, and high copper is largely an acute-phase / inflammation epiphenomenon, so a linear per-SD HR is not defensible.
- **Modifiability (low).** Copper status is rarely a clinical target; deficiency (rare, e.g., excess zinc, bariatric surgery, malabsorption) is corrected with copper repletion. Elevated serum copper reflects inflammation/disease and is not itself a modifiable nutritional target. — *effect:* Deficiency correction is reliable; elevated serum copper falls only when the underlying inflammatory or neoplastic cause resolves.; *timeframe:* Weeks (deficiency correction); *evidence:* mechanistic. No outcome trial shows benefit from lowering serum copper in the general population; high serum copper is an epiphenomenon of inflammation/malignancy, not a treatment target.
- **Outcome evidence.** All-cause mortality native HR 1.50 [categorical, tier B]; CVD native HR 1.40 [categorical, tier B]; Cancer native HR 1.45 [categorical, tier C]; Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)

### Ceruloplasmin

*Also: serum ceruloplasmin, ferroxidase*

- **Mechanism.** Ceruloplasmin is the main copper-carrying plasma protein and a ferroxidase essential for iron loading onto transferrin; it is also a positive acute-phase reactant. Elevated ceruloplasmin reflects inflammation, estrogen exposure or copper status; very low ceruloplasmin is a hallmark of Wilson disease and aceruloplasminemia. As a biomarker it largely overlaps with serum copper and inflammatory signaling.
- **Measurement.** Immunoturbidimetric/nephelometric assay (protein mass) or enzymatic ferroxidase activity assay. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** noisy — Ceruloplasmin behaves as an acute-phase reactant and rises with estrogen, pregnancy and inflammation; immunoassays measure both holo- and apo-protein, so values do not cleanly index copper status or any single pathway.
- **Major confounders.** inflammation / acute-phase response; estrogen / oral contraceptives / pregnancy; liver synthetic function; copper status; smoking
- **Testing cadence.** baseline; mainly used diagnostically (e.g., Wilson disease workup)
- **Modifiability (low).** Ceruloplasmin is not a stand-alone intervention target; it is used diagnostically (Wilson disease, copper deficiency). Treating an underlying inflammatory state lowers reactive elevation. — *effect:* Reactive elevation falls when the inflammatory cause resolves; there is no rationale for targeting ceruloplasmin per se.; *timeframe:* Weeks to months (with resolution of inflammation); *evidence:* mechanistic. No outcome trial targets ceruloplasmin; its elevation in the general population is an inflammatory epiphenomenon.
- **Outcome evidence.** All-cause mortality HR/SD 1.17 (tier C); CVD HR/SD 1.20 (tier B); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)

### Transferrin saturation

*Also: TSAT, transferrin saturation %, iron saturation*

- **Mechanism.** Transferrin saturation is the percentage of transferrin iron-binding sites occupied (serum iron / TIBC x 100); it indexes iron availability for erythropoiesis and tissue delivery. Low TSAT indicates iron deficiency; high TSAT indicates iron overload (e.g., hemochromatosis) with risk of catalytic free-iron-driven oxidative tissue damage. It is a more dynamic iron-status marker than ferritin and less inflammation-confounded.
- **Measurement.** Calculated from serum iron and total iron-binding capacity (or transferrin); colorimetric/immunoturbidimetric assays. — *sample:* venous blood; *approx. cost:* $15-40.
- **Signal quality.** moderate — TSAT has marked diurnal variation (serum iron peaks in the morning) and responds to recent intake and inflammation; a single measurement is noisy, but it is less acute-phase-confounded than ferritin and tracks functional iron availability.
- **Major confounders.** diurnal variation in serum iron; recent dietary iron / supplements; inflammation (lowers serum iron); fasting state; HFE genotype; menstrual blood loss
- **Testing cadence.** baseline (morning, fasting); part of an iron panel with ferritin and TIBC
- **Standardization (U-shaped).** Transferrin saturation has a U-shaped mortality relation: both low TSAT (iron deficiency) and high TSAT (iron overload) carry excess risk, so no single linear per-SD HR is defensible.
- **Modifiability (high).** Low TSAT (iron deficiency): oral or IV iron repletion and treatment of blood loss. High TSAT (iron overload/hemochromatosis): therapeutic phlebotomy; address dietary/supplemental iron excess. — *effect:* Iron repletion reliably raises TSAT and corrects deficiency; phlebotomy reliably lowers iron stores in overload.; *timeframe:* Weeks to months; *evidence:* RCT. IV iron RCTs in heart failure (e.g., FAIR-HF, AFFIRM-AHF) improve symptoms/hospitalization in iron-deficient patients; phlebotomy is standard for hemochromatosis, though a general-population iron-reduction trial (VA PAD trial) was largely null for hard outcomes.
- **Outcome evidence.** All-cause mortality native HR 1.20 [categorical, tier B]; CVD native HR 1.25 [categorical, tier B]; Cancer native HR 1.20 [categorical, tier C]; Dementia (documented, no pooled HR); Frailty native HR 1.30 [categorical, tier C]

### TIBC

*Also: total iron-binding capacity, transferrin (as TIBC)*

- **Mechanism.** Total iron-binding capacity reflects the iron-binding capacity of circulating transferrin and rises in iron deficiency (transferrin upregulation) and falls in iron overload, inflammation and malnutrition. As transferrin is a negative acute-phase reactant and a marker of hepatic synthetic/nutritional status, low TIBC also flags inflammation, liver disease and protein-energy malnutrition -- making it partly a general illness marker.
- **Measurement.** Direct colorimetric assay or calculated from serum iron and unsaturated iron-binding capacity; transferrin immunoassay as an equivalent. — *sample:* venous blood; *approx. cost:* $15-40.
- **Signal quality.** moderate — TIBC moves inversely with inflammation and tracks nutritional/hepatic status as well as iron stores, so a low value is a mixed signal of iron overload, inflammation or malnutrition; best interpreted within a full iron panel.
- **Major confounders.** inflammation (lowers TIBC/transferrin); protein-energy malnutrition; liver disease; iron stores; pregnancy and estrogen (raise transferrin); nephrotic syndrome
- **Testing cadence.** baseline; part of an iron panel with ferritin and transferrin saturation
- **Modifiability (low).** TIBC is not itself a treatment target. It rises with iron deficiency (and falls with iron repletion) and falls with inflammation/malnutrition, so it changes only as a downstream consequence of treating those states. — *effect:* TIBC normalizes when the underlying iron, nutritional or inflammatory state is corrected; there is no rationale for targeting TIBC directly.; *timeframe:* Weeks to months (with treatment of the underlying state); *evidence:* mechanistic. No outcome trial targets TIBC; it is interpreted only within an iron panel and as a nutritional/inflammatory marker.
- **Outcome evidence.** All-cause mortality HR/SD 1.14 (tier B); CVD HR/SD 1.11 (tier C); Cancer (documented, no pooled HR); Dementia (documented, no pooled HR); Frailty (documented, no pooled HR)

### Selenium

*Also: serum selenium, plasma selenium, selenoprotein P*

- **Mechanism.** Selenium is incorporated into ~25 selenoproteins, including glutathione peroxidases and thioredoxin reductases central to antioxidant defense, and selenoprotein P (the main transport form). Adequate selenium supports redox homeostasis, immune function and thyroid hormone metabolism. Both deficiency and excess are harmful, and the dose-response is narrow, giving a U-shaped relation with mortality.
- **Measurement.** ICP-MS or atomic-absorption spectroscopy for serum/plasma selenium; selenoprotein P immunoassay as a functional index. — *sample:* venous blood; *approx. cost:* $40-90.
- **Signal quality.** moderate — Serum selenium reflects recent intake and varies enormously by geography (soil selenium); plateaus once selenoproteins are saturated, so above the saturation point higher serum selenium no longer indexes 'better' status and may indicate excess.
- **Major confounders.** geographic/soil selenium and diet; baseline population selenium status; inflammation (lowers selenium); supplement use; selenoprotein saturation plateau
- **Testing cadence.** baseline; mainly relevant in low-selenium regions
- **Standardization (U-shaped).** Selenium has a U-shaped mortality relation: deficiency and excess both carry risk, with benefit largely confined to a narrow range and only in populations with low baseline status; a single linear per-SD HR is not defensible.
- **Modifiability (moderate).** Dietary selenium (Brazil nuts, seafood, organ meats) and supplements; chiefly relevant in low-selenium regions. Avoid excess (narrow therapeutic window). — *effect:* Supplementation reliably raises serum selenium and saturates selenoproteins; status is highly modifiable, but the safe-and-beneficial range is narrow.; *timeframe:* Weeks to months; *evidence:* RCT meta-analysis. SELECT (Lippman 2009, JAMA), NPC trial, and Cochrane reviews: selenium is highly modifiable but supplementation in replete populations did NOT reduce cancer, CVD or mortality and increased type 2 diabetes risk.
- **Outcome evidence.** All-cause mortality native HR 0.76 [categorical, tier B]; CVD HR/SD 1.00 (tier A); Cancer HR/SD 1.00 (tier A); Dementia (documented, no pooled HR); Frailty HR/SD 1.13 (tier C)

### TMAO (with choline)

*Also: trimethylamine N-oxide, TMAO, choline, betaine, gut-microbiome metabolite*

- **Mechanism.** Trimethylamine N-oxide (TMAO) is a gut-microbiome-dependent metabolite: dietary choline, phosphatidylcholine, betaine and L-carnitine are converted by gut bacteria to trimethylamine, then oxidized by hepatic FMO3 to TMAO. Elevated circulating TMAO is associated with atherosclerosis, thrombosis (enhanced platelet reactivity), renal impairment and adverse cardiovascular outcomes; choline is the dietary precursor and is itself an essential nutrient. The marker integrates diet, microbiome composition and renal clearance.
- **Measurement.** LC-MS/MS quantification of plasma TMAO (and, when paired, choline and betaine). — *sample:* venous blood; *approx. cost:* $80-200.
- **Signal quality.** noisy — TMAO has high intra-individual variability (diet- and microbiome-dependent), is strongly raised by reduced renal clearance, and a single measurement is an unstable index; the causal interpretation versus renal-function confounding is debated.
- **Major confounders.** renal function / eGFR (major); recent diet (red meat, fish, eggs); gut microbiome composition; FMO3 genotype/activity; fasting state; high intra-individual variability
- **Testing cadence.** not routine; research/specialty use
- **Modifiability (moderate).** Reduce dietary precursors (red meat, egg yolk, high-dose carnitine/choline supplements); plant-forward diets and microbiome modulation lower TMAO. Renal-function support reduces accumulation. Choline itself remains an essential nutrient and should not be over-restricted. — *effect:* Dietary change can lower plasma TMAO ~20-50% over weeks, and a vegan/vegetarian microbiome produces little TMAO even after a precursor challenge.; *timeframe:* Weeks (diet/microbiome-dependent); *evidence:* observational. Feeding/dietary-intervention studies show TMAO is readily lowered by reducing animal-product intake; however, NO outcome trial has shown that lowering TMAO reduces mortality, CVD or dementia -- modifiability of the marker is established, but hard-outcome benefit is not.
- **Outcome evidence.** All-cause mortality HR/SD 1.18 (tier A); CVD HR/SD 1.23 (tier A); Cancer (documented, no pooled HR); Dementia HR/SD 1.07 (tier C); Frailty HR/SD 1.13 (tier C)


## Cancer Screening & Risk

### Multi-cancer early detection test (Galleri)

*Also: Galleri, MCED, multi-cancer early detection, cfDNA methylation cancer test, GRAIL Galleri*

- **Mechanism.** Galleri analyzes cell-free DNA methylation patterns in plasma to detect a shared cancer signal across >50 cancer types and predict the tissue of origin. The intent is to detect cancers earlier — including cancers with no standard screening — so treatment begins at a more curable stage; the longevity rationale is stage-shift toward earlier diagnosis, not measurement of an aging axis.
- **Measurement.** Single blood draw; targeted bisulfite sequencing of cell-free DNA methylation with a machine-learning classifier returning a binary cancer-signal-detected/not-detected call plus predicted tissue of origin. — *sample:* venous blood; *approx. cost:* $750-950.
- **Signal quality.** moderate — High specificity (~99.5%) keeps false positives low, but overall sensitivity for early-stage disease is limited (~17% stage I, ~40% stage II) and varies widely by cancer type; performance is a test characteristic, not a graded biological signal.
- **Major confounders.** low sensitivity for early-stage and indolent cancers; tissue-of-origin misassignment; lead-time and overdiagnosis bias inflate apparent benefit; no completed evidence of mortality reduction; clonal hematopoiesis and benign conditions as signal sources
- **Testing cadence.** marketed as annual in adults at elevated risk (typically age 50+); optimal interval not established by trial
- **Modifiability (fixed).** Not a modifiable biomarker — Galleri is a diagnostic test. The actionable response to a cancer-signal-detected result is confirmatory diagnostic workup directed by the predicted tissue of origin; the test result itself is not a target to be improved. — *effect:* Not applicable — a test result, not a physiological quantity.; *timeframe:* not applicable; *evidence:* observational. PATHFINDER (Schrag 2023); NHS-Galleri trial protocol (Neal RD et al., ISRCTN91431511).
- **Outcome evidence.** Cancer (documented, no pooled HR)

### PSA (prostate-specific antigen)

*Also: prostate-specific antigen, total PSA, serum PSA*

- **Mechanism.** PSA is a kallikrein protease secreted by prostate epithelium; serum levels rise with prostate volume, inflammation and, importantly, neoplasia. A single PSA measured in mid-life is a strong long-range marker of underlying prostate-cancer biology — men in the top of the distribution at age 45-55 carry most of the long-term risk of clinically significant and fatal prostate cancer.
- **Measurement.** Immunoassay for total PSA on a standard chemistry/immunoassay analyzer; free-PSA and derived indices available; results in ng/mL. — *sample:* venous blood; *approx. cost:* $20-50.
- **Signal quality.** moderate — Reliable, standardized assay, but PSA is prostate- not cancer-specific: benign hyperplasia, prostatitis, recent ejaculation, instrumentation and 5-alpha-reductase inhibitors all shift levels. A single mid-life value is nonetheless a robust long-term risk marker; serial population screening carries overdiagnosis risk.
- **Major confounders.** benign prostatic hyperplasia and prostatitis; recent ejaculation, cycling, digital rectal exam or instrumentation; 5-alpha-reductase inhibitors (finasteride/dutasteride roughly halve PSA); prostate volume and age; obesity (hemodilution lowers measured PSA)
- **Testing cadence.** baseline mid-life value (age ~45-50) for risk stratification; risk-adapted interval screening thereafter (e.g., every 2-4 years if low, less often if very low)
- **Standardization (Quantile / cutpoint).** PSA evidence here is a long-term predictive contrast — a single mid-life PSA value stratified by percentile or clinical cutpoint against future prostate-cancer outcomes; PSA is highly right-skewed and the cells use extreme-group contrasts rather than a smooth per-SD slope.
- **Modifiability (fixed).** PSA is a diagnostic/risk marker rather than a longevity target. 5-alpha-reductase inhibitors lower PSA ~50% (a measurement artifact, not risk reduction); the actionable response to an elevated or rising mid-life PSA is risk-adapted surveillance, MRI and biopsy as indicated. — *effect:* Not a modifiable physiological quantity; 5-ARIs roughly halve measured PSA but this reflects assay shift, not changed prostate-cancer biology.; *timeframe:* not applicable; *evidence:* observational. Malmo Preventive Project (Vickers 2013, BMJ); ERSPC screening trial (Schroder FH et al., Lancet 2014).
- **Outcome evidence.** Cancer native HR 7.70 [categorical, tier B]

### Age-appropriate screening compliance

*Also: cancer screening adherence, guideline-concordant screening, up-to-date screening status*

- **Mechanism.** Being up to date with age-appropriate cancer screening (e.g., colorectal, breast, cervical, and lung screening where eligible) shifts cancer diagnosis to earlier, more curable stages and, for colorectal and cervical screening, removes precancerous lesions. The longevity-relevant effect is reduced cancer-specific mortality; effects on all-cause mortality are small and inconsistent across trials.
- **Measurement.** Self-report or claims/registry-based classification of whether an individual is up to date with each guideline-recommended screen for their age and risk; sometimes summarized as a composite adherence index. — *sample:* functional test; *approx. cost:* $0-50.
- **Signal quality.** noisy — A behavioural composite that varies by which screens are included and how 'up to date' is defined; heavily confounded by the healthy-screenee effect (people who screen also have healthier behaviours and better access), so observational mortality associations overstate causal benefit relative to randomized screening trials.
- **Major confounders.** healthy-screenee / healthy-adherer bias; socioeconomic status and healthcare access; general health-seeking behaviour and comorbidity; lead-time and overdiagnosis bias; heterogeneity in which screens are counted
- **Testing cadence.** ongoing — status reassessed against each guideline-recommended screening interval
- **Standardization (Behaviour).** Screening compliance is an adherence behaviour (up-to-date vs not for guideline-recommended cancer screens), not a continuous biomarker — its evidence is screening-trial cancer-mortality reductions and observational adherence contrasts, expressed as categorical risk ratios that cannot be standardized to a per-SD scale.
- **Modifiability (high).** Highly modifiable behaviour: patient reminders/recall systems, mailed FIT or self-sampling kits, navigation services, and removing cost barriers substantially raise screening uptake. — *effect:* Mailed FIT outreach and patient navigation raise colorectal screening completion by ~20-30 absolute percentage points versus usual care in underserved populations.; *timeframe:* weeks to months to complete an overdue screen; ongoing to sustain up-to-date status; *evidence:* RCT meta-analysis. Cochrane and USPSTF reviews of interventions to increase cancer-screening uptake; mailed-FIT and patient-navigation randomized trials.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); Cancer native HR 0.78 [categorical, tier A]


## Genetic & Pharmacogenomic

### MTHFR variants (C677T)

*Also: MTHFR C677T, rs1801133, methylenetetrahydrofolate reductase polymorphism, MTHFR 677 TT genotype, MTHFR A1298C*

- **Mechanism.** MTHFR encodes methylenetetrahydrofolate reductase, which converts 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate, the methyl donor for remethylation of homocysteine to methionine. The C677T (rs1801133) variant produces a thermolabile enzyme with reduced activity; TT homozygotes have ~25-30% higher plasma homocysteine, an effect amplified when folate intake is low. Mildly elevated homocysteine has been proposed as a vascular and neurodegenerative risk factor, though folate-lowering trials have not reduced clinical events.
- **Measurement.** Single-SNP genotyping (rs1801133, often paired with rs1801131/A1298C) by PCR, TaqMan assay, or array; performed once, results are lifelong. — *sample:* venous blood; *approx. cost:* $30-150.
- **Signal quality.** clean — Genotype is fixed and measured with essentially no biological or analytic variability. The downstream phenotype (homocysteine) is highly modifiable by folate status, so the genotype's clinical meaning depends heavily on context.
- **Major confounders.** none for the genotype itself (fixed); folate and B12 status strongly modify the homocysteine phenotype; mandatory folic acid fortification attenuates genotype effect; ancestry-dependent allele frequency
- **Testing cadence.** once in a lifetime (genotype is fixed)
- **Standardization (Genotype).** MTHFR C677T is a genotype (CC/CT/TT), not a continuous measurement — there is no standard deviation to standardize against; HRs are TT-vs-CC group contrasts and are folate-status dependent.
- **Modifiability (fixed).** Genotype is fixed and cannot be changed. The downstream homocysteine phenotype is modifiable: adequate folate (including 5-methyltetrahydrofolate), B12 and B6 normalise homocysteine in TT homozygotes. — *effect:* Folate repletion can lower elevated homocysteine by ~25% or more, but homocysteine-lowering has not translated into reduced cardiovascular or mortality events in RCTs.; *timeframe:* weeks for homocysteine normalisation; genotype itself never changes; *evidence:* RCT meta-analysis. B-vitamin homocysteine-lowering trial meta-analyses (e.g., VITATOPS, SEARCH, HOPE-2); MTHFR genotype biology.
- **Outcome evidence.** All-cause mortality native HR 0.79 [categorical, tier C]; CVD native HR 1.16 [categorical, tier B]

### Pharmacogenomic panel (CYP2D6, CYP2C19, etc.)

*Also: PGx panel, drug-metabolism genotyping, CYP2D6 genotype, CYP2C19 genotype, preemptive pharmacogenomic testing, CYP2C9 / VKORC1 / SLCO1B1 / TPMT panel*

- **Mechanism.** A pharmacogenomic panel genotypes drug-metabolizing enzymes and transporters (CYP2D6, CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, and others) to classify a person as poor, intermediate, normal, rapid, or ultrarapid metabolizer for specific drugs. It predicts DRUG RESPONSE and adverse-reaction risk — e.g. clopidogrel non-response in CYP2C19 poor metabolizers, codeine toxicity in CYP2D6 ultrarapid metabolizers — not aging biology. Its value is medication safety and dosing, realized only when a relevant drug is prescribed.
- **Measurement.** Multi-gene genotyping by array, targeted sequencing, or PCR with copy-number assessment (important for CYP2D6); reported as star-allele diplotypes and metabolizer phenotypes. Performed once, results are lifelong. — *sample:* venous blood; *approx. cost:* $150-500.
- **Signal quality.** clean — Genotype is fixed and measured with high analytic accuracy; CYP2D6 copy-number and hybrid alleles are the main technical challenge. The clinical signal is conditional — it only matters for the specific drugs each gene affects.
- **Major confounders.** none for the genotype itself (fixed); phenoconversion: co-medications and inflammation can shift effective metabolizer status; drug-gene relevance is conditional on which medications are prescribed; ancestry-dependent star-allele frequencies
- **Testing cadence.** once in a lifetime (genotype is fixed); consult at each new prescription
- **Modifiability (fixed).** Genotype is fixed and cannot be changed. The actionable output is drug selection and dosing: genotype-guided prescribing (e.g. avoiding codeine/tramadol in CYP2D6 ultrarapid metabolizers, alternative antiplatelet therapy in CYP2C19 poor metabolizers, dose reduction of fluoropyrimidines in DPYD variant carriers) reduces clinically relevant adverse drug reactions. — *effect:* The PREPARE randomized trial found preemptive panel-based genotype-guided prescribing reduced clinically relevant adverse drug reactions by ~30% across 12 drug-gene pairs.; *timeframe:* immediate at the point of prescribing; *evidence:* RCT. Swen JJ et al. A 12-gene pharmacogenetic panel to prevent adverse drug reactions: an open-label, multicentre, controlled, cluster-randomised crossover implementation study (PREPARE / Ubiquitous Pharmacogenomics). Lancet 2023;401:347-356. doi:10.1016/S0140-6736(22)01841-4
- **Outcome evidence.** none recorded

### Polygenic risk scores (CAD, breast cancer, Alzheimer's)

*Also: PRS, polygenic score, genome-wide polygenic score, CAD PRS, PRS313 breast cancer, Alzheimer's PRS*

- **Mechanism.** A polygenic risk score sums the genome-wide allelic burden for a disease, weighting hundreds to millions of common variants by their GWAS effect sizes, into a single continuous score that is approximately normally distributed. Because it integrates lifelong genetic liability, a validated PRS stratifies disease risk independently of, and additively to, conventional risk factors. PRS is fixed from birth and is best read as a measure of inherited susceptibility for a specific outcome (coronary artery disease, breast cancer, Alzheimer's disease).
- **Measurement.** Genome-wide genotyping array (optionally imputed) or sequencing, then computation of a published, externally validated score; reported as a continuous score, typically standardized to SD units. Performed once, results are lifelong. — *sample:* venous blood; *approx. cost:* $100-400.
- **Signal quality.** moderate — Genotype inputs are fixed and accurate, but PRS performance is ancestry-dependent (most scores were derived in European-ancestry cohorts and attenuate in other populations) and varies with the score version and the outcome definition. The per-SD HR is well defined for a validated score in a matched population.
- **Major confounders.** ancestry / population stratification (scores derived mainly in European-ancestry cohorts); score version and SNP-selection method; outcome definition and ascertainment; age at assessment (effect size declines with older inclusion age)
- **Testing cadence.** once in a lifetime (genotype is fixed); score may be recomputed as improved versions are released
- **Modifiability (fixed).** The score itself is fixed (germline genotype). It is actionable indirectly: a high PRS justifies earlier and more intensive control of modifiable risk factors (e.g. lipid lowering and lifestyle for high CAD PRS, earlier/enhanced screening for high breast-cancer PRS). Trial and observational evidence shows adherence to a favorable lifestyle can roughly halve the absolute excess CAD risk conferred by a high PRS. — *effect:* PRS is immutable; downstream disease risk is partially modifiable — favorable lifestyle attenuates ~46% of the relative CAD risk in high-PRS individuals.; *timeframe:* not applicable for the score; years for downstream risk-factor modification; *evidence:* observational. Khera AV et al. Genetic Risk, Adherence to a Healthy Lifestyle, and Coronary Disease. N Engl J Med 2016;375:2349-2358. doi:10.1056/NEJMoa1605086
- **Outcome evidence.** All-cause mortality HR/SD 1.08 (tier B); CVD HR/SD 1.73 (tier A); Cancer HR/SD 1.61 (tier A); Dementia HR/SD 1.49 (tier B)

