# Biomarker Profiles

Deliverable 1 (profiles) and Deliverable 3 (modifiability layer), generated from `data/biomarkers.json` by `tools/`.
Domains 2-4 of 14 are complete; remaining domains follow in subsequent passes.

Hazard ratios in the outcome line are standardized to **HR per +1 SD** where the
exposure is continuous; categorical exposures show the native HR (see `data/HR_STANDARDIZATION.md`).


## Lipids / Cardiovascular Risk

### ApoB

*Also: apolipoprotein B, apolipoprotein B-100, total apoB*

- **Mechanism.** ApoB is the structural protein of all atherogenic lipoproteins (LDL, VLDL, IDL, Lp(a)); each particle carries exactly one apoB molecule, so apoB concentration equals the number of atherogenic particles. Particle number drives arterial wall infiltration and atherosclerosis, making apoB a more direct causal proxy for atherogenic burden than cholesterol mass.
- **Measurement.** Immunoturbidimetric or immunonephelometric assay on a standard chemistry analyzer; non-fasting acceptable; standardized internationally. — *sample:* venous blood; *approx. cost:* $15-40.
- **Signal quality.** clean — Low biological and analytic variability; non-fasting measurement valid; well-standardized assays. More reproducible than calculated LDL-C.
- **Major confounders.** acute illness/inflammation; pregnancy; recent very high-fat meal (modest); lipid-lowering medication
- **Testing cadence.** baseline + annual; more frequent during lipid-lowering titration
- **Modifiability (high).** Statins, ezetimibe, and PCSK9 inhibitors lower apoB; PCSK9 inhibitors and high-intensity statins produce the largest reductions. — *effect:* High-intensity statin lowers apoB ~35-45%; adding ezetimibe +15-20%; PCSK9 inhibitors lower apoB 45-55%.; *timeframe:* 4-6 weeks to steady state; *evidence:* RCT meta-analysis. CTT Collaboration statin meta-analyses; FOURIER (Sabatine 2017) and ODYSSEY OUTCOMES (Schwartz 2018) PCSK9 RCTs.
- **Outcome evidence.** All-cause mortality HR/SD 1.11 (tier B); CVD HR/SD 1.43 (tier A); Cancer HR/SD 1.10 (tier C)

### LDL-C

*Also: LDL cholesterol, low-density lipoprotein cholesterol, calculated LDL*

- **Mechanism.** LDL-C measures the cholesterol mass carried by LDL particles, the dominant atherogenic lipoprotein. Cumulative lifetime LDL-C exposure causally drives atherosclerotic plaque formation; genetic, epidemiologic and trial evidence is concordant.
- **Measurement.** Calculated (Friedewald/Martin-Hopkins equation from total cholesterol, HDL-C, triglycerides) or direct enzymatic assay; fasting traditionally preferred for calculated values. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Calculated LDL-C is biased at high triglycerides or low LDL-C; cholesterol content per particle varies, causing discordance with particle number. Direct assays improve reliability.
- **Major confounders.** non-fasting state (affects calculated value); high triglycerides; acute illness; pregnancy; lipid-lowering drugs
- **Testing cadence.** baseline + annual; more frequent during treatment titration
- **Modifiability (high).** Statins (first-line), ezetimibe, PCSK9 inhibitors, bempedoic acid; diet (reduced saturated fat, fiber, plant sterols). — *effect:* High-intensity statin lowers LDL-C ~50%; ezetimibe ~15-20% additional; PCSK9 inhibitor ~50-60%; diet alone ~5-15%.; *timeframe:* 4-6 weeks to steady state; *evidence:* RCT meta-analysis. Cholesterol Treatment Trialists' (CTT) Collaboration meta-analyses, Lancet.
- **Outcome evidence.** All-cause mortality (documented, no pooled HR); CVD HR/SD 1.38 (tier A); Dementia (documented, no pooled HR)

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
- **Outcome evidence.** CVD HR/SD 1.37 (tier A)

### Triglyceride:HDL ratio

*Also: TG/HDL ratio, TG:HDL-C ratio, triglyceride to HDL cholesterol ratio*

- **Mechanism.** The TG/HDL-C ratio is a surrogate for insulin resistance and atherogenic dyslipidemia (high remnants, small dense LDL, low HDL). It captures the metabolic syndrome lipid pattern in a single number and predicts cardiovascular events and incident diabetes.
- **Measurement.** Calculated from standard lipid panel (triglycerides divided by HDL-C, both in mg/dL or both in mmol/L - units must match). — *sample:* venous blood; *approx. cost:* $10-25.
- **Signal quality.** noisy — Inherits high variability of triglycerides; ratio is unit-dependent and population-dependent (thresholds differ by ethnicity and mg/dL vs mmol/L).
- **Major confounders.** fasting state; alcohol; recent meal; insulin resistance status; ethnicity-specific thresholds; unit convention (mg/dL vs mmol/L)
- **Testing cadence.** baseline + periodic with lipid panel
- **Modifiability (moderate).** Weight loss, low-carbohydrate diet, exercise, glycemic control improve the ratio mainly by lowering triglycerides and raising HDL-C. — *effect:* Lifestyle interventions can lower the ratio substantially (often 25-50%) by improving insulin sensitivity.; *timeframe:* 8-12 weeks; *evidence:* observational. Low-carbohydrate diet and exercise intervention studies on TG/HDL-C.
- **Outcome evidence.** All-cause mortality native HR 2.00 [categorical, tier C]; CVD native HR 1.38 [categorical, tier C]; Dementia HR/SD 1.12 (tier C)

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
- **Outcome evidence.** CVD HR/SD 1.50 (tier A)

### Remnant cholesterol

*Also: remnant-C, remnant lipoprotein cholesterol, TRL-cholesterol*

- **Mechanism.** Remnant cholesterol is the cholesterol carried by triglyceride-rich lipoprotein remnants (VLDL, IDL, chylomicron remnants), calculated as total cholesterol minus LDL-C minus HDL-C. Remnant particles penetrate and are retained in the arterial wall and also drive inflammation; Mendelian randomization implicates remnant cholesterol as causal for ischemic heart disease.
- **Measurement.** Calculated from standard lipid panel (total cholesterol minus HDL-C minus LDL-C); can also be directly measured. Non-fasting samples capture remnant burden well. — *sample:* venous blood; *approx. cost:* $0-15 (derived from standard panel).
- **Signal quality.** moderate — Derived value; inherits triglyceride variability and depends on accuracy of the LDL-C estimate, but non-fasting measurement is a strength.
- **Major confounders.** fasting/non-fasting state; recent alcohol/meal; diabetes; method used to estimate LDL-C; acute illness
- **Testing cadence.** baseline + periodic with lipid panel
- **Modifiability (moderate).** Weight loss, reduced carbohydrate and alcohol intake, exercise, glycemic control; statins, fibrates, high-dose omega-3 and (investigationally) ANGPTL3/APOC3 inhibitors lower remnant cholesterol. — *effect:* Lifestyle can lower remnant-C substantially via triglyceride reduction; statins lower it ~20-30%.; *timeframe:* 8-12 weeks; *evidence:* observational. Lifestyle and lipid-lowering intervention studies; MR supports causal benefit of lowering.
- **Outcome evidence.** All-cause mortality HR/SD 1.17 (tier B); CVD HR/SD 1.67 (tier A)

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
- **Modifiability (fixed).** CAC reflects established calcified plaque and does not regress; the score itself cannot be lowered. Statins paradoxically increase calcified plaque density while stabilizing plaque. The actionable response is aggressive risk-factor and lipid management triggered by an elevated score. — *effect:* CAC score is essentially fixed/non-regressing; clinical value is in guiding intensity of preventive therapy.; *timeframe:* not applicable (score does not decrease); *evidence:* observational. Statin-CAC progression studies; guideline use of CAC to guide therapy intensity.
- **Outcome evidence.** All-cause mortality native HR 2.60 [categorical, tier A]; CVD native HR 6.84 [categorical, tier A]; Dementia native HR 1.71 [categorical, tier B]

### Carotid intima-media thickness (cIMT)

*Also: cIMT, carotid IMT, CCA-IMT, carotid intima-media thickness*

- **Mechanism.** cIMT measures the combined thickness of the intima and media layers of the carotid artery wall by ultrasound. Increased cIMT reflects subclinical atherosclerosis and arterial aging and is associated with future stroke and myocardial infarction, though it adds limited incremental predictive value over standard risk factors.
- **Measurement.** B-mode carotid ultrasound; standardized measurement of mean common carotid artery wall thickness; operator-dependent. — *sample:* imaging; *approx. cost:* $100-300.
- **Signal quality.** moderate — Operator- and protocol-dependent; measurement definitions vary (near vs far wall, segments included), limiting comparability. Reproducible within standardized protocols.
- **Major confounders.** age; operator/sonographer variability; measurement protocol differences (segment, wall); ultrasound equipment differences
- **Testing cadence.** one-time for risk assessment; serial measurement of progression has not proven clinically useful
- **Modifiability (moderate).** Statins, blood pressure control, and lifestyle changes slow or modestly reverse cIMT progression. — *effect:* Statins can slow progression and produce small regression of cIMT; absolute changes are small.; *timeframe:* 1-2 years for measurable change; *evidence:* RCT meta-analysis. Meta-analyses of statin and antihypertensive trials measuring cIMT progression (surrogate endpoint).
- **Outcome evidence.** CVD HR/SD 1.27 (tier A); Dementia HR/SD 1.08 (tier B)

### ApoE genotype

*Also: APOE genotype, apolipoprotein E genotype, APOE e4, ApoE4 carrier status*

- **Mechanism.** APOE has three common alleles (e2, e3, e4) encoding apolipoprotein E, which mediates clearance of triglyceride-rich and remnant lipoproteins. The e4 allele raises LDL-C and remnant levels and modestly increases coronary heart disease risk; it is also the strongest common genetic risk factor for late-onset Alzheimer's disease. e2 lowers LDL-C but can predispose to type III hyperlipoproteinemia.
- **Measurement.** Genotyping of the two APOE SNPs (rs429358, rs7412) by PCR or array; a single test, results are lifelong. — *sample:* venous blood; *approx. cost:* $50-200.
- **Signal quality.** clean — Genotype is fixed and measured with essentially no biological or analytic variability; the test is definitive once performed.
- **Major confounders.** none for the genotype itself (fixed); phenotypic expression modified by environment, sex, ancestry
- **Testing cadence.** once in a lifetime (genotype is fixed)
- **Modifiability (fixed).** Genotype is fixed and cannot be changed. e4 carriers benefit from earlier and more intensive management of modifiable risk factors (LDL-C lowering, blood pressure, exercise, sleep) for both cardiovascular and dementia risk. — *effect:* Not applicable - genotype is immutable; downstream lipid and cognitive risk is partially modifiable.; *timeframe:* not applicable; *evidence:* mechanistic. APOE biology; risk-factor modification trials in e4 carriers.
- **Outcome evidence.** All-cause mortality native HR 1.22 [categorical, tier B]; CVD native HR 1.09 [categorical, tier B]; Dementia native HR 3.20 [categorical, tier A]


## Mitochondrial / Cardiorespiratory Fitness

### VO2max (relative)

*Also: maximal oxygen uptake, cardiorespiratory fitness, peak VO2, CRF, aerobic capacity*

- **Mechanism.** VO2max integrates the capacity of the lungs, heart, vasculature and skeletal-muscle mitochondria to deliver and utilize oxygen. It is the single most robust functional predictor of all-cause and cardiovascular mortality, and declining VO2max tracks the loss of physiologic reserve that defines biological aging.
- **Measurement.** Gold standard: graded maximal exercise test (treadmill or cycle ergometer) with expired-gas analysis (cardiopulmonary exercise testing, CPET). Submaximal estimates and treadmill-time-derived METs are common surrogates in large cohorts. — *sample:* functional test; *approx. cost:* $150-400.
- **Signal quality.** clean — Directly measured CPET VO2max is highly reproducible (test-retest CV ~3-5%). Estimated/submaximal values introduce moderate noise; effort dependence and ceiling at true max are the main measurement issues.
- **Major confounders.** age; sex; body composition / adiposity; submaximal vs maximal effort; altitude; beta-blocker use; test modality (treadmill vs cycle); reverse causation from subclinical disease
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (high).** Structured aerobic / endurance training, especially moderate-to-vigorous continuous exercise and high-intensity interval training. — *effect:* Typically 10-25% increase in VO2max with 3-6 months of training; larger relative gains in deconditioned and older individuals, attenuated in already-fit and very old.; *timeframe:* 8-12 weeks for measurable gains; 6 months for substantial improvement.; *evidence:* RCT meta-analysis. Numerous RCT meta-analyses of aerobic and HIIT training (e.g., Bacon 2013 PLoS One; Milanovic 2015 Sports Med) consistently show ~10-25% VO2max improvement.
- **Outcome evidence.** All-cause mortality HR/SD 0.76 (tier A); CVD HR/SD 0.72 (tier A); Cancer HR/SD 0.64 (tier B); Dementia HR/SD 0.81 (tier B)

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
- **Modifiability (high).** Progressive resistance training, especially high-velocity / power-oriented training and explosive sit-to-stand exercise. — *effect:* Lower-body strength improves 25-100%+ and muscle power improves substantially in older adults with 8-12 weeks of resistance/power training; power-oriented training is more effective than slow heavy training for power gains.; *timeframe:* 8-12 weeks.; *evidence:* RCT meta-analysis. Power- and resistance-training RCT meta-analyses in older adults.
- **Outcome evidence.** All-cause mortality native HR 1.57 [categorical, tier B]; Frailty (documented, no pooled HR)

### Muscle mass (appendicular lean mass index)

*Also: appendicular lean mass, ALM, ALMI, skeletal muscle mass index, lean body mass, sarcopenia*

- **Mechanism.** Appendicular lean mass index (limb lean mass / height^2) quantifies skeletal-muscle quantity, the substrate for strength, metabolic glucose disposal and protein reserve. Low muscle mass (sarcopenia) reflects anabolic decline and predicts disability and mortality, though strength predicts outcomes more strongly than mass alone.
- **Measurement.** Dual-energy X-ray absorptiometry (DEXA) appendicular lean mass; bioelectrical impedance is a lower-cost surrogate. — *sample:* imaging; *approx. cost:* $75-250.
- **Signal quality.** moderate — DEXA lean mass is reproducible (CV ~1-2%), but hydration status affects readings and lean mass includes non-muscle tissue; muscle quality/strength is not captured.
- **Major confounders.** body size / height normalization; hydration status; fat mass (obesity); DEXA vs BIA method; age; sex; ethnicity
- **Testing cadence.** baseline + every 1-2 years
- **Modifiability (moderate).** Progressive resistance training, with adequate dietary protein; combined with anabolic stimulus in deficient states. — *effect:* Resistance training increases appendicular lean mass by roughly 0.5-2 kg over 3-6 months in older adults; gains are modest in magnitude relative to strength gains.; *timeframe:* 12-24 weeks.; *evidence:* RCT meta-analysis. Resistance-training RCT meta-analyses (e.g., Peterson 2011 Ageing Res Rev) show significant lean-mass gains in older adults.
- **Outcome evidence.** All-cause mortality HR/SD 0.50 (tier B); Frailty (documented, no pooled HR)

### Gait speed

*Also: walking speed, usual gait speed, habitual gait speed, 4-meter walk*

- **Mechanism.** Usual gait speed is an integrative 'vital sign' that summarizes the function of the cardiovascular, musculoskeletal, neurological and energetic systems. It declines with age, predicts survival and disability across populations, and is one of the most validated single functional measures in geriatrics.
- **Measurement.** Time to walk a measured course (typically 4 m or 6 m) at usual pace, often the better of two trials; can also be measured over longer distances or with instrumented walkways. — *sample:* functional test; *approx. cost:* $0-50.
- **Signal quality.** clean — Highly reproducible and standardized; main variability is course length, static vs dynamic start and usual-pace instruction, all controllable.
- **Major confounders.** age; height / leg length; orthopedic and neurological conditions; course length and start protocol; footwear; test environment
- **Testing cadence.** baseline + annual
- **Modifiability (moderate).** Multicomponent exercise (aerobic + resistance + balance/gait training); resistance and power training in deconditioned older adults. — *effect:* Exercise interventions improve gait speed by roughly 0.05-0.10+ m/s in older adults; ~0.05 m/s is a clinically meaningful change.; *timeframe:* 12-24 weeks.; *evidence:* RCT meta-analysis. Exercise-intervention RCT meta-analyses in older adults show significant gait-speed improvement.
- **Outcome evidence.** All-cause mortality HR/SD 0.71 (tier A); Dementia HR/SD 1.59 (tier B); Frailty (documented, no pooled HR)

### Sit-to-stand (chair-rise test)

*Also: chair-rise test, five-times sit-to-stand, 5xSTS, 30-second chair stand, sitting-rising test*

- **Mechanism.** The chair-rise / sit-to-stand test measures lower-body strength, power and balance during a functionally essential task. Poor performance signals loss of mobility reserve, fall risk and frailty, and predicts mortality independent of conventional risk factors.
- **Measurement.** Time for five chair rises (5xSTS), number of rises in 30 seconds, or the floor-based sitting-rising test scored 0-10. — *sample:* functional test; *approx. cost:* $0-25.
- **Signal quality.** moderate — Reproducible with standardized protocol; depends on chair height, arm use, technique and motivation, and different versions (timed, count, floor-based) are not interchangeable.
- **Major confounders.** chair height / arm use; joint pain and orthopedic limitation; test version; motivation / effort; age; balance impairment
- **Testing cadence.** baseline + annual
- **Modifiability (high).** Progressive resistance and power training targeting the lower body; practiced sit-to-stand exercise. — *effect:* Chair-rise time and repetitions improve substantially with 8-12 weeks of resistance/power training in older adults.; *timeframe:* 8-12 weeks.; *evidence:* RCT meta-analysis. Resistance- and power-training RCT meta-analyses in older adults.
- **Outcome evidence.** All-cause mortality native HR 3.84 [categorical, tier B]; CVD native HR 6.05 [categorical, tier C]; Frailty (documented, no pooled HR)

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
- **Outcome evidence.** All-cause mortality HR/SD 1.11 (tier A); CVD HR/SD 1.16 (tier B)

### Heart rate variability (HRV)

*Also: HRV, SDNN, RMSSD, vagal tone, autonomic balance*

- **Mechanism.** Heart rate variability quantifies beat-to-beat variation in heart rate, indexing autonomic (mainly parasympathetic) regulation. Higher HRV reflects healthy vagal tone and physiologic adaptability; HRV declines with age and disease and lower values predict mortality and cardiovascular events.
- **Measurement.** Time-domain (SDNN, RMSSD) and frequency-domain measures from ECG or photoplethysmography over short recordings (~5 min), 24-hour Holter, or consumer wearables. — *sample:* wearable; *approx. cost:* $0-300.
- **Signal quality.** noisy — HRV is highly state-dependent (posture, breathing, stress, sleep, caffeine, time of day) with large day-to-day variability; recording length and device differ across studies. Trends are more informative than single values.
- **Major confounders.** age; posture and breathing rate; physical activity and stress before recording; recording length; device / measurement method; medications (beta-blockers); alcohol and caffeine; time of day / circadian phase
- **Testing cadence.** trend monitoring (e.g., daily wearable) or baseline + annual clinical recording
- **Modifiability (moderate).** Aerobic exercise training; secondarily slow-paced breathing / HRV-biofeedback and improved sleep and stress management. — *effect:* Aerobic training produces modest increases in HRV (RMSSD/SDNN); biofeedback can acutely raise HRV but durable resting changes are smaller.; *timeframe:* 8-16 weeks.; *evidence:* RCT meta-analysis. Meta-analyses of exercise training and of HRV-biofeedback show modest HRV improvements.
- **Outcome evidence.** All-cause mortality HR/SD 1.19 (tier B); CVD HR/SD 1.21 (tier B)


## Inflammation

### hsCRP

*Also: high-sensitivity C-reactive protein, CRP, hs-CRP*

- **Mechanism.** Acute-phase protein synthesized by the liver in response to IL-6; serves as a downstream integrator of systemic low-grade inflammation ('inflammaging'). Elevated hsCRP marks the chronic inflammatory state that accompanies atherosclerosis, metabolic dysfunction, and biological aging, though Mendelian randomization indicates it is a marker rather than a cause of disease.
- **Measurement.** High-sensitivity immunoturbidimetric or immunonephelometric assay on a standard clinical chemistry analyzer. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Analytically robust and well-standardized, but high within-person biological variability: levels spike acutely with infection, trauma, or recent exercise. Single measurements are noisy; serial measurements or values >10 mg/L should be excluded as acute reactions. Tertile/quartile stability over years is moderate.
- **Major confounders.** acute infection or recent illness; obesity and adiposity; smoking; estrogen / hormone therapy; recent vigorous exercise or trauma; chronic conditions (RA, periodontal disease); statin and NSAID use
- **Testing cadence.** baseline + annual; repeat if >10 mg/L to exclude acute inflammation
- **Modifiability (high).** Statin therapy, weight loss, and regular aerobic exercise; canakinumab/colchicine for residual inflammatory risk. — *effect:* Rosuvastatin lowered hsCRP ~37% in JUPITER; diet-induced weight loss produces a comparable reduction (roughly proportional to kg lost); regular exercise lowers hsCRP by ~20-30%.; *timeframe:* Weeks for statins; 3-6 months for weight loss and exercise.; *evidence:* RCT. Ridker PM et al. JUPITER trial. N Engl J Med. 2008;359:2195-2207. POUNDS LOST trial (weight loss & hsCRP), Am J Clin Nutr 2012.
- **Outcome evidence.** All-cause mortality HR/SD 1.54 (tier A); CVD HR/SD 1.37 (tier A); Cancer native HR 1.25 [categorical, tier B]; Dementia HR/SD 1.37 (tier B); Frailty HR/SD 1.27 (tier C)

### IL-6

*Also: interleukin-6, IL6*

- **Mechanism.** Pleiotropic pro-inflammatory cytokine and the principal driver of hepatic acute-phase protein synthesis (including CRP and fibrinogen). A central mediator of 'inflammaging'; unlike CRP, human genetic (IL6R) evidence supports a causal role in atherosclerotic disease, making the IL-6 pathway a validated therapeutic target.
- **Measurement.** ELISA or high-sensitivity immunoassay (e.g., electrochemiluminescence); proteomic panels (Olink, SomaScan) also quantify it. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** moderate — Biologically meaningful but analytically more variable than CRP: short half-life, diurnal variation, sensitivity to recent exercise/infection, and assay-platform differences. Single measurements are noisy; cross-assay comparability is limited.
- **Major confounders.** acute infection / illness; obesity (adipose tissue secretes IL-6); recent exercise; smoking; diurnal rhythm; assay platform variability; chronic inflammatory disease
- **Testing cadence.** baseline + annual; not yet routine in standard care
- **Modifiability (high).** Direct pathway blockade (IL-6R antibody tocilizumab; IL-1beta antibody canakinumab; colchicine) lowers IL-6 axis activity; lifestyle (exercise, weight loss) produces smaller reductions. — *effect:* Canakinumab (CANTOS) reduced hsCRP ~37% and IL-6 substantially with a 15% cut in major cardiovascular events; tocilizumab strongly suppresses CRP/fibrinogen. Exercise and weight loss lower IL-6 modestly (~10-25%).; *timeframe:* Days to weeks for biologics; months for lifestyle change.; *evidence:* RCT. Ridker PM et al. CANTOS trial (canakinumab). N Engl J Med. 2017;377:1119-1131. IL6R MR consortium, Lancet 2012.
- **Outcome evidence.** All-cause mortality HR/SD 1.48 (tier A); CVD HR/SD 1.25 (tier A); Cancer HR/SD 1.20 (tier C); Dementia HR/SD 1.32 (tier B); Frailty native HR 1.19 [unconvertible, tier C]

### TNF-alpha

*Also: tumor necrosis factor alpha, TNF-a, TNFa*

- **Mechanism.** Pro-inflammatory cytokine central to innate immunity, apoptosis signaling, insulin resistance, and cachexia. Implicated in inflammaging and tissue catabolism; its soluble receptors (sTNFR1/sTNFR2) are more stable, prognostically informative analytes than free TNF-alpha itself.
- **Measurement.** High-sensitivity ELISA or multiplex immunoassay; soluble TNF receptors (sTNFR1/2) often preferred for stability. — *sample:* venous blood; *approx. cost:* $30-90.
- **Signal quality.** noisy — Free circulating TNF-alpha is present at very low concentrations with a short half-life, marked diurnal variation, and poor assay reproducibility near the detection limit. Soluble receptors are far more reliable; most robust prognostic data come from sTNFR measurements rather than TNF-alpha itself.
- **Major confounders.** assay sensitivity limits and platform differences; diurnal variation; acute infection; obesity / adipose secretion; renal function (affects soluble receptor clearance); smoking
- **Testing cadence.** research / specialist use; not routine
- **Modifiability (moderate).** Anti-TNF biologics (etanercept, infliximab, adalimumab) in inflammatory disease; lifestyle (weight loss, exercise) for modest reductions in the general population. — *effect:* Anti-TNF therapy strongly suppresses TNF signaling and, in observational meta-analysis of rheumatoid arthritis cohorts, was associated with reduced cardiovascular events (RR ~0.46); lifestyle change yields modest reductions.; *timeframe:* Days to weeks for biologics; months for lifestyle.; *evidence:* observational. Barnabe C et al. Systematic review and meta-analysis: anti-TNF therapy and cardiovascular events in rheumatoid arthritis. Arthritis Care Res. 2011;63(4):522-529. PMID: 20957658.
- **Outcome evidence.** All-cause mortality HR/SD 1.34 (tier C); CVD HR/SD 1.63 (tier C); Frailty (documented, no pooled HR)

### GlycA

*Also: glycoprotein acetyls, NMR glycoprotein acetylation signal, N-acetyl glycoprotein*

- **Mechanism.** Composite NMR signal arising from N-acetyl methyl groups of glycosylated acute-phase proteins (alpha-1-acid glycoprotein, haptoglobin, alpha-1-antitrypsin, alpha-1-antichymotrypsin, transferrin). Integrates both protein abundance and glycosylation state, capturing a stable, time-averaged measure of chronic systemic inflammation that predicts long-term disease risk.
- **Measurement.** Proton (1H) nuclear magnetic resonance spectroscopy of plasma/serum, typically as part of a multi-analyte NMR metabolomics panel. — *sample:* venous blood; *approx. cost:* $30-100.
- **Signal quality.** clean — Lower within-person biological variability than CRP because it averages several glycoproteins with longer half-lives; NMR quantification is highly reproducible. Less sensitive to transient acute spikes, making it a comparatively stable inflammation index. Newer marker with less standardization across labs.
- **Major confounders.** obesity / adiposity; smoking; acute infection (less than CRP); metabolic syndrome and diabetes; NMR platform / lab differences; statin use
- **Testing cadence.** baseline + every few years; not yet routine clinical use
- **Modifiability (moderate).** Weight loss, exercise, and statin therapy lower GlycA, paralleling effects on other acute-phase markers. — *effect:* Statins and lifestyle interventions reduce GlycA modestly (roughly 5-15%); magnitude tracks reductions in adiposity and CRP.; *timeframe:* Months.; *evidence:* observational. Observational and secondary-analysis data (e.g., statin trials with NMR substudies); no dedicated large RCT with GlycA as primary endpoint.
- **Outcome evidence.** All-cause mortality HR/SD 1.24 (tier B); CVD HR/SD 1.22 (tier B); Cancer HR/SD 1.08 (tier C); Frailty (documented, no pooled HR)

### Neutrophil:lymphocyte ratio (NLR)

*Also: NLR, neutrophil-to-lymphocyte ratio*

- **Mechanism.** Ratio derived from a standard complete blood count; integrates the innate (neutrophil) and adaptive (lymphocyte) immune compartments. A high ratio reflects active inflammatory/stress responses with relative lymphopenia, a pattern associated with immunosenescence, physiologic stress, and adverse aging trajectories.
- **Measurement.** Calculated from an automated complete blood count with differential (absolute neutrophil count divided by absolute lymphocyte count). — *sample:* venous blood; *approx. cost:* $5-20.
- **Signal quality.** noisy — Cheap and ubiquitous (derived from routine CBC), but highly non-specific and labile: shifts acutely with infection, physical/psychological stress, corticosteroids, and time of day. No standardized cut-points; tertile/quartile definitions vary widely across studies, limiting comparability.
- **Major confounders.** acute infection or illness; corticosteroid and other immunomodulatory drugs; physical / psychological stress; diurnal variation; smoking; absence of standardized thresholds; hematologic conditions
- **Testing cadence.** opportunistic -- available on any CBC; interpret with clinical context
- **Modifiability (low).** No NLR-specific intervention; it falls with resolution of underlying inflammation/illness and improves modestly with exercise and weight loss. — *effect:* Not well quantified; NLR normalizes when an acute or chronic inflammatory driver is treated.; *timeframe:* Days to weeks (acute resolution) to months (lifestyle).; *evidence:* observational. No dedicated RCT targeting NLR as an endpoint; inferred from CBC normalization in observational data.
- **Outcome evidence.** All-cause mortality HR/SD 1.05 (tier B); CVD native HR 1.62 [categorical, tier B]; Cancer native HR 1.27 [categorical, tier C]; Dementia HR/SD 1.16 (tier B)

### Fibrinogen

*Also: plasma fibrinogen, factor I*

- **Mechanism.** Hepatically synthesized acute-phase glycoprotein and the terminal substrate of the coagulation cascade. Sits at the intersection of inflammation and thrombosis: elevated levels increase blood viscosity, platelet aggregation, and clot formation, mechanistically linking systemic inflammation to atherothrombotic events.
- **Measurement.** Clauss clotting assay (functional) or immunoassay (antigen) on a coagulation analyzer. — *sample:* venous blood; *approx. cost:* $10-30.
- **Signal quality.** moderate — Reasonably standardized assay but, as an acute-phase reactant, levels rise with infection, smoking, pregnancy, and inflammation. Within-person variability is moderate; longer half-life than CRP gives somewhat more stable readings.
- **Major confounders.** acute infection / inflammation; smoking; pregnancy and estrogen / hormone therapy; obesity; diabetes; assay method (Clauss vs antigen); age
- **Testing cadence.** baseline + periodic; available within coagulation panels
- **Modifiability (moderate).** Smoking cessation, weight loss, and regular physical activity lower fibrinogen; no fibrinogen-specific drug therapy in primary prevention. — *effect:* Smoking cessation and exercise reduce fibrinogen modestly (roughly 0.2-0.5 g/L over months); fibrate drugs lower it somewhat but without proven outcome benefit attributable to fibrinogen.; *timeframe:* Months.; *evidence:* observational. Observational and secondary trial data; no RCT demonstrates outcome benefit from fibrinogen lowering per se.
- **Outcome evidence.** All-cause mortality HR/SD 1.70 (tier A); CVD HR/SD 1.94 (tier A); Dementia HR/SD 1.30 (tier B); Frailty HR/SD 1.31 (tier C)

### Homocysteine

*Also: total homocysteine, tHcy, hyperhomocysteinemia*

- **Mechanism.** Sulfur-containing amino acid intermediate of methionine metabolism. Elevated levels reflect impaired one-carbon metabolism (folate/B12/B6 status, renal function, MTHFR genotype) and have been linked to endothelial dysfunction, oxidative stress, and neurotoxicity. Often grouped with inflammatory/vascular risk markers though it is not a classic acute-phase reactant.
- **Measurement.** Immunoassay or LC-MS/MS on fasting plasma; sample must be processed promptly as levels rise with delayed serum separation. — *sample:* venous blood; *approx. cost:* $20-60.
- **Signal quality.** moderate — Analytically reliable but strongly determined by modifiable nutritional status (folate, B12, B6) and renal function; pre-analytic handling matters (levels rise if red cells are not separated quickly). Genetically influenced (MTHFR). Less labile day-to-day than cytokines.
- **Major confounders.** folate / vitamin B12 / B6 status; renal impairment (major determinant); MTHFR C677T genotype; age and male sex; hypothyroidism; delayed sample processing; certain drugs (methotrexate, antiepileptics)
- **Testing cadence.** baseline; recheck after B-vitamin repletion if elevated
- **Modifiability (high).** Folic acid plus vitamin B12 (and B6) supplementation reliably lowers plasma homocysteine. — *effect:* B-vitamin supplementation lowers homocysteine by roughly 25-30%; however this lowering did NOT translate into reduced cardiovascular events or mortality in large RCTs.; *timeframe:* Weeks for the biochemical change.; *evidence:* RCT meta-analysis. Marti-Carvajal AJ et al. Cochrane Database Syst Rev. 2017;8:CD006612. JAMA Intern Med 2010 meta-analysis of 8 RCTs.
- **Outcome evidence.** All-cause mortality HR/SD 1.26 (tier B); CVD native HR 0.98 [categorical, tier A]; Dementia HR/SD 1.12 (tier B)

### Ferritin (as inflammation marker)

*Also: serum ferritin, hyperferritinemia*

- **Mechanism.** Primary intracellular iron-storage protein; circulating ferritin reflects body iron stores but is ALSO a positive acute-phase reactant that rises with inflammation independent of iron. This dual identity makes elevated ferritin a mixed signal of iron overload, hepatocellular injury, and/or systemic inflammation.
- **Measurement.** Immunoassay (chemiluminescent or immunoturbidimetric) on a standard analyzer. — *sample:* venous blood; *approx. cost:* $15-40.
- **Signal quality.** noisy — Interpretation is genuinely ambiguous: an elevated value may indicate iron overload, an acute-phase response, liver disease, malignancy, or metabolic syndrome. Without concurrent CRP, transferrin saturation, and clinical context, a single ferritin is difficult to interpret as an inflammation marker specifically.
- **Major confounders.** iron stores / iron overload (primary determinant); acute-phase response / inflammation; liver disease and alcohol; malignancy; metabolic syndrome and fatty liver; recent transfusion; age and sex
- **Testing cadence.** baseline; interpret alongside CRP and transferrin saturation
- **Modifiability (low).** Phlebotomy / venesection for genuine iron overload; treatment of the underlying inflammatory or hepatic condition when ferritin elevation is reactive. — *effect:* Phlebotomy reliably lowers ferritin when iron-driven; reactive (inflammatory) hyperferritinemia falls only when the inflammatory cause resolves.; *timeframe:* Weeks to months.; *evidence:* mechanistic. No outcome RCT shows mortality benefit from lowering ferritin per se in the general population; iron-reduction trials (e.g., VA cooperative trial in PAD) were largely null.
- **Outcome evidence.** All-cause mortality HR/SD 1.05 (tier B); CVD (documented, no pooled HR); Cancer (documented, no pooled HR)

### suPAR

*Also: soluble urokinase plasminogen activator receptor, suPAR*

- **Mechanism.** Circulating cleaved form of the membrane-bound urokinase plasminogen activator receptor (uPAR), shed mainly from immune and endothelial cells. Reflects chronic immune activation; mechanistically implicated in atherosclerosis and as a pathogenic mediator of kidney injury. Marketed as a stable, broad 'biomarker of immune activation and aging'.
- **Measurement.** ELISA (suPARnostic) or proteomic immunoassay; quantified in plasma or serum. — *sample:* venous blood; *approx. cost:* $30-100.
- **Signal quality.** clean — Comparatively stable analyte with lower within-person and diurnal variability than CRP or IL-6, and less reactive to acute infection -- a key reason it is promoted as a marker of chronic inflammation and biological aging. Main limitations are cost and limited assay standardization / availability.
- **Major confounders.** renal function (suPAR rises as GFR falls); smoking; age; obesity; HIV and chronic infection; chronic inflammatory disease
- **Testing cadence.** baseline; specialist / research use, not yet routine
- **Modifiability (low).** No suPAR-specific therapy; levels fall with smoking cessation and treatment of underlying chronic inflammation/infection. Experimental anti-uPAR strategies exist for chronic kidney disease. — *effect:* Smoking cessation lowers suPAR modestly; no intervention has demonstrated outcome benefit via suPAR lowering.; *timeframe:* Months.; *evidence:* mechanistic. No RCT targets suPAR as a modifiable endpoint in the general population; modifiability inferred from observational determinants.
- **Outcome evidence.** All-cause mortality HR/SD 1.43 (tier B); CVD HR/SD 1.35 (tier B); Cancer (documented, no pooled HR)

