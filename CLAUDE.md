# MoTrPAC Redox Analysis Project

## What this is
Secondary analysis of MoTrPAC rat endurance training data (Nature 2024, Cell Metabolism 2024)
through the lens of the Neufer mitochondrial redox framework.

## Core hypothesis
Exercise training enhances not just oxidative capacity (ETS expansion) but also the
mitochondrial "electrical grid" — NNT-driven NADPH generation that maintains proteome
redox tone via glutathione and thioredoxin buffering circuits. This can be tested by
querying MoTrPAC multi-omic data for coordinated temporal responses across:
1. ETS components (ΔGredox)
2. ATP synthase (ΔGpmf → ΔGATP)
3. NNT (ΔGpmf → ΔGNADPH)
4. Glutathione/thioredoxin systems (ΔGNADPH → proteome redox tone)
5. Q-pool feeder enzymes (reductive stress entry points)

## Script execution order
```
00_install_packages.R    — install dependencies (run once)
01_explore_data.R        — data structure exploration
02_query_neufer.R        — build data/neufer_da_all.csv (required by 03–07)
03_figures.R             — tissue heatmap + BAT deep dive
04_sensitivity.R         — GPx4 specificity, trx-prot concordance, ETS vs buffering regression
05_permutation.R         — permutation test (R² specificity, ~5–10 min)
06_gene_by_gene.R        — 32×14 ETS×buffering correlation heatmap
07_human_data.R          — human MoTrPAC PASS1A cross-reference
```

Supplementary analyses (supp_*.R / supp_*.py) extend specific results.

## Key files
- `scripts/neufer_geneset.R` — 85-gene set (15 categories): ETS, NNT, GSH, TRX, SOD, Q-pool, etc.
- `scripts/mito_dynamics_geneset.R` — 72-gene set (16 categories) for mito dynamics analysis
- `data/neufer_da_all.csv` — 42,770 rows: all Neufer genes × all DA tables (output of 02)
- `data/supplementary_table1_gene_set.csv` — 85-gene set with full metadata

## Data access
MoTrPAC R packages:
- `MotrpacRatTraining6moData` — differential analysis results + normalized data
- `MotrpacRatTraining6mo` — analysis functions

Data objects follow: `{ASSAY}_{TISSUE}_DA` naming convention.
Assays: TRNSCRPT, PROT, PHOSPHO, ACETYL, UBIQ, METAB, ATAC, METHYL, IMMUNO
Tissues: SKMGN (not SKM_GN), HEART, LIVER, BAT, KIDNEY, LUNG, etc.

## Key DA columns
- TRNSCRPT: logFC, shrunk_logFC, zscore, adj_p_value, p_value, selection_fdr
- PROT/PHOSPHO/ACETYL/UBIQ: logFC, tscore, adj_p_value, p_value, selection_fdr
- shrunk_logFC only exists for TRNSCRPT

## Key findings
- 84/85 Neufer genes map (Prodh missing from FEATURE_TO_GENE)
- NNT: 0 significant hits at adj_p < 0.05 in any tissue/assay at 8w
- ETS vs buffering regression: sub-linear in both sexes (β_males_conservative = 0.55,
  β_females = 0.12); conservative estimate excludes BAT (Cook's D = 16.7) and BLOOD
- BAT artifact: multi-feature mapping in FEATURE_TO_GENE produces negative ETS mean;
  treat BAT as sensitivity case throughout
- Txnrd2: positively correlated with ETS index (r = +0.72, p = 0.001, BAT excluded, males)
- NNT: negative cross-tissue trend (r = -0.40, p = 0.11, BAT excluded, males)
- Permutation R² p = 0.006 (curated framework exceeds 99.4% of random gene sets)

## R environment
Use `Rscript` to execute scripts. All packages installed via 00_install_packages.R.
PDF output: use `device = cairo_pdf` (handles Unicode characters).
