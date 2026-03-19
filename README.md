# MoTrPAC Redox Scaling Analysis

**Mitochondrial oxidative capacity outpaces redox buffering during endurance training across rat tissues: sex-dimorphic scaling revealed by a thermodynamic framework analysis of MoTrPAC**

*Mislinski et al., in preparation*

---

## Summary

This repository contains the full analysis pipeline for a secondary analysis of the [MoTrPAC rat endurance training dataset](https://www.motrpac-data.org/) (Nature 2024; Cell Metabolism 2024) using the Vandiver–Neufer thermodynamic framework for mitochondrial redox biology.

We constructed an 85-gene set organized by thermodynamic function (ETS electron carriers, ATP synthase, NNT, glutathione system, thioredoxin system, superoxide dismutases, Q-pool feeders) and queried 42,770 observations across 19 tissues, 5 omic layers, 4 training time points, and both sexes. The central finding: after 8 weeks of endurance training, redox buffering capacity scales significantly sub-linearly with ETS expansion in both males (β = 0.55 conservative, p = 0.0017) and females (β = 0.12, p = 4.7 × 10⁻¹²), with sex-dimorphic quality control architectures — males coordinating ETS with damage-sensing programs (mtUPR, mitophagy), females with energy-sensing AMPK.

---

## Data source

All rat multi-omic data are accessed via the **MotrpacRatTraining6moData** Bioconductor package — no raw data download required.

```r
BiocManager::install("MotrpacRatTraining6moData")
```

For more information on the MoTrPAC dataset:
- Data portal: https://www.motrpac-data.org/
- Primary paper: MoTrPAC Study Group. *Nature* 629, 174–183 (2024)
- Multi-omic analysis: Amar et al. *Cell Metabolism* 36, 1411–1429 (2024)

Human data (script 07) uses supplementary tables from:
- Keshishian et al. bioRxiv 10.64898/2026.03.04.705181 — download S3.xlsx manually (see notes in `07_human_data.R`)

---

## Reproducing the analysis

### Prerequisites

- R ≥ 4.2
- Python ≥ 3.9 (for sensitivity table PDF generation only)

### Step-by-step

```bash
# 1. Clone the repository
git clone https://github.com/jackmis610/motrpac-redox.git
cd motrpac-redox

# 2. Install R dependencies (run once — takes ~5 min first time)
Rscript scripts/00_install_packages.R

# 3. Explore data structure (optional but recommended)
Rscript scripts/01_explore_data.R

# 4. Build the core query table — REQUIRED before steps 5–9
Rscript scripts/02_query_neufer.R
# → produces data/neufer_da_all.csv (~42,770 rows)

# 5. Generate tissue × category heatmap and BAT deep dive
Rscript scripts/03_figures.R

# 6. ETS vs. buffering regression + rigor analyses
Rscript scripts/04_sensitivity.R

# 7. Permutation test (~5–10 min)
Rscript scripts/05_permutation.R

# 8. Gene-by-gene correlation heatmap
Rscript scripts/06_gene_by_gene.R

# 9. Human MoTrPAC PASS1A cross-reference
Rscript scripts/07_human_data.R
```

All figures are written to `figures/` and all data outputs to `data/`.

---

## Directory structure

```
motrpac-redox/
├── scripts/
│   ├── 00_install_packages.R       # install dependencies
│   ├── 01_explore_data.R           # data structure exploration
│   ├── 02_query_neufer.R           # core query → data/neufer_da_all.csv
│   ├── 03_figures.R                # tissue heatmap + BAT figures
│   ├── 04_sensitivity.R            # regression + rigor analyses
│   ├── 05_permutation.R            # permutation test (R² specificity)
│   ├── 06_gene_by_gene.R           # 32×14 ETS×buffering correlation heatmap
│   ├── 07_human_data.R             # human MoTrPAC cross-reference
│   ├── neufer_geneset.R            # 85-gene set definition (sourced by 02)
│   ├── mito_dynamics_geneset.R     # 72-gene mito dynamics set
│   └── supp_*.R / supp_*.py        # supplementary analyses
│
├── data/
│   ├── neufer_da_all.csv           # 42,770-row primary query output
│   ├── supplementary_table1_gene_set.csv  # 85-gene set with metadata
│   ├── ets_buffering_correlation_data.csv
│   ├── regression_sensitivity_table.csv
│   ├── gene_gene_ets_buffering_cor_noBat.csv
│   ├── txnrd2_nnt_correlation.csv
│   └── human/                      # human MoTrPAC derived outputs
│
├── figures/                        # all generated figures (PNG)
├── CLAUDE.md                       # project context for AI-assisted development
└── README.md
```

---

## Dependencies

**R packages:**

| Package | Source | Purpose |
|---|---|---|
| `MotrpacRatTraining6moData` | Bioconductor | MoTrPAC data |
| `MotrpacRatTraining6mo` | Bioconductor | Analysis helpers |
| `ComplexHeatmap` | Bioconductor | Gene-by-gene heatmap |
| `circlize` | Bioconductor | Color scales |
| `tidyverse` | CRAN | Data wrangling + ggplot2 |
| `patchwork` | CRAN | Multi-panel figures |
| `ggrepel` | CRAN | Non-overlapping labels |
| `readxl` | CRAN | Human supplement XLSX |

**Python (optional, for supplementary table PDF):**
```
pip3 install reportlab Pillow
```

---

## Key results

| Finding | Value | Script |
|---|---|---|
| Male β (conservative, −BAT−BLOOD) | 0.551, R²=0.574, p=0.0017 | `04_sensitivity.R` |
| Female β (−BAT) | 0.115, R²=0.280, p=4.7×10⁻¹² | `04_sensitivity.R` |
| Framework R² vs. null (permutation) | p=0.006, curated R²=0.855 | `05_permutation.R` |
| Top positive pair (GPx4–Atp5f1b) | r=0.979 | `06_gene_by_gene.R` |
| NNT vs. ETS index (males, −BAT) | r=−0.397, p=0.11 | `supp_txnrd2_analysis.R` |
| Txnrd2 vs. ETS index (males, −BAT) | r=+0.719, p=0.001 | `supp_txnrd2_analysis.R` |

---

## Citation

If you use this analysis or code, please cite:

> Mislinski J, Subudhi AW, Jacobs RA. Mitochondrial oxidative capacity outpaces redox buffering during exercise training across rat tissues: sex-dimorphic scaling revealed by a thermodynamic framework analysis of MoTrPAC. *In preparation* (2026).

And the primary MoTrPAC papers:

> MoTrPAC Study Group. Temporal dynamics of the multi-omic response to endurance exercise training. *Nature* 629, 174–183 (2024).

> Amar D, et al. The mitochondrial multi-omic response to exercise training across rat tissues. *Cell Metabolism* 36, 1411–1429 (2024).

---

## Author

**Jack Mislinski**
Department of Human Physiology & Nutrition
University of Colorado Colorado Springs
Advisors: Dr. Andrew Subudhi, Dr. Robert Jacobs

---

## License

Code in this repository is released under the MIT License. The MoTrPAC dataset itself is governed by the [MoTrPAC data use policies](https://www.motrpac-data.org/).
