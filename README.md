# Air Pollution in Indian Cities (2022–2025): A Statistical Study Using R

[![R Version](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Report-Complete-brightgreen.svg)](reports/lab_report.md)

This repository contains the complete laboratory investigation, statistical codebase, outputs, and master report for **"What Drives Air Pollution in Indian Cities? A Statistical Study Using R"**.

The study evaluates multi-year atmospheric dynamics and particulate matter ($PM_{2.5}$) across **29 major Indian cities** (balanced panel of 842,160 hourly records / 35,090 daily observations) from **August 2022 to November 2025**.

---

## 📄 Full Master Report
The complete, publication-grade lab report with all mathematical models, hypotheses, 21 embedded figures, detailed tables, and complete R scripts is available at:
👉 **[reports/lab_report.md](reports/lab_report.md)**

---

## 🔬 Experiments Conducted

| # | Script | Experiment Title | Tested Hypothesis / Focus | Key Finding |
|---|---|---|---|---|
| **1** | [`scripts/01_import_and_clean.R`](scripts/01_import_and_clean.R) | **Data Ingestion & Hygiene** | Classification, missingness audit, daily aggregation | Dropped 8 unpopulated upper-air columns; generated clean daily panel |
| **2** | [`scripts/02_summary_statistics.R`](scripts/02_summary_statistics.R) | **Summary Statistics & Dispersion** | Central tendency, variance, skewness, kurtosis | National mean $PM_{2.5} = 34.74\ \mu\text{g/m}^3$ (skewness = $+2.16$); Gurugram & Delhi worst |
| **3** | [`scripts/03_correlation_simple_regression.R`](scripts/03_correlation_simple_regression.R) | **Correlation & Simple Regression** | $H_0: \rho = 0$ for $PM_{2.5}$ vs $PM_{10}$, wind, humidity | $PM_{2.5}$ strongly linked with $CO$ ($r = 0.765$) and $PM_{10}$ ($r = 0.678, \rho = 0.951$) |
| **4** | [`scripts/04_multiple_linear_regression.R`](scripts/04_multiple_linear_regression.R) | **Multiple Linear Regression** | Predictors of $PM_{2.5}$; first-principles VIF diagnostics | Combined $R^2_{\text{adj}} = 0.6868$; $CO$ is dominant driver ($\beta^* = +0.576$); all VIF $< 4.4$ |
| **5** | [`scripts/05_distribution_fitting.R`](scripts/05_distribution_fitting.R) | **Theoretical Distribution Fitting** | Normal vs Log-Normal MLE; Poisson exceedance | Log-Normal decisively beats Normal ($\Delta\text{AIC} = -23,536$); Poisson overdispersed ($D = 14.93$) |
| **6** | [`scripts/06_one_sample_tests.R`](scripts/06_one_sample_tests.R) | **One-Sample Hypothesis Tests** | City means vs NAAQS ($60\ \mu\text{g/m}^3$) and WHO ($15\ \mu\text{g/m}^3$) | National mean exceeds WHO ($t = 131.1, p < 10^{-15}$); Gurugram & Delhi fail NAAQS on $>50\%$ days |
| **7** | [`scripts/07_ttests.R`](scripts/07_ttests.R) | **Two-Sample Hypothesis Tests** | Rain effect, paired festival spike, weekend effect | Rain scrubs PM2.5 by 47% ($d = 0.833$); Festival surge $+23.3\%$ ($p = 0.00037$); No weekend effect |
| **8** | [`scripts/08_chisquare_tests.R`](scripts/08_chisquare_tests.R) | **Chi-Square Tests** | Goodness-of-Fit across seasons; Contingency tests | Winter & Post-Monsoon $= 79\%$ of bad days; Crop burning doubles exceedance odds ($\text{OR} = 2.04$) |
| **9** | [`scripts/09_anova_crd_rbd.R`](scripts/09_anova_crd_rbd.R) | **ANOVA (CRD & RBD Designs)** | CRD One-Way vs RBD Two-Way blocking by City | **Relative Efficiency (RE) = 3.698**: Blocking by city yields a 270% gain in precision |

---

## 📁 Repository Structure

```
AQI_report/
├── README.md                            # Project overview & quickstart
├── .gitignore                           # Excludes >100MB raw files and temporary data
├── datasets/
│   └── README.md                        # Provenance and instructions to download raw Kaggle dataset
├── data/
│   └── aqi_daily_clean.csv              # Processed daily panel dataset (35,090 rows, 38 columns)
├── scripts/                             # Complete, self-contained R scripts
│   ├── 01_import_and_clean.R
│   ├── 02_summary_statistics.R
│   ├── 03_correlation_simple_regression.R
│   ├── 04_multiple_linear_regression.R
│   ├── 05_distribution_fitting.R
│   ├── 06_one_sample_tests.R
│   ├── 07_ttests.R
│   ├── 08_chisquare_tests.R
│   └── 09_anova_crd_rbd.R
├── outputs/                             # 21 Publication-grade PNG figures and 22 CSV tables
│   ├── fig01_seasonal_pm25_boxplot.png ... fig21_anova_rbd_interaction_profile.png
│   └── table_summary_stats_overall.csv ... table_anova_crd_rbd_summary.csv
└── reports/
    └── lab_report.md                    # Master comprehensive lab report document
```

---

## 🚀 Quickstart: Reproducing the Results

To reproduce all analyses, clone this repository and execute the scripts in order using Base R:

```bash
git clone https://github.com/shaik-hasan-AS/AirPollutionReport.git
cd AirPollutionReport

# Run experiments sequentially (uses base R and MASS)
Rscript scripts/02_summary_statistics.R
Rscript scripts/03_correlation_simple_regression.R
Rscript scripts/04_multiple_linear_regression.R
Rscript scripts/05_distribution_fitting.R
Rscript scripts/06_one_sample_tests.R
Rscript scripts/07_ttests.R
Rscript scripts/08_chisquare_tests.R
Rscript scripts/09_anova_crd_rbd.R
```

All plots will be written directly to `outputs/`.
