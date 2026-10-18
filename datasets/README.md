# Datasets Information

Due to GitHub's 100 MB file size limit, the raw hourly dataset files are excluded from git tracking via `.gitignore`:
- `INDIA_AQI_COMPLETE_20251126.csv` (~282.6 MB, 842,160 rows, 71 columns)
- `aqi_india_38cols_knn_final.csv` (~177.7 MB, 842,160 rows, 31 columns)

### Where to Download the Raw Data
The raw dataset can be downloaded directly from Kaggle:
- **Dataset**: [Air Quality Dataset: Indian Cities (2022-2025)](https://www.kaggle.com/datasets/bhautikvekariya21/air-quality-dataset-indian-cities-2022-2025)
- **Primary file to place in `datasets/`**: `INDIA_AQI_COMPLETE_20251126.csv`

### Cleaned Analytical Dataset
The processed daily dataset is version-controlled and available directly at:
- `data/aqi_daily_clean.csv` (12 MB, 35,090 rows, 38 columns)

Generated via:
```bash
Rscript scripts/01_import_and_clean.R
```
All statistical scripts (`02_...` through `09_...`) run directly on `data/aqi_daily_clean.csv`.
