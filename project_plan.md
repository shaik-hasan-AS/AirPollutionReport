# Probability and Statistics Lab Final Project - Plan

Working title: **What drives air pollution in Indian cities? A statistical study using R**

(Title is flexible. Change it once the findings are in.)

---

## 1. What we are trying to do

Build a big, finding-heavy lab report that applies 6-8 of the lab experiments in R to one real, recent, India-related dataset, and ends with conclusions that are actually useful (not just "p < 0.05").

### Hard requirements (from the assignment and our own rules)

- Real-world data, not synthetic, and recently published
- Related to India
- Useful: the findings should mean something to a reader
- Niche but not too niche
- Imported/downloaded by us, analysed in R
- 6-8 experiments out of the 11 on the lab list
- A pretty big report with strong findings

### The 11 lab experiments (reference)

1. Data types, importing/exporting data
2. Summary statistics, tables, graphs
3. Correlation and simple linear regression (R squared)
4. Multiple linear regression (multiple R squared)
5. Binomial distribution fitting
6. Normal and Poisson distribution fitting
7. Hypothesis tests, one sample mean and proportion
8. Hypothesis tests, two sample means and proportions
9. t-test, independent and dependent samples
10. Chi-square, goodness of fit and contingency
11. ANOVA: CRD, RBD, Latin square

---

## 2. What we are doing (decisions so far)

### Topic
Air quality in Indian cities: pollution levels, seasons, weather, festivals, and what predicts bad air.

### Dataset (chosen)
**Air Quality Dataset: Indian Cities (2022-2025)** on Kaggle
(kaggle.com/datasets/bhautikvekariya21/air-quality-dataset-indian-cities-2022-2025)

What we know from the listing and the sample rows:

- About 842,160 hourly records, 29 cities (state capitals plus Delhi, includes Chennai), 2022 to 2025
- 63 columns
- Identifiers: City, State, Latitude, Longitude, Datetime
- Time features: Year, Month, Day, Hour, Day_of_Week, Day_Name, Week_of_Year, Is_Weekend, Quarter, Season, Time_of_Day
- Weather: temperature at 2m/80m/120m/180m, humidity, dew point, wind speed (10m/80m/120m), wind direction, gusts, precipitation, rain flags, pressure, radiation, UV, cloud cover (total/low/mid/high), sunshine seconds
- Pollution: PM2.5, PM10, PM_Ratio, CO, NO2, SO2, O3, NH3, Dust, AOD
- AQI columns: US AQI (overall and per pollutant), EU AQI, AQI_Category, PM25_Category_India
- Extra flags: Temp_Inversion, Inversion_Strength_C, Festival_Period, Crop_Burning_Season

### Known caveats (we state these openly in the report)

- **It is most likely model-derived data** (the variable set looks like the Open-Meteo / CAMS atmospheric model output), not readings from ground sensors. It is real, not synthetic, but the methods section must say plainly that these are modelled air quality values.
- Some columns have blanks (AQI columns in the first ~24 hours, height-based temperature and wind columns). This is useful for the data cleaning part of experiment 1.
- The extra flags (Festival_Period, Crop_Burning_Season, Temp_Inversion, categories) were created by the dataset uploader. We must check how they are defined before building findings on them.
- Hourly data from one city is strongly autocorrelated. Classical tests assume independence, so we aggregate before testing.
- With 800k+ rows, almost everything becomes "significant". We report effect sizes and confidence intervals, not just p-values.

### Decision: single dataset
We are using only this dataset (no CPCB cross-check for now). See the flexibility section if we want to add one later.

---

## 3. Approach

### Data handling
1. Read with `data.table::fread()` (the file is large, around 400 MB)
2. Check structure, missing values per column and per city, date ranges per city
3. Create a **daily dataset**: one row per city per day (mean PM2.5, PM10, gases, weather; sum of rain; max of AQI if needed)
4. Keep the hourly data for the few analyses that need it (day vs night, time of day)
5. Export the cleaned daily file as CSV (this also covers the "exporting data" part of experiment 1)

### Standards to compare against
- India NAAQS PM2.5: 40 ug/m3 annual, 60 ug/m3 for 24 hours
- WHO 2021 guideline PM2.5: 5 ug/m3 annual, 15 ug/m3 for 24 hours

(Double check these numbers before they go in the final report.)

### R packages (likely)
`data.table`, `tidyverse` (dplyr, ggplot2, lubridate), `MASS` (fitdistr), `car` (VIF, Levene), `corrplot`, `broom`, `knitr`/`rmarkdown` for the report

---

## 4. Experiment plan

Target: **8 experiments**. The list below is the default. Anything marked "swap" can be replaced (see section 6).

| # | Experiment | Analysis | Main variables | Possible finding |
|---|-----------|----------|----------------|------------------|
| 1 | Data types, import/export | Classify every column, import, clean blanks, export daily file | All | Data cleaning story |
| 2 | Summary stats and plots | Mean, median, SD, skewness by city/season; histograms, boxplots, city x month heatmap | PM2.5, PM10, AQI, Season | Worst cities and worst months |
| 3 | Correlation and simple regression | PM2.5 vs PM10, wind, humidity; interpret R squared | PM2.5, PM10, weather | How much one variable explains |
| 4 | Multiple regression | PM2.5 ~ weather + gases; multiple and adjusted R squared; VIF check | Weather, NO2, SO2, O3, Dust | Strongest predictor of PM2.5 |
| 6 | Normal / Poisson fitting | Normal and log-normal fit of daily PM2.5; Poisson fit of "bad air days per month" | PM2.5, categories | Raw PM2.5 is skewed, log fits better |
| 7 | One-sample tests | Is a city's mean PM2.5 above the India/WHO limits? Proportion of days over the daily limit | PM2.5 | Which cities fail the standard |
| 9 | t-tests (independent and paired) | Weekday vs weekend, rainy vs dry, festival vs normal (paired by city) | PM2.5, Is_Weekend, Festival_Period | Diwali effect, rain effect |
| 10 | Chi-square | Goodness of fit: bad days across seasons. Contingency: AQI category vs season or crop burning season | Category columns | Seasonality of bad air |
| 11 | ANOVA | CRD: PM2.5 across seasons; RBD: season as treatment, city as block; Tukey post-hoc | PM2.5, Season, City | Monsoon vs winter gap |

Optional extras (swap candidates):

- **Experiment 5 (binomial):** number of "Good" days out of 30 per month, fitted with a binomial model
- **Experiment 8 (two-sample proportions):** coastal vs inland share of poor days
- **Latin square (part of 11):** hard to fit naturally, so treat as optional

### Findings we want to go after

1. Which city has the worst air and in which month
2. Is the Diwali effect real, and does it hold in every city
3. Does rain actually clean the air
4. Do temperature inversions trap pollution
5. Is there a weekend effect
6. What best predicts PM2.5 (weather vs gases)
7. Do coastal cities breathe better than inland ones
8. Is crop burning season visible in the data
9. Do the AQI systems (US, EU, India) agree with each other

---

## 5. Next steps

### Step 1 - Get context (before any code)
- [ ] Confirm: solo or group
- [ ] Confirm: deadline
- [ ] Confirm: report format (Word, PDF, or R Markdown output) and expected length
- [ ] Download the Kaggle CSV
- [ ] In R, run and share: `str(df)`, `colSums(is.na(df))`, `table(df$City)`
- [ ] Check the Kaggle page for the data source notes and license

### Step 2 - Validate the data
- [ ] Missing values per column and per city
- [ ] Date range and row count per city (are all 29 cities complete?)
- [ ] Check how Festival_Period, Crop_Burning_Season and Temp_Inversion are defined
- [ ] Spot check a few values (for example PM2.5 range, negative or impossible values)

### Step 3 - Lock the experiments
- [ ] Pick the final 6-8 experiments from the table
- [ ] Write a one-line hypothesis for each (H0 and H1)
- [ ] Decide which cities to feature (for example Chennai, Delhi, Mumbai, Kolkata plus a few more)

### Step 4 - Code, one experiment at a time
- [ ] One R script per experiment (`01_import_clean.R`, `02_summary.R`, and so on)
- [ ] Each script saves its plots and tables to an `outputs/` folder
- [ ] Check assumptions every time (normality, equal variance, independence)

### Step 5 - Write the report
- [ ] Draft: introduction, data description, methods, one section per experiment, discussion, limitations, conclusion
- [ ] Add plots and tables from the outputs folder
- [ ] Write in plain, natural language. No filler.
- [ ] Final read-through for consistency of numbers

### Suggested folder structure

```
project/
  data/          raw and cleaned CSVs
  scripts/       01_import_clean.R ... 09_anova.R
  outputs/       figures and tables
  report/        final report
  project_plan.md
```

---

## 6. Flexibility

This plan is a starting point, not a contract. Things we can change without starting over:

**If the dataset has problems**
- Some cities have big gaps: drop those cities or restrict to the cities with full coverage
- Flags (Festival_Period, Crop_Burning_Season) are unreliable: define our own using known festival dates
- Model-based data is not acceptable to the professor: add official CPCB station data for Chennai and Delhi (OpenCity files) and use it as the main source or as validation. A bonus finding would be a regression of CPCB PM2.5 against this dataset's PM2.5.
- The dataset turns out unusable: fall back to the Kaggle 2010-2023 CPCB city data, or switch topic to NFHS district data or VAHAN EV registrations

**If time is short**
- Cut to 6 experiments: keep 2, 3, 4, 7, 9, 11
- Use fewer cities (for example 6 to 8)
- Use daily data only and skip the hourly analyses

**If we have extra time**
- Add experiments 5 and 8
- Add the CPCB cross-check
- Add a map of mean PM2.5 by city
- Add a short "policy takeaways" section

**If a test does not give a clear result**
- Report it anyway. A null result (for example no weekend effect) is still a finding, as long as the test is done properly and we explain it.

**If the focus should change**
- We can narrow to a Tamil Nadu / Chennai angle, or widen to all 29 cities. The experiments stay the same, only the sample changes.

---

## 7. Open questions

- Solo or group?
- Deadline?
- Report format and length?
- How many experiments do we commit to (6, 7 or 8)?
- Is modelled air quality data acceptable for the "real data" rule, or do we add CPCB data?
- Which cities do we feature?

---

## 8. Definition of done

- 6-8 experiments completed in R with working, re-runnable scripts
- Every test has stated hypotheses, checked assumptions, and a clear interpretation
- At least 5 findings that read as real insights
- Limitations (model-based data, aggregation choices) written out honestly
- A clean, big report with figures and tables, ready to submit
