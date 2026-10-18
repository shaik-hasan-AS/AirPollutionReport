# ==============================================================================
# Experiment 1: Data Types, Importing, Cleaning, and Aggregation
# Script: scripts/01_import_and_clean.R
# Purpose: Import raw 71-column hourly AQI dataset, classify data types, 
#          audit missing values, drop unpopulated upper-air variables,
#          aggregate to daily granularity, and export clean analytical dataset.
# ==============================================================================

cat(">>> Starting Experiment 1: Data Import and Cleaning...\n")
start_time <- Sys.time()

# 1. Define paths
raw_path   <- "datasets/INDIA_AQI_COMPLETE_20251126.csv"
daily_path <- "data/aqi_daily_clean.csv"
audit_path <- "outputs/column_audit_exp01.csv"

# Ensure output directories exist
dir.create("data", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# 2. Ingest raw hourly dataset
cat("Reading raw dataset:", raw_path, "...\n")
raw_df <- read.csv(raw_path, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d rows and %d columns in %.2f seconds.\n", 
            nrow(raw_df), ncol(raw_df), as.numeric(difftime(Sys.time(), start_time, units="secs"))))

# 3. Classify data types and audit missingness
cat("Auditing column types and missing value counts...\n")
col_names  <- names(raw_df)
col_types  <- sapply(raw_df, class)
na_counts  <- sapply(raw_df, function(x) sum(is.na(x) | x == ""))
na_percent <- round((na_counts / nrow(raw_df)) * 100, 2)

audit_df <- data.frame(
  Column_Name = col_names,
  Data_Type   = col_types,
  Missing_Count = na_counts,
  Missing_Pct = na_percent,
  stringsAsFactors = FALSE
)

# Identify completely unpopulated columns (100% missing)
unpopulated_cols <- audit_df$Column_Name[audit_df$Missing_Pct == 100]
cat("Columns with 100% missing values identified:\n")
print(unpopulated_cols)

# Write column audit report
write.csv(audit_df, audit_path, row.names = FALSE)
cat("Audit saved to:", audit_path, "\n")

# 4. Clean and drop 100% missing columns
cat("Dropping 100% unpopulated upper-air and unrecorded columns...\n")
clean_df <- raw_df[, !(names(raw_df) %in% unpopulated_cols)]
cat(sprintf("Retained %d active columns.\n", ncol(clean_df)))

# 5. Extract Date string for daily grouping (Format: YYYY-MM-DD)
clean_df$Date <- substr(clean_df$Datetime, 1, 10)

# 6. Aggregation to Daily Granularity
# Aggregation strategy:
# - Pollutants and continuous weather: Arithmetic daily mean
# - Cumulative weather (Precipitation, Rain): Daily sum
# - Extreme weather / AQI: Daily maximum
# - Metadata & Flags: First/representative daily value
cat("Aggregating 842,160 hourly records into daily city-level summaries...\n")

# Numeric columns to average
mean_cols <- c(
  "PM2_5_ugm3", "PM10_ugm3", "PM_Ratio", "CO_ugm3", "NO2_ugm3", 
  "SO2_ugm3", "O3_ugm3", "Dust_ugm3", "AOD",
  "Temp_2m_C", "Humidity_Percent", "Dew_Point_C", "Wind_Speed_10m_kmh", 
  "Wind_Dir_10m", "Pressure_MSL_hPa", "Cloud_Cover_Percent", "Solar_Radiation_Wm2"
)

# Aggregate means by City, State, Date
means_agg <- aggregate(clean_df[, mean_cols], 
                       by = list(City = clean_df$City, State = clean_df$State, Date = clean_df$Date), 
                       FUN = function(x) mean(as.numeric(x), na.rm = TRUE))

# Numeric columns to sum (rain)
sum_cols <- c("Precipitation_mm", "Rain_mm")
sums_agg <- aggregate(clean_df[, sum_cols], 
                      by = list(City = clean_df$City, State = clean_df$State, Date = clean_df$Date), 
                      FUN = function(x) sum(as.numeric(x), na.rm = TRUE))

# Numeric columns to take daily maximum
max_cols <- c("Wind_Gusts_kmh", "US_AQI")
max_agg  <- aggregate(clean_df[, max_cols], 
                      by = list(City = clean_df$City, State = clean_df$State, Date = clean_df$Date), 
                      FUN = function(x) {
                        v <- suppressWarnings(as.numeric(x[!is.na(x) & x != ""]))
                        if (length(v) == 0) return(NA) else max(v)
                      })

# Categorical / contextual flags (take modal or first observation of day)
cat_cols <- c("Year", "Month", "Day", "Day_Name", "Is_Weekend", "Quarter", "Season", 
              "Festival_Period", "Crop_Burning_Season", "Latitude", "Longitude")
first_agg <- aggregate(clean_df[, cat_cols],
                       by = list(City = clean_df$City, State = clean_df$State, Date = clean_df$Date),
                       FUN = function(x) x[1])

# Merge aggregated frames
daily_clean <- merge(means_agg, sums_agg, by = c("City", "State", "Date"))
daily_clean <- merge(daily_clean, max_agg, by = c("City", "State", "Date"))
daily_clean <- merge(daily_clean, first_agg, by = c("City", "State", "Date"))

# 7. Add derived regulatory compliance and meteorological indicator columns
# India NAAQS 24-hr standard for PM2.5 = 60 ug/m3
# WHO 2021 24-hr guideline for PM2.5   = 15 ug/m3
daily_clean$Exceed_NAAQS <- ifelse(daily_clean$PM2_5_ugm3 > 60, 1, 0)
daily_clean$Exceed_WHO   <- ifelse(daily_clean$PM2_5_ugm3 > 15, 1, 0)
daily_clean$Has_Rain     <- ifelse(daily_clean$Precipitation_mm > 0.1, 1, 0)

# Sort chronologically by City and Date
daily_clean <- daily_clean[order(daily_clean$City, daily_clean$Date), ]

# 8. Export clean daily dataset
cat("Exporting cleaned daily dataset to:", daily_path, "...\n")
write.csv(daily_clean, daily_path, row.names = FALSE)

# 9. Verification Summary
elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
cat("\n================================================================\n")
cat("EXPERIMENT 1 COMPLETE: DATA INGESTION, AUDIT & EXPORT\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("Raw Hourly Observations : %d across 29 cities\n", nrow(raw_df)))
cat(sprintf("Raw Columns             : %d\n", ncol(raw_df)))
cat(sprintf("Dropped 100%% NA Columns : %d (%s)\n", length(unpopulated_cols), paste(unpopulated_cols, collapse=", ")))
cat(sprintf("Active Columns Retained : %d\n", ncol(clean_df)))
cat(sprintf("Aggregated Daily Records: %d (Expected: 29 cities * 1,210 days = 35,090)\n", nrow(daily_clean)))
cat(sprintf("Daily Export Features   : %d columns\n", ncol(daily_clean)))
cat(sprintf("Total Elapsed Runtime   : %.2f seconds\n", elapsed))
cat("================================================================\n")
