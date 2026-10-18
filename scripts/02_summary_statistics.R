# ==============================================================================
# Experiment 2: Summary Statistics, Dispersion Metrics, and Visualizations
# Script: scripts/02_summary_statistics.R
# Purpose: Calculate parametric and non-parametric summary statistics (mean,
#          median, trimmed mean, variance, SD, IQR, CV, skewness, kurtosis) 
#          overall, by city, and by season; generate publication-grade figures.
# ==============================================================================

cat(">>> Starting Experiment 2: Summary Statistics & Visualizations...\n")

# 1. Load clean daily dataset
daily_path <- "data/aqi_daily_clean.csv"
if (!file.exists(daily_path)) {
  stop("Daily cleaned dataset not found. Please run 01_import_and_clean.R first.")
}

df <- read.csv(daily_path, stringsAsFactors = FALSE)
cat(sprintf("Loaded daily dataset: %d observations across %d variables.\n", nrow(df), ncol(df)))

# 2. Statistical Helper Functions for Distribution Moments
calc_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m3 <- sum((x - mean(x))^3) / n
  s3 <- (sd(x))^3
  return(m3 / s3)
}

calc_kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m4 <- sum((x - mean(x))^4) / n
  s4 <- (sd(x))^4
  return(m4 / s4) # Pearson kurtosis (Normal = 3)
}

compute_moments <- function(x, var_name = "Variable") {
  x <- x[!is.na(x)]
  n <- length(x)
  mean_val   <- mean(x)
  sd_val     <- sd(x)
  median_val <- median(x)
  q25        <- quantile(x, 0.25)
  q75        <- quantile(x, 0.75)
  iqr_val    <- q75 - q25
  trimmed_05 <- mean(x, trim = 0.05)
  cv_val     <- (sd_val / mean_val) * 100
  skew_val   <- calc_skewness(x)
  kurt_val   <- calc_kurtosis(x)
  
  data.frame(
    Variable       = var_name,
    N              = n,
    Mean           = round(mean_val, 2),
    Trimmed_Mean_5 = round(trimmed_05, 2),
    Median         = round(median_val, 2),
    Variance       = round(var(x), 2),
    Std_Dev        = round(sd_val, 2),
    IQR            = round(iqr_val, 2),
    CV_Percent     = round(cv_val, 2),
    Skewness       = round(skew_val, 2),
    Kurtosis       = round(kurt_val, 2),
    Min            = round(min(x), 2),
    Max            = round(max(x), 2),
    stringsAsFactors = FALSE
  )
}

# 3. Overall Descriptive Statistics for Key Environmental Indicators
vars_to_summarize <- c("PM2_5_ugm3", "PM10_ugm3", "CO_ugm3", "NO2_ugm3", "SO2_ugm3", "O3_ugm3", 
                       "Temp_2m_C", "Humidity_Percent", "Wind_Speed_10m_kmh", "Precipitation_mm")

overall_list <- lapply(vars_to_summarize, function(v) compute_moments(df[[v]], v))
overall_table <- do.call(rbind, overall_list)

cat("\n--- OVERALL DESCRIPTIVE SUMMARY STATISTICS ---\n")
print(overall_table)
write.csv(overall_table, "outputs/table_summary_stats_overall.csv", row.names = FALSE)

# 4. City-Wise PM2.5 Summary Table
cities <- sort(unique(df$City))
city_list <- lapply(cities, function(cty) {
  sub_x <- df$PM2_5_ugm3[df$City == cty]
  res <- compute_moments(sub_x, cty)
  res$Pct_Exceed_NAAQS <- round(mean(sub_x > 60) * 100, 2)
  res$Pct_Exceed_WHO   <- round(mean(sub_x > 15) * 100, 2)
  res
})
city_table <- do.call(rbind, city_list)
names(city_table)[1] <- "City"
city_table <- city_table[order(-city_table$Mean), ]

cat("\n--- TOP 10 POLLUTED CITIES BY MEAN PM2.5 ---\n")
print(head(city_table[, c("City", "Mean", "Median", "Std_Dev", "Skewness", "Pct_Exceed_NAAQS")], 10))
write.csv(city_table, "outputs/table_summary_stats_by_city.csv", row.names = FALSE)

# 5. Seasonal Summary Table
seasons <- c("Winter", "Summer", "Monsoon", "Post_Monsoon")
season_list <- lapply(seasons, function(sea) {
  sub_x <- df$PM2_5_ugm3[df$Season == sea]
  res <- compute_moments(sub_x, sea)
  res$Pct_Exceed_NAAQS <- round(mean(sub_x > 60) * 100, 2)
  res
})
season_table <- do.call(rbind, season_list)
names(season_table)[1] <- "Season"
rownames(season_table) <- NULL

cat("\n--- SEASONAL SUMMARY STATISTICS FOR PM2.5 ---\n")
print(season_table[, c("Season", "N", "Mean", "Median", "Std_Dev", "IQR", "Pct_Exceed_NAAQS")])
write.csv(season_table, "outputs/table_summary_stats_by_season.csv", row.names = FALSE)

# 6. Generate Publication-Quality Visualizations
cat("\nGenerating publication-grade plots in outputs/ ...\n")

# FIG 01: Seasonal Boxplot
png("outputs/fig01_seasonal_pm25_boxplot.png", width = 900, height = 650, res = 120)
par(mar = c(5, 5, 4, 2), bg = "#ffffff")
season_colors <- c("#4575b4", "#fee090", "#91bfdb", "#d73027")
boxplot(PM2_5_ugm3 ~ Season, data = df,
        col = season_colors,
        border = "#2c3e50",
        outline = FALSE, # focus on distribution IQR without extreme tails
        notch = TRUE,
        ylab = expression(paste("Daily PM"[2.5], " (", mu, "g/m"^3, ")")),
        xlab = "Season",
        main = expression(paste("Distribution of Daily PM"[2.5], " Across Seasons (2022-2025)")),
        font.main = 2, cex.main = 1.2, cex.lab = 1.1)
abline(h = 60, col = "#d9534f", lwd = 2, lty = 2) # Indian NAAQS
abline(h = 15, col = "#5cb85c", lwd = 2, lty = 3) # WHO Guideline
legend("topleft", 
       legend = c("India NAAQS 24-hr Limit (60 ug/m3)", "WHO 24-hr Guideline (15 ug/m3)"),
       col = c("#d9534f", "#5cb85c"), lty = c(2, 3), lwd = 2, bty = "n", cex = 0.9)
dev.off()
cat("Created: outputs/fig01_seasonal_pm25_boxplot.png\n")

# FIG 02: City-wise Mean PM2.5 Ranking
png("outputs/fig02_city_mean_pm25_barchart.png", width = 1000, height = 800, res = 120)
par(mar = c(5, 8, 4, 2), bg = "#ffffff")
ord_cities <- city_table[order(city_table$Mean), ]
bar_cols <- ifelse(ord_cities$Mean > 60, "#d9534f", "#337ab7")
bp <- barplot(ord_cities$Mean, horiz = TRUE, names.arg = ord_cities$City, las = 1,
              col = bar_cols, border = NA, cex.names = 0.8,
              xlab = expression(paste("Mean Daily PM"[2.5], " (", mu, "g/m"^3, ")")),
              main = expression(paste("Ranking of 29 Indian Cities by Mean PM"[2.5], " (2022-2025)")),
              xlim = c(0, max(ord_cities$Mean) * 1.15))
abline(v = 60, col = "#c9302c", lwd = 2, lty = 2)
abline(v = 15, col = "#4cae4c", lwd = 2, lty = 3)
text(ord_cities$Mean + 2, bp, labels = sprintf("%.1f", ord_cities$Mean), cex = 0.75, pos = 4)
legend("bottomright", 
       legend = c("Exceeds NAAQS Annual/Daily Limit (>60)", "Complies with Daily NAAQS (<=60)"),
       fill = c("#d9534f", "#337ab7"), border = NA, bty = "n", cex = 0.9)
dev.off()
cat("Created: outputs/fig02_city_mean_pm25_barchart.png\n")

# FIG 03: Histogram & Distribution Skewness
png("outputs/fig03_pm25_histogram_density.png", width = 900, height = 650, res = 120)
par(mar = c(5, 5, 4, 2), bg = "#ffffff")
h <- hist(df$PM2_5_ugm3, breaks = 60, freq = FALSE,
          col = "#ecf0f1", border = "#bdc3c7",
          xlim = c(0, 250),
          xlab = expression(paste("Daily PM"[2.5], " Concentration (", mu, "g/m"^3, ")")),
          main = expression(paste("Empirical Distribution of Daily PM"[2.5], " (Positive Skewness)")),
          cex.lab = 1.1, font.main = 2)
# Overlaid kernel density
lines(density(df$PM2_5_ugm3, na.rm = TRUE), col = "#e74c3c", lwd = 2.5)
# Theoretical normal distribution with same mean & sd
x_seq <- seq(0, 250, length.out = 500)
lines(x_seq, dnorm(x_seq, mean = mean(df$PM2_5_ugm3), sd = sd(df$PM2_5_ugm3)), 
      col = "#2980b9", lwd = 2, lty = 2)
abline(v = mean(df$PM2_5_ugm3), col = "#2c3e50", lwd = 2, lty = 1)
abline(v = median(df$PM2_5_ugm3), col = "#8e44ad", lwd = 2, lty = 4)
legend("topright",
       legend = c(
         sprintf("Empirical Kernel Density (Skew = %.2f)", calc_skewness(df$PM2_5_ugm3)),
         "Theoretical Normal Curve",
         sprintf("Mean = %.1f ug/m3", mean(df$PM2_5_ugm3)),
         sprintf("Median = %.1f ug/m3", median(df$PM2_5_ugm3))
       ),
       col = c("#e74c3c", "#2980b9", "#2c3e50", "#8e44ad"),
       lty = c(1, 2, 1, 4), lwd = c(2.5, 2, 2, 2), bty = "n", cex = 0.9)
dev.off()
cat("Created: outputs/fig03_pm25_histogram_density.png\n")

# FIG 04: City x Season Heatmap Matrix
png("outputs/fig04_city_season_heatmap.png", width = 800, height = 900, res = 120)
par(mar = c(6, 9, 4, 3), bg = "#ffffff")
mat <- tapply(df$PM2_5_ugm3, list(df$City, df$Season), mean)
# Order seasons logically
mat <- mat[, c("Winter", "Summer", "Monsoon", "Post_Monsoon")]
colnames(mat) <- c("Winter", "Summer", "Monsoon", "Post-Monsoon")
# Order cities by Winter pollution
mat <- mat[order(mat[, "Winter"]), ]

color_palette <- colorRampPalette(c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c"))(50)
image(1:ncol(mat), 1:nrow(mat), t(mat), 
      col = color_palette, axes = FALSE,
      xlab = "", ylab = "",
      main = expression(paste("Mean Daily PM"[2.5], " by City and Season (", mu, "g/m"^3, ")")))
axis(1, at = 1:ncol(mat), labels = colnames(mat), las = 1, font = 2)
axis(2, at = 1:nrow(mat), labels = rownames(mat), las = 1, cex.axis = 0.8)
box()

# Add text values in cells
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    val <- round(mat[i, j], 1)
    txt_col <- ifelse(val > 70, "#ffffff", "#000000")
    text(j, i, labels = val, cex = 0.65, col = txt_col, font = 2)
  }
}
dev.off()
cat("Created: outputs/fig04_city_season_heatmap.png\n")

cat("\n>>> Experiment 2 Complete! Generated 3 tables and 4 publication-quality figures.\n")
