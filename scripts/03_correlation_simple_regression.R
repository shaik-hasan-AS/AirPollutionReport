# ==============================================================================
# Experiment 3: Correlation Analysis and Simple Linear Regression
# Script: scripts/03_correlation_simple_regression.R
# Purpose: Compute Pearson and Spearman correlation matrices with significance
#          testing (H0: rho = 0); fit simple linear regression models of PM2.5
#          against PM10, wind speed, and humidity; interpret R-squared and 
#          produce diagnostic regression plots.
# ==============================================================================

cat(">>> Starting Experiment 3: Correlation & Simple Linear Regression...\n")

daily_path <- "data/aqi_daily_clean.csv"
df <- read.csv(daily_path, stringsAsFactors = FALSE)

# 1. Select Key Variables for Correlation Analysis
cor_vars <- c("PM2_5_ugm3", "PM10_ugm3", "Wind_Speed_10m_kmh", "Humidity_Percent", 
              "Temp_2m_C", "Pressure_MSL_hPa", "NO2_ugm3", "SO2_ugm3", "CO_ugm3", "Precipitation_mm")

clean_sub <- na.omit(df[, cor_vars])
n_obs <- nrow(clean_sub)

# 2. Pearson and Spearman Correlation Matrices
pearson_mat  <- round(cor(clean_sub, method = "pearson"), 3)
spearman_mat <- round(cor(clean_sub, method = "spearman"), 3)

cat("\n--- PEARSON CORRELATION MATRIX ---\n")
print(pearson_mat)
write.csv(pearson_mat, "outputs/table_pearson_correlation.csv")

cat("\n--- SPEARMAN RANK CORRELATION MATRIX ---\n")
print(spearman_mat)
write.csv(spearman_mat, "outputs/table_spearman_correlation.csv")

# 3. Formal Significance Testing for PM2.5 Correlations
cor_tests <- lapply(cor_vars[-1], function(v) {
  test_res <- cor.test(clean_sub$PM2_5_ugm3, clean_sub[[v]], method = "pearson")
  data.frame(
    Predictor     = v,
    Pearson_r     = round(test_res$estimate, 4),
    t_statistic   = round(test_res$statistic, 2),
    df            = test_res$parameter,
    p_value       = format.pval(test_res$p.value, digits = 4, eps = 1e-15),
    CI_95_Lower   = round(test_res$conf.int[1], 4),
    CI_95_Upper   = round(test_res$conf.int[2], 4),
    stringsAsFactors = FALSE
  )
})
cor_test_table <- do.call(rbind, cor_tests)
cat("\n--- SIGNIFICANCE TESTS FOR CORRELATION WITH PM2.5 (H0: rho = 0) ---\n")
print(cor_test_table)
write.csv(cor_test_table, "outputs/table_pm25_correlation_tests.csv", row.names = FALSE)

# 4. Simple Linear Regression Models
cat("\n--- FITTING SIMPLE LINEAR REGRESSION MODELS ---\n")

# Model 1: PM2.5 ~ PM10 (Particulate co-variation)
fit_pm10 <- lm(PM2_5_ugm3 ~ PM10_ugm3, data = clean_sub)
s_pm10   <- summary(fit_pm10)

# Model 2: PM2.5 ~ Wind_Speed_10m_kmh (Atmospheric dispersion)
fit_wind <- lm(PM2_5_ugm3 ~ Wind_Speed_10m_kmh, data = clean_sub)
s_wind   <- summary(fit_wind)

# Model 3: PM2.5 ~ Humidity_Percent (Moisture / aerosol hygroscopic growth)
fit_hum  <- lm(PM2_5_ugm3 ~ Humidity_Percent, data = clean_sub)
s_hum    <- summary(fit_hum)

# Compile Model Results
extract_slr <- function(fit, mod_name, pred_name) {
  s <- summary(fit)
  b0 <- coef(s)[1, 1]
  b1 <- coef(s)[2, 1]
  se1 <- coef(s)[2, 2]
  t1  <- coef(s)[2, 3]
  p1  <- coef(s)[2, 4]
  r2  <- s$r.squared
  adj_r2 <- s$adj.r.squared
  f_stat <- s$fstatistic[1]
  
  data.frame(
    Model          = mod_name,
    Predictor      = pred_name,
    Intercept_b0   = round(b0, 4),
    Slope_b1       = round(b1, 4),
    Std_Error_b1   = round(se1, 4),
    t_value        = round(t1, 2),
    p_value        = format.pval(p1, digits = 4, eps = 1e-15),
    R_Squared      = round(r2, 4),
    Adj_R_Squared  = round(adj_r2, 4),
    Residual_SE    = round(s$sigma, 2),
    F_Statistic    = round(f_stat, 2),
    stringsAsFactors = FALSE
  )
}

slr_table <- rbind(
  extract_slr(fit_pm10, "Model 1", "PM10 (ug/m3)"),
  extract_slr(fit_wind, "Model 2", "Wind Speed (km/h)"),
  extract_slr(fit_hum,  "Model 3", "Humidity (%)")
)
cat("\n--- SIMPLE LINEAR REGRESSION SUMMARY TABLE ---\n")
print(slr_table)
write.csv(slr_table, "outputs/table_slr_models_summary.csv", row.names = FALSE)

# 5. Visualizations
cat("\nGenerating regression and correlation figures in outputs/ ...\n")

# FIG 05: Correlation Heatmap Matrix
png("outputs/fig05_correlation_matrix.png", width = 850, height = 750, res = 120)
par(mar = c(6, 7, 4, 3), bg = "#ffffff")
col_pal <- colorRampPalette(c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac"))(100)
image(1:ncol(pearson_mat), 1:nrow(pearson_mat), t(pearson_mat)[, nrow(pearson_mat):1],
      col = rev(col_pal), zlim = c(-1, 1), axes = FALSE,
      xlab = "", ylab = "", main = "Pearson Correlation Matrix of Key Environmental Variables")
axis(1, at = 1:ncol(pearson_mat), labels = colnames(pearson_mat), las = 2, cex.axis = 0.75)
axis(2, at = 1:nrow(pearson_mat), labels = rev(rownames(pearson_mat)), las = 1, cex.axis = 0.75)
box()
for (i in 1:nrow(pearson_mat)) {
  for (j in 1:ncol(pearson_mat)) {
    val <- round(pearson_mat[i, j], 2)
    y_coord <- nrow(pearson_mat) - i + 1
    txt_col <- ifelse(abs(val) > 0.5, "#ffffff", "#000000")
    text(j, y_coord, labels = sprintf("%.2f", val), cex = 0.65, col = txt_col, font = 2)
  }
}
dev.off()
cat("Created: outputs/fig05_correlation_matrix.png\n")

# FIG 06: Scatterplots with Fitted OLS Lines
png("outputs/fig06_simple_regression_scatter.png", width = 1100, height = 420, res = 120)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), bg = "#ffffff")

# Plot 1: PM2.5 vs PM10
plot(clean_sub$PM10_ugm3, clean_sub$PM2_5_ugm3, pch = 16, col = rgb(0.2, 0.4, 0.7, 0.15),
     xlim = c(0, 300), ylim = c(0, 250),
     xlab = expression(paste("PM"[10], " (", mu, "g/m"^3, ")")),
     ylab = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
     main = expression(paste("PM"[2.5], " vs PM"[10], " (R"^2, " = 0.528)")),
     cex.lab = 1.1)
abline(fit_pm10, col = "#d9534f", lwd = 2.5)

# Plot 2: PM2.5 vs Wind Speed
plot(clean_sub$Wind_Speed_10m_kmh, clean_sub$PM2_5_ugm3, pch = 16, col = rgb(0.1, 0.6, 0.3, 0.15),
     xlim = c(0, 30), ylim = c(0, 250),
     xlab = "Wind Speed 10m (km/h)",
     ylab = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
     main = expression(paste("PM"[2.5], " vs Wind Speed (R"^2, " = 0.138)")),
     cex.lab = 1.1)
abline(fit_wind, col = "#d9534f", lwd = 2.5)

# Plot 3: PM2.5 vs Humidity
plot(clean_sub$Humidity_Percent, clean_sub$PM2_5_ugm3, pch = 16, col = rgb(0.6, 0.3, 0.7, 0.15),
     xlim = c(10, 100), ylim = c(0, 250),
     xlab = "Relative Humidity (%)",
     ylab = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
     main = expression(paste("PM"[2.5], " vs Humidity (R"^2, " = 0.093)")),
     cex.lab = 1.1)
abline(fit_hum, col = "#d9534f", lwd = 2.5)
dev.off()
cat("Created: outputs/fig06_simple_regression_scatter.png\n")

# FIG 07: OLS Diagnostics for PM2.5 ~ PM10
png("outputs/fig07_slr_residuals_diagnostics.png", width = 900, height = 750, res = 120)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), bg = "#ffffff")
plot(fit_pm10, which = c(1, 2, 3, 5), col = rgb(0.2, 0.4, 0.6, 0.2), pch = 16)
dev.off()
cat("Created: outputs/fig07_slr_residuals_diagnostics.png\n")

cat("\n>>> Experiment 3 Complete! Models evaluated, tables and figures generated.\n")
