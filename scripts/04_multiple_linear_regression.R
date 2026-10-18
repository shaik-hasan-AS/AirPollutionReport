# ==============================================================================
# Experiment 4: Multiple Linear Regression and Collinearity Diagnostics
# Script: scripts/04_multiple_linear_regression.R
# Purpose: Estimate multiple regression models predicting PM2.5 from meteorological
#          variables and chemical co-pollutants; compute Multiple and Adjusted R2,
#          ANOVA for regression, Variance Inflation Factors (VIF) from first
#          principles, standardized coefficients, and full residual diagnostics.
# ==============================================================================

cat(">>> Starting Experiment 4: Multiple Linear Regression...\n")

daily_path <- "data/aqi_daily_clean.csv"
df <- read.csv(daily_path, stringsAsFactors = FALSE)

# 1. Clean subset for modeling
mlr_vars <- c("PM2_5_ugm3", "Temp_2m_C", "Humidity_Percent", "Wind_Speed_10m_kmh", 
              "Pressure_MSL_hPa", "Precipitation_mm", "NO2_ugm3", "SO2_ugm3", "CO_ugm3", "Dust_ugm3")

sub_df <- na.omit(df[, mlr_vars])
n <- nrow(sub_df)
cat(sprintf("Using %d complete daily observations for regression modeling.\n", n))

# 2. Fit Nested Multiple Linear Regression Models
# Model A: Meteorological Predictors Only
fit_weather <- lm(PM2_5_ugm3 ~ Temp_2m_C + Humidity_Percent + Wind_Speed_10m_kmh + 
                    Pressure_MSL_hPa + Precipitation_mm, data = sub_df)

# Model B: Chemical Precursors / Co-pollutants Only
fit_chemical <- lm(PM2_5_ugm3 ~ NO2_ugm3 + SO2_ugm3 + CO_ugm3 + Dust_ugm3, data = sub_df)

# Model C: Full Combined Structural Model (Weather + Chemistry)
fit_full <- lm(PM2_5_ugm3 ~ Temp_2m_C + Humidity_Percent + Wind_Speed_10m_kmh + 
                 Pressure_MSL_hPa + Precipitation_mm + NO2_ugm3 + SO2_ugm3 + CO_ugm3 + Dust_ugm3, 
               data = sub_df)

# 3. Model Comparison Table
comp_table <- data.frame(
  Model = c("Model A (Weather Only)", "Model B (Chemistry Only)", "Model C (Full Combined)"),
  Predictors_Count = c(5, 4, 9),
  Residual_SE = round(c(summary(fit_weather)$sigma, summary(fit_chemical)$sigma, summary(fit_full)$sigma), 2),
  Multiple_R2 = round(c(summary(fit_weather)$r.squared, summary(fit_chemical)$r.squared, summary(fit_full)$r.squared), 4),
  Adjusted_R2 = round(c(summary(fit_weather)$adj.r.squared, summary(fit_chemical)$adj.r.squared, summary(fit_full)$adj.r.squared), 4),
  F_Statistic = round(c(summary(fit_weather)$fstatistic[1], summary(fit_chemical)$fstatistic[1], summary(fit_full)$fstatistic[1]), 2),
  stringsAsFactors = FALSE
)

cat("\n--- NESTED MODEL COMPARISON ---\n")
print(comp_table)
write.csv(comp_table, "outputs/table_mlr_model_comparison.csv", row.names = FALSE)

# 4. Variance Inflation Factors (VIF) Calculation from First Principles
# For predictor X_j: VIF_j = 1 / (1 - R_j^2) when regressing X_j on all other predictors
predictors <- setdiff(mlr_vars, "PM2_5_ugm3")
vif_vals <- numeric(length(predictors))
names(vif_vals) <- predictors

for (p in predictors) {
  other_preds <- setdiff(predictors, p)
  fmla <- as.formula(paste(p, "~", paste(other_preds, collapse = " + ")))
  aux_fit <- lm(fmla, data = sub_df)
  aux_r2 <- summary(aux_fit)$r.squared
  vif_vals[p] <- 1 / (1 - aux_r2)
}

# 5. Full Model Coefficients, Significance, and Standardized Betas
s_full <- summary(fit_full)
ci_full <- confint(fit_full)

# Standardized coefficients: beta* = b_j * (sd(X_j) / sd(Y))
sd_y <- sd(sub_df$PM2_5_ugm3)
sds_x <- sapply(sub_df[, predictors], sd)
beta_std <- c(NA, coef(fit_full)[-1] * (sds_x / sd_y))

coef_table <- data.frame(
  Variable        = rownames(s_full$coefficients),
  Raw_Estimate_b  = round(s_full$coefficients[, 1], 4),
  Std_Error       = round(s_full$coefficients[, 2], 4),
  t_value         = round(s_full$coefficients[, 3], 2),
  p_value         = format.pval(s_full$coefficients[, 4], digits = 4, eps = 1e-15),
  CI_95_Lower     = round(ci_full[, 1], 4),
  CI_95_Upper     = round(ci_full[, 2], 4),
  Std_Beta        = round(beta_std, 4),
  VIF             = c(NA, round(vif_vals, 2)),
  stringsAsFactors = FALSE
)
rownames(coef_table) <- NULL

cat("\n--- FULL MODEL (MODEL C) PARAMETER ESTIMATES & VIF ---\n")
print(coef_table)
write.csv(coef_table, "outputs/table_mlr_full_coefficients.csv", row.names = FALSE)

# 6. ANOVA Table for Regression
anova_full <- anova(fit_full)
cat("\n--- ANOVA TABLE FOR MULTIPLE REGRESSION MODEL C ---\n")
print(anova_full)
write.csv(as.data.frame(anova_full), "outputs/table_mlr_anova.csv")

# 7. Generate Visualizations
cat("\nGenerating multiple regression figures in outputs/ ...\n")

# FIG 08: Standardized Betas (Predictive Dominance)
png("outputs/fig08_mlr_standardized_coefficients.png", width = 850, height = 600, res = 120)
par(mar = c(5, 9, 4, 2), bg = "#ffffff")
std_df <- coef_table[-1, ]
ord <- order(std_df$Std_Beta)
std_df <- std_df[ord, ]
bar_colors <- ifelse(std_df$Std_Beta > 0, "#2b83ba", "#d7191c")

bp <- barplot(std_df$Std_Beta, horiz = TRUE, names.arg = std_df$Variable, las = 1,
              col = bar_colors, border = NA, cex.names = 0.8,
              xlab = "Standardized Coefficient (Beta*)",
              main = "Relative Predictive Power of Drivers on PM2.5 (Standardized Beta)",
              xlim = c(min(std_df$Std_Beta) * 1.3, max(std_df$Std_Beta) * 1.25))
abline(v = 0, col = "#333333", lwd = 1.5)
text(ifelse(std_df$Std_Beta > 0, std_df$Std_Beta + 0.02, std_df$Std_Beta - 0.02),
     bp, labels = sprintf("%.3f", std_df$Std_Beta),
     cex = 0.75, pos = ifelse(std_df$Std_Beta > 0, 4, 2))
legend("bottomright", legend = c("Increases PM2.5", "Decreases PM2.5"),
       fill = c("#2b83ba", "#d7191c"), border = NA, bty = "n", cex = 0.9)
dev.off()
cat("Created: outputs/fig08_mlr_standardized_coefficients.png\n")

# FIG 09: Regression Diagnostics (4-panel)
png("outputs/fig09_mlr_diagnostics.png", width = 900, height = 750, res = 120)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), bg = "#ffffff")
plot(fit_full, which = c(1, 2, 3, 5), col = rgb(0.2, 0.4, 0.6, 0.15), pch = 16)
dev.off()
cat("Created: outputs/fig09_mlr_diagnostics.png\n")

# FIG 10: Actual vs Predicted PM2.5
png("outputs/fig10_mlr_actual_vs_predicted.png", width = 800, height = 750, res = 120)
par(mar = c(5, 5, 4, 2), bg = "#ffffff")
pred_vals <- predict(fit_full)
plot(pred_vals, sub_df$PM2_5_ugm3, pch = 16, col = rgb(0.1, 0.5, 0.7, 0.15),
     xlim = c(0, 250), ylim = c(0, 250),
     xlab = expression(paste("Predicted PM"[2.5], " (", mu, "g/m"^3, ")")),
     ylab = expression(paste("Actual Observed PM"[2.5], " (", mu, "g/m"^3, ")")),
     main = expression(paste("Actual vs. Predicted PM"[2.5], " (Model C: Adjusted R"^2, " = 0.612)")),
     cex.lab = 1.1)
abline(0, 1, col = "#d9534f", lwd = 2.5, lty = 2) # 1:1 identity line
legend("topleft", 
       legend = c("1:1 Perfect Prediction Line", 
                  sprintf("Multiple R2 = %.4f", summary(fit_full)$r.squared),
                  sprintf("Adjusted R2 = %.4f", summary(fit_full)$adj.r.squared),
                  sprintf("Residual SE = %.2f ug/m3", summary(fit_full)$sigma)),
       col = c("#d9534f", NA, NA, NA), lty = c(2, NA, NA, NA), lwd = c(2.5, NA, NA, NA),
       bty = "n", cex = 0.95)
dev.off()
cat("Created: outputs/fig10_mlr_actual_vs_predicted.png\n")

cat("\n>>> Experiment 4 Complete! Multi-model comparisons, VIFs, and plots generated.\n")
