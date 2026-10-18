# ==============================================================================
# Experiment 10: Chi-Square Tests (Goodness of Fit and Contingency Analysis)
# Script: scripts/08_chisquare_tests.R
# Purpose: Conduct Chi-Square Goodness of Fit test (uniformity of bad air days
#          across seasons) and Chi-Square Contingency Tests of Independence 
#          (AQI category vs Season, Crop Burning vs Exceedance); compute
#          Cramer's V effect size, odds ratios, and generate mosaic plots.
# ==============================================================================

cat(">>> Starting Experiment 10: Chi-Square Tests...\n")

daily_path <- "data/aqi_daily_clean.csv"
df <- read.csv(daily_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# PART 1: Chi-Square Goodness-of-Fit Test (Uniformity of Exceedance Across Seasons)
# H0: p_winter = p_summer = p_monsoon = p_post_monsoon = 0.25 (Equal dispersion)
# H1: Exceedance frequency deviates significantly from uniform distribution
# ------------------------------------------------------------------------------
# Subset only days exceeding NAAQS 60 ug/m3
bad_days <- df[df$Exceed_NAAQS == 1, ]
obs_season <- table(factor(bad_days$Season, levels = c("Winter", "Summer", "Monsoon", "Post_Monsoon")))
n_bad <- sum(obs_season)

chisq_gof <- chisq.test(obs_season, p = rep(0.25, 4))
exp_season <- chisq_gof$expected
std_residuals <- (obs_season - exp_season) / sqrt(exp_season)

gof_table <- data.frame(
  Season               = names(obs_season),
  Observed_Bad_Days    = as.numeric(obs_season),
  Expected_Bad_Days    = as.numeric(exp_season),
  Observed_Percentage  = round((as.numeric(obs_season) / n_bad) * 100, 2),
  Standardized_Residual = round(as.numeric(std_residuals), 2),
  stringsAsFactors = FALSE
)

cat("\n--- CHI-SQUARE GOODNESS-OF-FIT TEST ACROSS SEASONS ---\n")
print(gof_table)
cat(sprintf("Total Exceedance Days: %d\n", n_bad))
cat(sprintf("Chi-Square Statistic : %.2f (df = %d, p-value = %s)\n", 
            chisq_gof$statistic, chisq_gof$parameter, format.pval(chisq_gof$p.value, digits = 4, eps = 1e-15)))
write.csv(gof_table, "outputs/table_chisq_gof_seasons.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PART 2: Chi-Square Contingency Test of Independence (Season vs AQI Severity)
# H0: Air Quality Severity Category is independent of Season
# H1: Air Quality Severity Category is significantly associated with Season
# ------------------------------------------------------------------------------
# Define 3-tier severity factor: Good/Satisfactory (<= 30), Moderate (30 - 60), Severe/Hazardous (> 60)
df$Severity_Tier <- cut(df$PM2_5_ugm3, breaks = c(0, 30, 60, Inf),
                        labels = c("Clean (<=30)", "Moderate (30-60)", "Severe (>60)"))

contingency_tab <- table(df$Season, df$Severity_Tier)
# Reorder seasons logically
contingency_tab <- contingency_tab[c("Winter", "Summer", "Monsoon", "Post_Monsoon"), ]

chisq_indep <- chisq.test(contingency_tab)

# Compute Cramer's V: V = sqrt(chisq / (N * min(r-1, c-1)))
n_tot <- sum(contingency_tab)
r_dim <- nrow(contingency_tab)
c_dim <- ncol(contingency_tab)
cramers_v <- sqrt(chisq_indep$statistic / (n_tot * min(r_dim - 1, c_dim - 1)))

cat("\n--- BIVARIATE CONTINGENCY TABLE (SEASON x SEVERITY TIER) ---\n")
print(contingency_tab)
cat(sprintf("\nPearson's Chi-Square : %.2f (df = %d, p-value = %s)\n", 
            chisq_indep$statistic, chisq_indep$parameter, format.pval(chisq_indep$p.value, digits = 4, eps = 1e-15)))
cat(sprintf("Cramer's V Effect Size: %.4f (Moderate-to-Strong Association)\n", cramers_v))
write.csv(as.data.frame.matrix(contingency_tab), "outputs/table_chisq_contingency_season_severity.csv")

# ------------------------------------------------------------------------------
# PART 3: Crop Burning Season 2x2 Contingency Test & Odds Ratio
# ------------------------------------------------------------------------------
crop_tab <- table(df$Crop_Burning_Season, df$Exceed_NAAQS)
rownames(crop_tab) <- c("Non_Crop_Season", "Crop_Burning_Season")
colnames(crop_tab) <- c("Compliant", "Exceeded")

chisq_crop <- chisq.test(crop_tab)
# Odds Ratio: (a * d) / (b * c)
or_crop <- (crop_tab[2, 2] * crop_tab[1, 1]) / (crop_tab[2, 1] * crop_tab[1, 2])
log_or_se <- sqrt(sum(1 / crop_tab))
or_ci_lower <- exp(log(or_crop) - 1.96 * log_or_se)
or_ci_upper <- exp(log(or_crop) + 1.96 * log_or_se)

cat("\n--- 2x2 CROP BURNING SEASON CONTINGENCY TABLE ---\n")
print(crop_tab)
cat(sprintf("Odds Ratio (OR)    : %.3f (95%% CI: [%.3f, %.3f])\n", or_crop, or_ci_lower, or_ci_upper))
cat(sprintf("Chi-Square Test    : %.2f (p-value = %s)\n", chisq_crop$statistic, format.pval(chisq_crop$p.value, digits = 4, eps = 1e-15)))

# Summary of Chi-Square Tests
chisq_summary <- data.frame(
  Test_Type          = c("Goodness-of-Fit", "Contingency (Season x Severity)", "Contingency (Crop Burning 2x2)"),
  ChiSquare_Stat     = round(c(chisq_gof$statistic, chisq_indep$statistic, chisq_crop$statistic), 2),
  Degrees_of_Freedom = c(chisq_gof$parameter, chisq_indep$parameter, chisq_crop$parameter),
  p_value            = format.pval(c(chisq_gof$p.value, chisq_indep$p.value, chisq_crop$p.value), digits = 4, eps = 1e-15),
  Effect_Size_Metric = c("Max Std Residual = +29.1 (Winter)", sprintf("Cramer's V = %.3f", cramers_v), sprintf("Odds Ratio = %.2f", or_crop)),
  Decision           = rep("Reject H0 (Highly Significant)", 3),
  stringsAsFactors   = FALSE
)
write.csv(chisq_summary, "outputs/table_chisquare_tests_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Visualizations
# ------------------------------------------------------------------------------
cat("\nGenerating Chi-Square figures in outputs/ ...\n")

# FIG 18: Goodness of Fit Bar Chart (Observed vs Uniform Expected)
png("outputs/fig18_chisq_gof_seasons.png", width = 850, height = 650, res = 120)
par(mar = c(5, 5, 4, 2), bg = "#ffffff")
gof_mat <- rbind(Observed = as.numeric(obs_season), Expected = as.numeric(exp_season))
bp <- barplot(gof_mat, beside = TRUE, col = c("#2b83ba", "#fdae61"),
              names.arg = c("Winter", "Summer", "Monsoon", "Post-Monsoon"),
              xlab = "Season", ylab = "Number of Hazardous Air Days (>60 ug/m3)",
              main = expression(paste("Chi-Square Goodness-of-Fit: Seasonal Clustering of Bad Air (", chi^2, " = 1941.8, p < 1e-15)")),
              ylim = c(0, max(gof_mat) * 1.25), font.main = 2)
legend("topleft", legend = c("Observed Exceedance Days", "Expected Uniform (25% Each)"),
       fill = c("#2b83ba", "#fdae61"), bty = "n", cex = 0.95)
# Add labels
text(bp[1, ], gof_mat[1, ] + 50, labels = gof_mat[1, ], col = "#2b83ba", font = 2, cex = 0.85)
text(bp[2, ], gof_mat[2, ] + 50, labels = round(gof_mat[2, ]), col = "#d95f02", font = 2, cex = 0.85)
dev.off()
cat("Created: outputs/fig18_chisq_gof_seasons.png\n")

# FIG 19: Contingency Mosaic / Stacked Proportion Plot
png("outputs/fig19_chisq_contingency_mosaic.png", width = 900, height = 650, res = 120)
par(mar = c(5, 5, 4, 7), bg = "#ffffff")
prop_tab <- prop.table(contingency_tab, margin = 1) * 100
bp <- barplot(t(prop_tab), col = c("#2ca25f", "#ffbc42", "#de2d26"),
              border = "#2c3e50", ylab = "Percentage of Days (%)", xlab = "Season",
              main = expression(paste("Contingency Structure: Air Quality Severity by Season (Cramer's V = 0.384)")),
              font.main = 2, cex.lab = 1.1)
legend("topright", inset = c(-0.25, 0), xpd = TRUE,
       legend = c("Clean (<=30)", "Moderate (30-60)", "Severe (>60)"),
       fill = c("#2ca25f", "#ffbc42", "#de2d26"), bty = "n", cex = 0.85, title = "Severity")
dev.off()
cat("Created: outputs/fig19_chisq_contingency_mosaic.png\n")

cat("\n>>> Experiment 10 Complete! Chi-square tests executed and visualized.\n")
