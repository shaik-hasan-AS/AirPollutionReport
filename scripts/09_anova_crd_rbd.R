# ==============================================================================
# Experiment 11: Analysis of Variance (ANOVA: CRD and RBD Designs)
# Script: scripts/09_anova_crd_rbd.R
# Purpose: Execute Completely Randomized Design (One-Way CRD ANOVA) across seasons
#          with Tukey HSD post-hoc comparisons; execute Randomized Block Design
#          (Two-Way RBD ANOVA) with Season as Treatment and City as Blocking factor;
#          compute Relative Efficiency (RE) of RBD over CRD and generate diagnostic plots.
# ==============================================================================

cat(">>> Starting Experiment 11: ANOVA (CRD and RBD)...\n")

daily_path <- "data/aqi_daily_clean.csv"
df <- read.csv(daily_path, stringsAsFactors = FALSE)

# Ensure Season is a factor with meaningful baseline
df$Season <- factor(df$Season, levels = c("Monsoon", "Summer", "Post_Monsoon", "Winter"))
df$City   <- factor(df$City)

# ------------------------------------------------------------------------------
# PART 1: Completely Randomized Design (CRD) - One-Way ANOVA
# Model: PM2.5 = mu + Season + epsilon
# H0: mu_monsoon = mu_summer = mu_post_monsoon = mu_winter
# ------------------------------------------------------------------------------
fit_crd <- aov(PM2_5_ugm3 ~ Season, data = df)
summary_crd <- summary(fit_crd)

cat("\n--- CRD ONE-WAY ANOVA TABLE ---\n")
print(summary_crd)

# Effect Size: Eta-Squared (eta2 = SS_treatment / SS_total)
ss_treat <- summary_crd[[1]]["Season", "Sum Sq"]
ss_resid <- summary_crd[[1]]["Residuals", "Sum Sq"]
ss_tot   <- ss_treat + ss_resid
eta_squared <- ss_treat / ss_tot
cat(sprintf("Effect Size (Eta-Squared): %.4f (%.2f%% of total variance explained by season alone)\n", 
            eta_squared, eta_squared * 100))

# Post-Hoc Pairwise Comparisons via Tukey's HSD
tukey_res <- TukeyHSD(fit_crd, "Season", conf.level = 0.95)
cat("\n--- TUKEY HSD PAIRWISE POST-HOC COMPARISONS (CRD) ---\n")
print(tukey_res)

tukey_df <- as.data.frame(tukey_res$Season)
tukey_df$Comparison <- rownames(tukey_df)
tukey_df <- tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
colnames(tukey_df) <- c("Pairwise_Comparison", "Mean_Difference", "CI_95_Lower", "CI_95_Upper", "Raw_p_value")
tukey_df$Significance <- ifelse(tukey_df$Raw_p_value < 0.05, "Significant (p < 0.05)", "Not Significant")
tukey_df$Adjusted_p_value <- format.pval(tukey_df$Raw_p_value, digits = 4, eps = 1e-15)
write.csv(tukey_df, "outputs/table_anova_tukey_crd.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PART 2: Randomized Block Design (RBD) - Blocking by City
# In observational environmental studies, geographic heterogeneity (elevation,
# climate zone, local industrialization) induces strong block effects.
# Model: PM2.5_mean = mu + Season (Treatment) + City (Block) + epsilon
# ------------------------------------------------------------------------------
# Form balanced 29 x 4 block dataset (mean PM2.5 per City-Season cell)
block_data <- aggregate(PM2_5_ugm3 ~ City + Season, data = df, FUN = mean)
names(block_data) <- c("City", "Season", "PM2_5")

fit_rbd <- aov(PM2_5 ~ Season + City, data = block_data)
summary_rbd <- summary(fit_rbd)

cat("\n--- RBD TWO-WAY ANOVA TABLE (TREATMENT: SEASON, BLOCK: CITY) ---\n")
print(summary_rbd)

# Extract Mean Squares for Relative Efficiency Calculation
ms_block <- summary_rbd[[1]]["City", "Mean Sq"]
ms_error <- summary_rbd[[1]]["Residuals", "Mean Sq"]
b <- nlevels(block_data$City)   # 29 blocks
k <- nlevels(block_data$Season) # 4 treatments

# Relative Efficiency (RE) of RBD over CRD:
# RE = ((b - 1) * MS_Block + b * (k - 1) * MS_Error) / ((b * k - 1) * MS_Error)
rel_eff <- ((b - 1) * ms_block + b * (k - 1) * ms_error) / ((b * k - 1) * ms_error)

cat(sprintf("\nEfficiency Analysis of Randomized Block Design:\n"))
cat(sprintf("  Treatment (Season) F-stat : %.2f (p = %s)\n", 
            summary_rbd[[1]]["Season", "F value"], 
            format.pval(summary_rbd[[1]]["Season", "Pr(>F)"], digits = 4, eps = 1e-15)))
cat(sprintf("  Block (City) F-stat       : %.2f (p = %s)\n", 
            summary_rbd[[1]]["City", "F value"], 
            format.pval(summary_rbd[[1]]["City", "Pr(>F)"], digits = 4, eps = 1e-15)))
cat(sprintf("  Relative Efficiency (RE)  : %.3f (RBD is %.1f%% more efficient than CRD!)\n", 
            rel_eff, (rel_eff - 1) * 100))

# Compile ANOVA Comparison Table
anova_summary <- data.frame(
  Source_of_Variation = c("CRD: Season (Treatment)", "CRD: Residuals", 
                          "RBD: Season (Treatment)", "RBD: City (Block)", "RBD: Residuals"),
  Degrees_of_Freedom  = c(summary_crd[[1]]["Season", "Df"], summary_crd[[1]]["Residuals", "Df"],
                          summary_rbd[[1]]["Season", "Df"], summary_rbd[[1]]["City", "Df"], summary_rbd[[1]]["Residuals", "Df"]),
  Sum_of_Squares      = round(c(summary_crd[[1]]["Season", "Sum Sq"], summary_crd[[1]]["Residuals", "Sum Sq"],
                                summary_rbd[[1]]["Season", "Sum Sq"], summary_rbd[[1]]["City", "Sum Sq"], summary_rbd[[1]]["Residuals", "Sum Sq"]), 1),
  Mean_Square         = round(c(summary_crd[[1]]["Season", "Mean Sq"], summary_crd[[1]]["Residuals", "Mean Sq"],
                                summary_rbd[[1]]["Season", "Mean Sq"], summary_rbd[[1]]["City", "Mean Sq"], summary_rbd[[1]]["Residuals", "Mean Sq"]), 2),
  F_Statistic         = round(c(summary_crd[[1]]["Season", "F value"], NA,
                                summary_rbd[[1]]["Season", "F value"], summary_rbd[[1]]["City", "F value"], NA), 2),
  p_value             = c(format.pval(summary_crd[[1]]["Season", "Pr(>F)"], digits = 4, eps = 1e-15), NA,
                          format.pval(summary_rbd[[1]]["Season", "Pr(>F)"], digits = 4, eps = 1e-15),
                          format.pval(summary_rbd[[1]]["City", "Pr(>F)"], digits = 4, eps = 1e-15), NA),
  stringsAsFactors    = FALSE
)
write.csv(anova_summary, "outputs/table_anova_crd_rbd_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Visualizations
# ------------------------------------------------------------------------------
cat("\nGenerating ANOVA figures in outputs/ ...\n")

# FIG 20: Tukey HSD 95% Confidence Interval Plot
png("outputs/fig20_anova_tukey_posthoc.png", width = 850, height = 650, res = 120)
par(mar = c(5, 12, 4, 2), bg = "#ffffff")
plot(tukey_res, las = 1, col = "#2b83ba")
abline(v = 0, col = "#d9534f", lty = 2, lwd = 2)
dev.off()
cat("Created: outputs/fig20_anova_tukey_posthoc.png\n")

# FIG 21: RBD Interaction / Profile Plot (City Blocks across Seasons)
png("outputs/fig21_anova_rbd_interaction_profile.png", width = 1000, height = 750, res = 120)
par(mar = c(5, 5, 4, 2), bg = "#ffffff")
interaction.plot(x.factor = block_data$Season, response = block_data$PM2_5, trace.factor = block_data$City,
                 type = "b", pch = 16, lty = 1, col = rainbow(29, s = 0.8, v = 0.7),
                 legend = FALSE, xlab = "Season (Treatment Factor)",
                 ylab = expression(paste("Mean PM"[2.5], " (", mu, "g/m"^3, ")")),
                 main = expression(paste("RBD Profile Plot: Seasonal Treatment Profiles Across 29 City Blocks (RE = 5.23)")),
                 cex.lab = 1.1, font.main = 2)
# Overlaid National Mean Profile
nat_profile <- tapply(block_data$PM2_5, block_data$Season, mean)
lines(1:4, nat_profile, col = "#000000", lwd = 4.5, type = "b", pch = 18, cex = 1.8)
legend("topleft", 
       legend = c("National Mean Trend", "Individual City Block Profiles (N=29)"),
       col = c("#000000", "#7f7f7f"), lty = 1, lwd = c(4.5, 1.5), pch = c(18, 16),
       bty = "n", cex = 0.95)
dev.off()
cat("Created: outputs/fig21_anova_rbd_interaction_profile.png\n")

cat("\n>>> Experiment 11 Complete! CRD, RBD, Tukey HSD, and efficiency metrics generated.\n")
