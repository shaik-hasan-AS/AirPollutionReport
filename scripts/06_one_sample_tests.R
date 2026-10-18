# ==============================================================================
# Experiment 7: One-Sample Hypothesis Tests (Means and Proportions)
# Script: scripts/06_one_sample_tests.R
# Purpose: Perform one-sample t-tests comparing city and national PM2.5 means
#          against regulatory standards (India NAAQS: 60 ug/m3, WHO: 15 ug/m3);
#          perform one-sample proportion tests for exceedance rates; generate
#          publication-grade confidence interval forest plots.
# ==============================================================================

cat(">>> Starting Experiment 7: One-Sample Hypothesis Tests...\n")

daily_path <- "data/aqi_daily_clean.csv"
df <- read.csv(daily_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. National Level One-Sample t-Tests
# ------------------------------------------------------------------------------
# H0: mu = 60 ug/m3 (India 24-hr NAAQS) vs H1: mu != 60 ug/m3
t_nat_naaqs <- t.test(df$PM2_5_ugm3, mu = 60)

# H0: mu = 15 ug/m3 (WHO 2021 Guideline) vs H1: mu > 15 ug/m3
t_nat_who   <- t.test(df$PM2_5_ugm3, mu = 15, alternative = "greater")

cat("\n--- NATIONAL ONE-SAMPLE t-TEST: AGAINST INDIA NAAQS (60 ug/m3) ---\n")
print(t_nat_naaqs)

cat("\n--- NATIONAL ONE-SAMPLE t-TEST: AGAINST WHO GUIDELINE (15 ug/m3) ---\n")
print(t_nat_who)

# ------------------------------------------------------------------------------
# 2. City-by-City One-Sample t-Tests Against NAAQS (60 ug/m3)
# ------------------------------------------------------------------------------
cities <- sort(unique(df$City))

city_ttest_list <- lapply(cities, function(cty) {
  x <- df$PM2_5_ugm3[df$City == cty]
  n <- length(x)
  t_res <- t.test(x, mu = 60)
  
  # Decision at alpha = 0.05
  decision <- if (t_res$p.value < 0.05) {
    if (t_res$estimate > 60) "Significantly Exceeds NAAQS" else "Significantly Below NAAQS"
  } else {
    "Not Significantly Different from NAAQS"
  }
  
  data.frame(
    City         = cty,
    N            = n,
    Sample_Mean  = round(mean(x), 2),
    Std_Dev      = round(sd(x), 2),
    Std_Error    = round(sd(x) / sqrt(n), 3),
    t_statistic  = round(t_res$statistic, 2),
    df           = t_res$parameter,
    p_value      = format.pval(t_res$p.value, digits = 4, eps = 1e-15),
    CI_95_Lower  = round(t_res$conf.int[1], 2),
    CI_95_Upper  = round(t_res$conf.int[2], 2),
    Status       = decision,
    stringsAsFactors = FALSE
  )
})

city_ttest_table <- do.call(rbind, city_ttest_list)
city_ttest_table <- city_ttest_table[order(-city_ttest_table$Sample_Mean), ]
rownames(city_ttest_table) <- NULL

cat("\n--- CITY-WISE ONE-SAMPLE t-TEST RESULTS (AGAINST 60 ug/m3) ---\n")
print(head(city_ttest_table, 12))
write.csv(city_ttest_table, "outputs/table_one_sample_ttests_by_city.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. City-by-City One-Sample Proportion Tests
# H0: p <= 0.10 (Less than 10% of days exceed NAAQS) vs H1: p > 0.10
# ------------------------------------------------------------------------------
city_prop_list <- lapply(cities, function(cty) {
  sub_df <- df[df$City == cty, ]
  x_exceed <- sum(sub_df$Exceed_NAAQS)
  n_total  <- nrow(sub_df)
  
  p_test <- prop.test(x_exceed, n_total, p = 0.10, alternative = "greater", correct = TRUE)
  
  data.frame(
    City             = cty,
    Total_Days       = n_total,
    Exceed_Days      = x_exceed,
    Observed_Prop    = round(x_exceed / n_total, 4),
    Observed_Pct     = round((x_exceed / n_total) * 100, 2),
    ChiSq_Stat       = round(p_test$statistic, 2),
    p_value          = format.pval(p_test$p.value, digits = 4, eps = 1e-15),
    CI_95_Lower_Pct  = round(p_test$conf.int[1] * 100, 2),
    Reject_H0        = ifelse(p_test$p.value < 0.05, "Yes (>10% bad days)", "No"),
    stringsAsFactors = FALSE
  )
})

city_prop_table <- do.call(rbind, city_prop_list)
city_prop_table <- city_prop_table[order(-city_prop_table$Observed_Prop), ]
rownames(city_prop_table) <- NULL

cat("\n--- CITY-WISE ONE-SAMPLE PROPORTION TESTS (TESTING IF EXCEEDANCE > 10%) ---\n")
print(head(city_prop_table, 12))
write.csv(city_prop_table, "outputs/table_one_sample_proportions_by_city.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Visualizations
# ------------------------------------------------------------------------------
cat("\nGenerating one-sample test figures in outputs/ ...\n")

# FIG 13: 95% Confidence Interval Forest Plot of City Means
png("outputs/fig13_one_sample_ci_forest_plot.png", width = 1000, height = 850, res = 120)
par(mar = c(5, 9, 4, 2), bg = "#ffffff")
ord_df <- city_ttest_table[order(city_ttest_table$Sample_Mean), ]
y_pos <- 1:nrow(ord_df)

pt_cols <- ifelse(ord_df$Status == "Significantly Exceeds NAAQS", "#d7191c",
           ifelse(ord_df$Status == "Significantly Below NAAQS", "#2b83ba", "#fdae61"))

plot(ord_df$Sample_Mean, y_pos, pch = 16, cex = 1.1, col = pt_cols,
     xlim = c(0, max(ord_df$CI_95_Upper) * 1.1), ylim = c(0.5, nrow(ord_df) + 0.5),
     axes = FALSE, xlab = expression(paste("Mean Daily PM"[2.5], " (", mu, "g/m"^3, ") with 95% CI")),
     ylab = "", main = "One-Sample t-Test: 95% Confidence Intervals vs. Regulatory Limits")

axis(1, cex.axis = 0.9)
axis(2, at = y_pos, labels = ord_df$City, las = 1, cex.axis = 0.75)
box()

# Confidence interval whiskers
arrows(ord_df$CI_95_Lower, y_pos, ord_df$CI_95_Upper, y_pos, 
       angle = 90, code = 3, length = 0.03, col = pt_cols, lwd = 1.8)

# Reference Thresholds
abline(v = 60, col = "#c9302c", lwd = 2, lty = 2) # NAAQS
abline(v = 15, col = "#4cae4c", lwd = 2, lty = 3) # WHO

legend("bottomright", 
       legend = c(
         "Significantly Exceeds NAAQS (p < 0.05)",
         "Significantly Complies with NAAQS (p < 0.05)",
         "Statistically Indistinguishable from 60",
         "India 24-hr NAAQS (60 ug/m3)",
         "WHO 24-hr Guideline (15 ug/m3)"
       ),
       col = c("#d7191c", "#2b83ba", "#fdae61", "#c9302c", "#4cae4c"),
       pch = c(16, 16, 16, NA, NA), lty = c(NA, NA, NA, 2, 3), lwd = c(NA, NA, NA, 2, 2),
       bty = "n", cex = 0.85)
dev.off()
cat("Created: outputs/fig13_one_sample_ci_forest_plot.png\n")

# FIG 14: Proportion of Exceedance Days with 95% CI
png("outputs/fig14_one_sample_proportions_ci.png", width = 1000, height = 750, res = 120)
par(mar = c(5, 9, 4, 2), bg = "#ffffff")
ord_p <- city_prop_table[order(city_prop_table$Observed_Prop), ]
bp <- barplot(ord_p$Observed_Pct, horiz = TRUE, names.arg = ord_p$City, las = 1,
              col = ifelse(ord_p$Reject_H0 == "Yes (>10% bad days)", "#d9534f", "#337ab7"),
              border = NA, cex.names = 0.75,
              xlab = "Percentage of Days Exceeding Daily NAAQS (>60 ug/m3)",
              main = "Exceedance Frequency by City (One-Sample Proportion Tests)",
              xlim = c(0, max(ord_p$Observed_Pct) * 1.15))
abline(v = 10, col = "#333333", lwd = 2, lty = 2)
text(ord_p$Observed_Pct + 1.5, bp, labels = sprintf("%.1f%%", ord_p$Observed_Pct), cex = 0.7, pos = 4)
legend("bottomright", 
       legend = c("Exceeds 10% Benchmark (p < 0.05)", "Complies (<= 10% Exceedance)", "10% Regulatory Policy Target"),
       fill = c("#d9534f", "#337ab7", NA), border = NA,
       lty = c(NA, NA, 2), lwd = c(NA, NA, 2), col = c(NA, NA, "#333333"),
       bty = "n", cex = 0.9)
dev.off()
cat("Created: outputs/fig14_one_sample_proportions_ci.png\n")

cat("\n>>> Experiment 7 Complete! National and city-wise one-sample tests finished.\n")
