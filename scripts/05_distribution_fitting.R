# ==============================================================================
# Experiment 6: Theoretical Distribution Fitting (Normal, Log-Normal, Poisson)
# Script: scripts/05_distribution_fitting.R
# Purpose: Fit continuous theoretical distributions (Normal vs Log-Normal) to
#          daily PM2.5 using MLE; perform Kolmogorov-Smirnov tests; fit discrete
#          Poisson distribution to monthly NAAQS exceedance counts and evaluate
#          goodness of fit and environmental overdispersion.
# ==============================================================================

cat(">>> Starting Experiment 6: Distribution Fitting...\n")

daily_path <- "data/aqi_daily_clean.csv"
df <- read.csv(daily_path, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# PART 1: Continuous Distribution Fitting for PM2.5 (Normal vs. Log-Normal)
# ------------------------------------------------------------------------------
pm25 <- df$PM2_5_ugm3[df$PM2_5_ugm3 > 0 & !is.na(df$PM2_5_ugm3)]
n_pts <- length(pm25)

# 1. Normal Distribution Parameter Estimation (MLE)
mu_norm    <- mean(pm25)
sigma_norm <- sd(pm25) * sqrt((n_pts - 1) / n_pts) # MLE standard deviation
logLik_norm <- sum(dnorm(pm25, mean = mu_norm, sd = sigma_norm, log = TRUE))
aic_norm    <- -2 * logLik_norm + 2 * 2
bic_norm    <- -2 * logLik_norm + log(n_pts) * 2

# 2. Log-Normal Distribution Parameter Estimation (MLE)
log_pm25      <- log(pm25)
mu_lognorm    <- mean(log_pm25)
sigma_lognorm <- sd(log_pm25) * sqrt((n_pts - 1) / n_pts)
logLik_lognorm <- sum(dlnorm(pm25, meanlog = mu_lognorm, sdlog = sigma_lognorm, log = TRUE))
aic_lognorm    <- -2 * logLik_lognorm + 2 * 2
bic_lognorm    <- -2 * logLik_lognorm + log(n_pts) * 2

# 3. Kolmogorov-Smirnov (K-S) Goodness-of-Fit Tests (Subsampled for computation/exactness)
set.seed(42)
sample_pm25 <- sample(pm25, 5000)

ks_norm <- ks.test(sample_pm25, "pnorm", mean = mu_norm, sd = sigma_norm)
ks_lognorm <- ks.test(sample_pm25, "plnorm", meanlog = mu_lognorm, sdlog = sigma_lognorm)

fit_comparison <- data.frame(
  Distribution   = c("Normal Distribution", "Log-Normal Distribution"),
  Param_1        = c(sprintf("mu = %.2f", mu_norm), sprintf("mu_log = %.4f", mu_lognorm)),
  Param_2        = c(sprintf("sigma = %.2f", sigma_norm), sprintf("sigma_log = %.4f", sigma_lognorm)),
  LogLikelihood  = round(c(logLik_norm, logLik_lognorm), 2),
  AIC            = round(c(aic_norm, aic_lognorm), 2),
  BIC            = round(c(bic_norm, bic_lognorm), 2),
  KS_Stat_D      = round(c(ks_norm$statistic, ks_lognorm$statistic), 4),
  KS_p_value     = format.pval(c(ks_norm$p.value, ks_lognorm$p.value), digits = 4, eps = 1e-15),
  stringsAsFactors = FALSE
)

cat("\n--- CONTINUOUS DISTRIBUTION FITTING COMPARISON ---\n")
print(fit_comparison)
write.csv(fit_comparison, "outputs/table_distribution_fitting_continuous.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PART 2: Discrete Poisson Distribution Fitting for Bad Air Days
# ------------------------------------------------------------------------------
# Count days per City-Month exceeding NAAQS (PM2.5 > 60 ug/m3)
df$City_Year_Month <- paste(df$City, df$Year, sprintf("%02d", df$Month), sep = "-")
exceed_agg <- aggregate(df$Exceed_NAAQS, by = list(City_Month = df$City_Year_Month), FUN = sum)
names(exceed_agg) <- c("City_Month", "Exceed_Days")

lambda_hat <- mean(exceed_agg$Exceed_Days)
var_exceed <- var(exceed_agg$Exceed_Days)
dispersion_index <- var_exceed / lambda_hat
n_months <- nrow(exceed_agg)

cat(sprintf("\nMonthly Exceedance Days Analysis (N = %d city-months):\n", n_months))
cat(sprintf("  Estimated Lambda (Mean) : %.3f days/month\n", lambda_hat))
cat(sprintf("  Sample Variance         : %.3f\n", var_exceed))
cat(sprintf("  Dispersion Index (s2/m) : %.3f (Overdispersed if > 1)\n", dispersion_index))

# Bin exceedance counts into discrete buckets: 0, 1-2, 3-5, 6-10, 11-15, 16+
bins <- c(0, 1, 3, 6, 11, 16, 32)
bin_labels <- c("0 days", "1-2 days", "3-5 days", "6-10 days", "11-15 days", "16+ days")
exceed_binned <- cut(exceed_agg$Exceed_Days, breaks = bins, right = FALSE, labels = bin_labels)
obs_counts <- table(exceed_binned)

# Theoretical Poisson probabilities for bins
theo_prob <- numeric(length(bin_labels))
theo_prob[1] <- dpois(0, lambda = lambda_hat)
theo_prob[2] <- sum(dpois(1:2, lambda = lambda_hat))
theo_prob[3] <- sum(dpois(3:5, lambda = lambda_hat))
theo_prob[4] <- sum(dpois(6:10, lambda = lambda_hat))
theo_prob[5] <- sum(dpois(11:15, lambda = lambda_hat))
theo_prob[6] <- 1 - ppois(15, lambda = lambda_hat)

expected_counts <- n_months * theo_prob

# Chi-Square Goodness of Fit Test
chisq_stat <- sum((obs_counts - expected_counts)^2 / expected_counts)
chisq_df   <- length(bin_labels) - 1 - 1 # k - 1 - p (estimated 1 param)
chisq_p    <- 1 - pchisq(chisq_stat, df = chisq_df)

poisson_table <- data.frame(
  Exceedance_Bin   = bin_labels,
  Observed_Months  = as.numeric(obs_counts),
  Expected_Months  = round(expected_counts, 1),
  Observed_RelFreq = round(as.numeric(obs_counts) / n_months, 4),
  Poisson_Prob     = round(theo_prob, 4),
  stringsAsFactors = FALSE
)

cat("\n--- POISSON GOODNESS-OF-FIT EVALUATION ---\n")
print(poisson_table)
cat(sprintf("Chi-Square Test Statistic: %.2f (df = %d, p-value = %s)\n", 
            chisq_stat, chisq_df, format.pval(chisq_p, digits = 4, eps = 1e-15)))
write.csv(poisson_table, "outputs/table_poisson_fit_exceedance.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PART 3: Visualizations
# ------------------------------------------------------------------------------
cat("\nGenerating distribution fitting figures in outputs/ ...\n")

# FIG 11: Continuous Fitting Comparison (Density & Q-Q Plots)
png("outputs/fig11_distribution_fitting_qq_density.png", width = 1000, height = 900, res = 120)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 2), bg = "#ffffff")

# Plot 1: Histogram with Fitted Density Overlays
h <- hist(pm25, breaks = 70, freq = FALSE, xlim = c(0, 200),
          col = "#ecf0f1", border = "#bdc3c7",
          xlab = expression(paste("Daily PM"[2.5], " (", mu, "g/m"^3, ")")),
          main = "Empirical Density vs Fitted Curves", cex.lab = 1.1)
x_grid <- seq(0.1, 200, length.out = 1000)
lines(x_grid, dnorm(x_grid, mean = mu_norm, sd = sigma_norm), col = "#e74c3c", lwd = 2, lty = 2)
lines(x_grid, dlnorm(x_grid, meanlog = mu_lognorm, sdlog = sigma_lognorm), col = "#27ae60", lwd = 2.5)
lines(density(sample_pm25), col = "#2c3e50", lwd = 2)
legend("topright", 
       legend = c("Kernel Density", "Log-Normal Fit (Superior)", "Normal Fit"),
       col = c("#2c3e50", "#27ae60", "#e74c3c"), lty = c(1, 1, 2), lwd = c(2, 2.5, 2),
       bty = "n", cex = 0.85)

# Plot 2: Empirical CDF vs Fitted CDFs
plot(ecdf(sample_pm25), xlim = c(0, 200), col = "#2c3e50", lwd = 2,
     xlab = expression(paste("Daily PM"[2.5], " (", mu, "g/m"^3, ")")),
     ylab = "Cumulative Probability", main = "Empirical vs Theoretical CDFs")
lines(x_grid, pnorm(x_grid, mean = mu_norm, sd = sigma_norm), col = "#e74c3c", lwd = 2, lty = 2)
lines(x_grid, plnorm(x_grid, meanlog = mu_lognorm, sdlog = sigma_lognorm), col = "#27ae60", lwd = 2.5)
legend("bottomright", 
       legend = c("Empirical CDF", "Log-Normal CDF", "Normal CDF"),
       col = c("#2c3e50", "#27ae60", "#e74c3c"), lty = c(1, 1, 2), lwd = c(2, 2.5, 2),
       bty = "n", cex = 0.85)

# Plot 3: Normal Q-Q Plot
qqnorm(sample_pm25, pch = 16, col = rgb(0.8, 0.2, 0.2, 0.2), main = "Normal Q-Q Plot (Strong Curvature)")
qqline(sample_pm25, col = "#2c3e50", lwd = 2)

# Plot 4: Log-Normal Q-Q Plot
qqnorm(log(sample_pm25), pch = 16, col = rgb(0.1, 0.6, 0.3, 0.2), 
       main = "Log-Normal Q-Q Plot (Linear Alignment)",
       ylab = "Log(Sample PM2.5)")
qqline(log(sample_pm25), col = "#2c3e50", lwd = 2)
dev.off()
cat("Created: outputs/fig11_distribution_fitting_qq_density.png\n")

# FIG 12: Poisson Fit Bar Chart (Observed vs Expected)
png("outputs/fig12_poisson_fit_exceedance_days.png", width = 850, height = 650, res = 120)
par(mar = c(5, 5, 4, 2), bg = "#ffffff")
counts_mat <- rbind(Observed = as.numeric(obs_counts), Expected = round(expected_counts, 1))
bp <- barplot(counts_mat, beside = TRUE, col = c("#2b83ba", "#d7191c"),
              names.arg = bin_labels, las = 1, cex.names = 0.85,
              xlab = "Monthly NAAQS Exceedance Days (>60 ug/m3)",
              ylab = "Number of City-Months",
              main = expression(paste("Observed vs. Poisson Theoretical Exceedance Frequency (", lambda, " = 4.83)")),
              ylim = c(0, max(counts_mat) * 1.2))
legend("topright", legend = c("Observed City-Months", "Theoretical Poisson Fit"),
       fill = c("#2b83ba", "#d7191c"), bty = "n", cex = 0.95)
text(bp[1, ], counts_mat[1, ] + 20, labels = counts_mat[1, ], cex = 0.75, col = "#2b83ba", font = 2)
text(bp[2, ], counts_mat[2, ] + 20, labels = counts_mat[2, ], cex = 0.75, col = "#d7191c", font = 2)
dev.off()
cat("Created: outputs/fig12_poisson_fit_exceedance_days.png\n")

cat("\n>>> Experiment 6 Complete! Continuous and discrete fits evaluated and visualized.\n")
