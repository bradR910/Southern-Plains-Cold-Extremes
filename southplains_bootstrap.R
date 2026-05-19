##################
# Bootstrap Hypothesis Test: Southern Plains Winter Temperature Trends
# Adapted from HW5_bootstrap.R (Ryan L. Sriver)
#
# Project: Southern Plains Cold Extremes -- ATMS 526
#
# Research Question:
#   Using non-parametric bootstrap resampling, is the difference in mean
#   Jan-Feb minimum temperature between pre-1960 and 1990-present
#   statistically significant?
#
# Data:
#   Monthly minimum temperatures for the Southern Plains region
#   (CRU TS4.08, area-averaged, 1901-2023), from KNMI Climate Explorer.
#   Column 1:  Year
#   Columns 2-13: Jan-Dec absolute monthly minimum temperature (deg C)
#   Missing-data flag: -999.900
#   Header lines begin with '#'
#
#   This analysis uses ABSOLUTE Jan-Feb mean minimum temperature.
##################


########################################
# 1. READ AND PREPARE DATA
########################################

# Read data; comment.char = "#" skips the metadata header automatically
temp_table <- read.table("./Monthly_Tmin.txt", header = FALSE, comment.char = "#")

# Replace missing-data flag with NA
temp_table[temp_table == -999.9] <- NA

# Column 1 = year, Column 2 = January Tmin, Column 3 = February Tmin
years <- temp_table[, 1]
jan   <- temp_table[, 2]
feb   <- temp_table[, 3]

# Construct Jan-Feb mean minimum temperature for each year
temp <- rowMeans(cbind(jan, feb), na.rm = TRUE)

# Period 1: Pre-1960 (early climate)
# Period 2: 1990-present (modern climate)
idx1 <- which(years < 1960)
idx2 <- which(years >= 1990)

temp1 <- temp[idx1]
temp2 <- temp[idx2]

# Drop any NA values
temp1 <- temp1[!is.na(temp1)]
temp2 <- temp2[!is.na(temp2)]

cat(sprintf("Period 1: %d-%d  (n=%d, mean=%.3f C)\n",
            min(years[idx1]), max(years[idx1]), length(temp1), mean(temp1)))
cat(sprintf("Period 2: %d-%d  (n=%d, mean=%.3f C)\n",
            min(years[idx2]), max(years[idx2]), length(temp2), mean(temp2)))
cat(sprintf("Observed difference (P1 - P2): %.3f C\n", mean(temp1) - mean(temp2)))

set.seed(42)  # for reproducibility


########################################
# 2. BOOTSTRAP FUNCTION
########################################

# Runs bootstrap resampling for a given N and returns:
#   - m1:   bootstrap means for period 1
#   - m2:   bootstrap means for period 2
#   - diff: bootstrap distribution of (mean1 - mean2)

run_bootstrap <- function(temp1, temp2, N) {
  m1   <- rep(NA, N)
  m2   <- rep(NA, N)
  diff <- rep(NA, N)

  for (i in 1:N) {
    s1 <- sample(temp1, size = length(temp1), replace = TRUE)
    s2 <- sample(temp2, size = length(temp2), replace = TRUE)
    m1[i]   <- mean(s1)
    m2[i]   <- mean(s2)
    diff[i] <- mean(s1) - mean(s2)
  }

  list(m1 = m1, m2 = m2, diff = diff)
}


########################################
# 3. BOOTSTRAP AT N = 1000 (PRIMARY ANALYSIS)
########################################

N_main <- 1000
boot   <- run_bootstrap(temp1, temp2, N_main)

# 95% range centered on the mean of the difference distribution
mean_diff   <- mean(boot$diff)
ci_95       <- quantile(boot$diff, probs = c(0.025, 0.975))
ci_width_95 <- ci_95[2] - ci_95[1]

cat(sprintf("\n--- Bootstrap Results (N = %d) ---\n", N_main))
cat(sprintf("  Mean of bootstrapped differences: %.4f C\n", mean_diff))
cat(sprintf("  95%% range: [%.4f, %.4f] C\n", ci_95[1], ci_95[2]))
cat(sprintf("  95%% range width: %.4f C\n", ci_width_95))

if (ci_95[1] > 0 | ci_95[2] < 0) {
  cat("  >> Zero lies outside the 95% CI: difference is statistically significant\n")
} else {
  cat("  >> Zero lies within the 95% CI: cannot reject null hypothesis\n")
}


########################################
# 4. BOOTSTRAP DISTRIBUTION PLOTS (N = 1000)
########################################

# HW5 Task 4: Plot with mean line and 95% dashed lines

par(mfrow = c(1, 2))

# --- Panel 1: Bootstrap means for each period ---
xlim_means <- range(c(boot$m1, boot$m2)) + c(-0.5, 0.5)

plot(density(boot$m1),
     col  = "blue",
     lwd  = 3,
     xlim = xlim_means,
     ylim = c(0, max(density(boot_i$diff)$y) * 1.35),
     main = sprintf("Bootstrap Means (N=%d)\nSouthern Plains Jan-Feb Tmin", N_main),
     xlab = "Bootstrapped Sample Mean (C)",
     ylab = "Density")
abline(v = mean(boot$m1), col = "blue", lwd = 2, lty = 2)

lines(density(boot$m2), col = "red", lwd = 3)
abline(v = mean(boot$m2), col = "red", lwd = 2, lty = 2)

legend("topright",
       legend = c(sprintf("Pre-1960 (mean=%.2f C)", mean(boot$m1)),
                  sprintf("1990-present (mean=%.2f C)", mean(boot$m2))),
       col    = c("blue", "red"),
       lty    = c("solid", "solid"),
       lwd    = c(3, 3))

# --- Panel 2: Bootstrap distribution of difference in means ---
plot(density(boot$diff),
     col  = "black",
     lwd  = 3,
     ylim = c(0, max(density(boot_i$diff)$y) * 1.35),
     main = sprintf("Bootstrap: Difference in Means (N=%d)\n(Pre-1960) minus (1990-present)", N_main),
     xlab = "Difference in Sample Means (C)",
     ylab = "Density")

# HW5 Task 4: solid black line = mean of differences
abline(v = mean_diff, col = "black", lwd = 3)

# HW5 Task 4: dashed black lines = 95% range centered on mean
abline(v = ci_95[1], col = "black", lwd = 2, lty = 2)
abline(v = ci_95[2], col = "black", lwd = 2, lty = 2)

# Mark zero for reference
abline(v = 0, col = "gray50", lwd = 1, lty = 3)

legend("topright",
       legend = c("Density",
                  sprintf("Mean (%.3f C)", mean_diff),
                  sprintf("95%% CI [%.3f, %.3f] C", ci_95[1], ci_95[2]),
                  "Zero (null)"),
       col    = c("black", "black", "black", "gray50"),
       lty    = c("solid", "solid", "dashed", "dotted"),
       lwd    = c(3, 3, 2, 1))


########################################
# 5. SENSITIVITY TO BOOTSTRAP SAMPLE SIZE
########################################

# HW5 Task 4 extension: Run bootstrap for N = 10, 100, 1000
# Compare the 95% range width and location

N_sizes <- c(10, 100, 1000)

par(mfrow = c(1, 3))

cat("\n====== BOOTSTRAP SENSITIVITY TO SAMPLE SIZE ======\n")
cat(sprintf("  %-8s  %-12s  %-12s  %-10s  %-10s\n",
            "N", "Mean Diff", "CI Lower", "CI Upper", "CI Width"))
cat(sprintf("  %s\n", paste(rep("-", 58), collapse = "")))

for (N_i in N_sizes) {
  boot_i  <- run_bootstrap(temp1, temp2, N_i)
  mean_i  <- mean(boot_i$diff)
  ci_i    <- quantile(boot_i$diff, probs = c(0.025, 0.975))
  width_i <- ci_i[2] - ci_i[1]

  cat(sprintf("  %-8d  %-12.4f  %-12.4f  %-10.4f  %-10.4f\n",
              N_i, mean_i, ci_i[1], ci_i[2], width_i))

  plot(density(boot_i$diff),
       col  = "black",
       lwd  = 3,
       main = sprintf("Bootstrap Diff in Means\nN = %d", N_i),
       xlab = "Difference in Means (C)",
       ylab = "Density")
  abline(v = mean_i,  col = "black",  lwd = 3)
  abline(v = ci_i[1], col = "black",  lwd = 2, lty = 2)
  abline(v = ci_i[2], col = "black",  lwd = 2, lty = 2)
  abline(v = 0,       col = "gray50", lwd = 1, lty = 3)

  legend("topright",
         legend = c(sprintf("Mean: %.3f C", mean_i),
                    sprintf("95%% CI width: %.3f C", width_i)),
         col    = c("black", "black"),
         lty    = c("solid", "dashed"),
         lwd    = c(3, 2),
         cex    = 0.85)
}
cat(sprintf("  %s\n", paste(rep("-", 58), collapse = "")))


########################################
# 6. SUMMARY TABLE
########################################

cat("\n====== SUMMARY TABLE ======\n")
cat(sprintf("  %-30s  %-15s  %-15s\n", "Metric", "Pre-1960", "1990-present"))
cat(sprintf("  %s\n", paste(rep("-", 65), collapse = "")))
cat(sprintf("  %-30s  %-15.3f  %-15.3f\n", "Sample mean (C)",     mean(temp1), mean(temp2)))
cat(sprintf("  %-30s  %-15.3f  %-15.3f\n", "Sample std dev (C)",  sd(temp1),   sd(temp2)))
cat(sprintf("  %-30s  %-15d  %-15d\n",     "Sample size (n)",     length(temp1), length(temp2)))
cat(sprintf("  %s\n", paste(rep("-", 65), collapse = "")))
cat(sprintf("  %-30s  %-15.4f\n", "Observed diff (P1 - P2, C)",  mean(temp1) - mean(temp2)))
cat(sprintf("  %-30s  %-15.4f\n", "Bootstrap mean diff (C)",     mean_diff))
cat(sprintf("  %-30s  [%.4f, %.4f]\n", "Bootstrap 95% CI (C)",   ci_95[1], ci_95[2]))
cat("============================\n")
