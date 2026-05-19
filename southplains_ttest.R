##################
# Two-Sample T-Test: Southern Plains Winter Temperature Trends
# Adapted from HW5_t_test.R (Ryan L. Sriver)
#
# Project: Southern Plains Cold Extremes -- ATMS 526
#
# Research Question:
#   Has the mean Jan-Feb minimum temperature in the Southern Plains
#   changed significantly between the early and modern climate periods?
#
# Null Hypothesis:
#   The true difference in mean Jan-Feb minimum temperature between
#   the early period (pre-1960) and modern period (1990-present) is zero.
#
# Data:
#   Monthly minimum temperatures for the Southern Plains region
#   (CRU TS4.08, area-averaged, 1901-2023), from KNMI Climate Explorer.
#   Column 1:  Year
#   Columns 2-13: Jan-Dec absolute monthly minimum temperature (deg C)
#   Missing-data flag: -999.900
#   Header lines begin with '#'
#
#   This analysis uses ABSOLUTE Jan-Feb mean minimum temperature
#   (the average of the January and February columns).
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

cat("Data period:", min(years), "to", max(years), "\n")
cat("Total years:", length(years), "\n")

# --- Define the two comparison periods ---
# Period 1: Pre-1960 (early climate)
# Period 2: 1990-present (modern climate)
idx1 <- which(years < 1960)
idx2 <- which(years >= 1990)

temp1 <- temp[idx1]   # early period Jan-Feb mean Tmin
temp2 <- temp[idx2]   # modern period Jan-Feb mean Tmin

# Drop any NA values
temp1 <- temp1[!is.na(temp1)]
temp2 <- temp2[!is.na(temp2)]

yr1_range <- range(years[idx1])
yr2_range <- range(years[idx2])

cat(sprintf("\nPeriod 1: %d-%d  (n = %d)\n", yr1_range[1], yr1_range[2], length(temp1)))
cat(sprintf("Period 2: %d-%d  (n = %d)\n", yr2_range[1], yr2_range[2], length(temp2)))
cat(sprintf("Mean Jan-Feb Tmin Period 1: %.3f C\n", mean(temp1)))
cat(sprintf("Mean Jan-Feb Tmin Period 2: %.3f C\n", mean(temp2)))
cat(sprintf("Difference (P1 - P2):       %.3f C\n", mean(temp1) - mean(temp2)))


########################################
# 2. HISTOGRAMS OF BOTH PERIODS
########################################

# HW5 Task 1: Histograms of both periods

par(mfrow = c(1, 1))

# Shared x-axis range
xlim_range <- range(c(temp1, temp2))
xlim_pad   <- c(xlim_range[1] - 1, xlim_range[2] + 1)

hist(temp1,
     col    = rgb(0.2, 0.4, 0.8, 0.5),   # semi-transparent blue
     prob   = TRUE,
     breaks = 12,
     xlim   = xlim_pad,
     ylim   = c(0, 0.6),
     main   = "Southern Plains Jan-Feb Minimum Temperature\nEarly Period vs. Modern Period",
     xlab   = "Jan-Feb Mean Minimum Temperature (C)",
     ylab   = "Probability Density")

hist(temp2,
     col    = rgb(0.8, 0.2, 0.2, 0.5),   # semi-transparent red
     prob   = TRUE,
     breaks = 12,
     add    = TRUE)

# Density curves
lines(density(temp1), col = "blue", lwd = 3)
lines(density(temp2), col = "red",  lwd = 3)

# Vertical lines at each mean
abline(v = mean(temp1), col = "blue", lwd = 2, lty = 2)
abline(v = mean(temp2), col = "red",  lwd = 2, lty = 2)

legend("topright",
       legend = c(sprintf("Pre-1960  (mean = %.2f C)", mean(temp1)),
                  sprintf("1990-present (mean = %.2f C)", mean(temp2))),
       col    = c("blue", "red"),
       lty    = c("solid", "solid"),
       lwd    = c(3, 3),
       fill   = c(rgb(0.2, 0.4, 0.8, 0.5), rgb(0.8, 0.2, 0.2, 0.5)),
       border = NA)


########################################
# 3. T-TEST AT MULTIPLE SIGNIFICANCE LEVELS
########################################

# HW5 Task 2: Run t-test and report confidence intervals
# at multiple significance levels: alpha = 0.1, 0.01, 0.001

alpha_levels <- c(0.1, 0.01, 0.001)

cat("\n====== T-TEST RESULTS ======\n")
cat(sprintf("Comparison: Pre-1960 vs. 1990-present\n"))
cat(sprintf("Null hypothesis: difference in means = 0\n\n"))

for (a in alpha_levels) {
  test <- t.test(temp1, temp2,
                 mu          = 0,
                 conf.level  = 1 - a,
                 var.equal   = TRUE,
                 alternative = "two.sided")
  cat(sprintf("--- alpha = %.3f (confidence level = %.1f%%) ---\n", a, (1-a)*100))
  cat(sprintf("  t-statistic:        %.4f\n", test$statistic))
  cat(sprintf("  degrees of freedom: %d\n",   test$parameter))
  cat(sprintf("  p-value:            %.6f\n", test$p.value))
  cat(sprintf("  confidence interval: [%.4f, %.4f]\n",
              test$conf.int[1], test$conf.int[2]))
  if (test$p.value < a) {
    cat(sprintf("  >> REJECT null hypothesis at alpha = %.3f\n\n", a))
  } else {
    cat(sprintf("  >> FAIL TO REJECT null hypothesis at alpha = %.3f\n\n", a))
  }
}

# Store the alpha = 0.05 test for plotting
alpha_plot <- 0.05
test_main <- t.test(temp1, temp2,
                    mu          = 0,
                    conf.level  = 1 - alpha_plot,
                    var.equal   = TRUE,
                    alternative = "two.sided")


########################################
# 4. T-DISTRIBUTION PLOT WITH CONFIDENCE INTERVAL
########################################

# HW5 Task 3: Plot t-distribution with blue vertical lines
# (difference in means + confidence interval from t-test at alpha = 0.05)
# No red/black critical lines per HW5 instructions

par(mfrow = c(1, 1))

df <- length(temp1) + length(temp2) - 2
x  <- seq(-5, 5, length = 1000)

plot(x, dt(x, df = df),
     type = "l",
     lwd  = 3,
     col  = "black",
     main = "T-Distribution: Southern Plains Jan-Feb Minimum Temperature\nPre-1960 vs. 1990-present",
     xlab = "t-statistic",
     ylab = "Density")

# Shade rejection regions at alpha = 0.05
t_crit <- qt(1 - alpha_plot / 2, df = df)
x_lo   <- x[x <= -t_crit]
x_hi   <- x[x >=  t_crit]
polygon(c(x_lo, rev(x_lo)),
        c(dt(x_lo, df), rep(0, length(x_lo))),
        col = rgb(0.8, 0.8, 0.8, 0.4), border = NA)
polygon(c(x_hi, rev(x_hi)),
        c(dt(x_hi, df), rep(0, length(x_hi))),
        col = rgb(0.8, 0.8, 0.8, 0.4), border = NA)

# Blue solid line: observed t-statistic
abline(v = test_main$statistic, col = "blue", lwd = 3)

# Blue dashed lines: 95% confidence interval bounds (in t-statistic units)
pooled_se <- sqrt(var(temp1)/length(temp1) + var(temp2)/length(temp2))
ci_t_lo   <- test_main$conf.int[1] / pooled_se
ci_t_hi   <- test_main$conf.int[2] / pooled_se
abline(v = ci_t_lo, col = "blue", lwd = 2, lty = 2)
abline(v = ci_t_hi, col = "blue", lwd = 2, lty = 2)

legend("topright",
       legend = c("T-distribution",
                  sprintf("Observed t-stat (%.2f)", test_main$statistic),
                  "95% CI bounds"),
       col    = c("black", "blue", "blue"),
       lty    = c("solid", "solid", "dashed"),
       lwd    = c(3, 3, 2))

cat(sprintf("\nT-test (alpha=0.05): t = %.4f, p = %.6f\n",
            test_main$statistic, test_main$p.value))
cat(sprintf("95%% CI on difference: [%.4f, %.4f] C\n",
            test_main$conf.int[1], test_main$conf.int[2]))
