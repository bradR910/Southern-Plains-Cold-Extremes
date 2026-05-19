########################################
# Analyze DFW Daily Winter Minimum Temperatures
#  - Block Minima / GEV Extreme Value Analysis
#  - Return Periods and Risk Curves
#  - Bootstrap Uncertainty Quantification
#
# Project: Southern Plains Cold Extremes
# Course:  ATMS 526, Risk Analysis
# Adapted from: qq_temp_extremes.R (Ryan Sriver)
#
# NOTE ON APPROACH:
#   GEV theory is built for block maxima. For cold extremes (minima),
#   we negate the data so that the coldest values become the largest.
#   We then fit the GEV to these negated values and negate all return
#   levels back to recover actual temperatures.
#
#   MLE proved numerically unstable for this dataset (~75 block minima).
#   We use L-moments estimation instead, which is more robust for
#   small samples.
########################################

### Load Necessary Libraries
### Install if needed: install.packages(c("extRemes", "MASS"))
library(extRemes)
library(MASS)


########################################
# 1. READ AND PREPARE DATA
########################################

# Read DFW station daily temperature data
# Format: STATION, NAME, DATE (YYYY-MM-DD), TMIN (degrees Fahrenheit)
dfw_table <- read.csv("./DFWMinTemps.csv", header = TRUE, stringsAsFactors = FALSE)

# Parse the DATE column into year and month
dfw_table$DATE  <- as.Date(dfw_table$DATE, format = "%Y-%m-%d")
dfw_table$YEAR  <- as.integer(format(dfw_table$DATE, "%Y"))
dfw_table$MONTH <- as.integer(format(dfw_table$DATE, "%m"))

# Convert TMIN from Fahrenheit to Celsius
dfw_table$TMIN_C <- (as.numeric(dfw_table$TMIN) - 32) * (5/9)

# Remove rows with missing TMIN
dfw_table <- dfw_table[!is.na(dfw_table$TMIN_C), ]

# Extract working vectors
tmin_raw  <- dfw_table$TMIN_C
year_col  <- dfw_table$YEAR
month_col <- dfw_table$MONTH

cat("Data period:", min(year_col), "to", max(year_col), "\n")
cat("Total daily records:", nrow(dfw_table), "\n")

# Filter to January (1) and February (2) only
winter_idx  <- which(month_col %in% c(1, 2) & year_col <= 2025)
tmin_winter <- tmin_raw[winter_idx]
year_winter <- year_col[winter_idx]

# Winter year = calendar year (no Dec crossover needed for JF only)
winter_year <- year_winter

# Unique winter years
all_winter_years <- sort(unique(winter_year))
n_years <- length(all_winter_years)

cat("Number of winter seasons:", n_years, "\n")
cat("Period:", min(all_winter_years), "to", max(all_winter_years), "\n")
cat("Overall Jan-Feb Tmin range:",
    min(tmin_winter, na.rm = TRUE), "to",
    max(tmin_winter, na.rm = TRUE), "C\n")


########################################
# 2. CALCULATE BLOCK MINIMA
########################################

# Annual winter minimum temperature for each season
block_min <- sapply(all_winter_years, function(yr) {
  idx <- which(winter_year == yr)
  min(tmin_winter[idx], na.rm = TRUE)
})

cat("\nBlock minima summary (actual temperatures, C):\n")
print(summary(block_min))

# Remove positive outliers: a positive Jan-Feb block minimum indicates a year
# with no meaningful cold days and is not representative of cold tail behavior
n_before        <- length(block_min)
block_min_clean <- block_min[block_min < 0]
n_removed       <- n_before - length(block_min_clean)
cat(sprintf("\nRemoved %d outlier block minima (positive values)\n", n_removed))
cat(sprintf("Clean sample size: %d years\n", length(block_min_clean)))
cat("Clean block minima range:", range(block_min_clean), "C\n")

# Negate for GEV fitting: coldest value becomes the largest
block_max_neg <- -block_min_clean

# Feb 2021 event: DFW observed minimum was -2 F on Feb 16, 2021 = -18.9 C
feb2021_tmin <- -18.9
feb2021_neg  <- -feb2021_tmin   # = 18.9 in negated space


########################################
# 3. FIT GEV USING L-MOMENTS
########################################

# L-moments estimation is more robust than MLE for small samples (~75 points)
# and avoids optimizer instability. We fit to negated block minima.
print(sort(block_max_neg))
# Check which years have how many Jan-Feb observations
obs_per_year <- sapply(all_winter_years, function(yr) {
  sum(winter_year == yr, na.rm = TRUE)
})

# Print years with suspiciously few observations
cat("Years with fewer than 30 Jan-Feb observations:\n")
print(data.frame(
  year = all_winter_years[obs_per_year < 30],
  n_obs = obs_per_year[obs_per_year < 30]
))

# Also check the actual block minima alongside the year
cat("\nAll block minima by year:\n")
print(data.frame(
  year = all_winter_years,
  n_obs = obs_per_year,
  block_min_F = block_min * (9/5) + 32,  # back to F to make it easier to check
  block_min_C = round(block_min, 2)
))
gev_est <- fevd(block_max_neg, type = "GEV", method = "Lmoments")

# Extract GEV parameters
mu    <- gev_est$results["location"]
sigma <- gev_est$results["scale"]
xi    <- gev_est$results["shape"]

cat("\n--- Best-Estimate GEV Parameters (L-moments, negated block minima) ---\n")
cat(sprintf("  Location (mu):  %.4f\n", mu))
cat(sprintf("  Scale (sigma):  %.4f\n", sigma))
cat(sprintf("  Shape (xi):     %.4f\n", xi))

# Determine GEV family
if (abs(xi) < 0.05) {
  cat("  GEV family: Gumbel (xi ~ 0) -- light tail\n")
} else if (xi > 0) {
  cat("  GEV family: Frechet (xi > 0) -- heavy/unbounded tail\n")
} else {
  cat("  GEV family: Weibull (xi < 0) -- bounded tail\n")
}

# Return levels (negate back to recover actual temperatures in C)
rl_50  <- -qevd(1 - 1/50,  loc = mu, scale = sigma, shape = xi, type = "GEV")
rl_100 <- -qevd(1 - 1/100, loc = mu, scale = sigma, shape = xi, type = "GEV")

cat(sprintf("\n1-in-50  year return level: %.2f C\n", rl_50))
cat(sprintf("1-in-100 year return level: %.2f C\n", rl_100))

# Return period of Feb 2021 event
feb2021_p  <- pevd(feb2021_neg, loc = mu, scale = sigma, shape = xi, type = "GEV")
feb2021_rp <- 1 / (1 - feb2021_p)
cat(sprintf("Estimated return period of Feb 2021 (%.1f C): %.1f years\n",
            feb2021_tmin, feb2021_rp))


########################################
# 4. EXPLORATORY PLOTS
########################################

par(mfrow = c(2, 2))

# --- Panel 1: Histogram of all Jan-Feb daily Tmin ---
hist(tmin_winter,
     col    = "gray",
     prob   = TRUE,
     breaks = 30,
     main   = "DFW Daily Jan-Feb Tmin (Full Dist)",
     xlab   = "Daily Minimum Temperature (C)",
     ylab   = "Probability Density")
lines(density(tmin_winter, na.rm = TRUE), col = "red", lwd = 3)
abline(v = feb2021_tmin, col = "blue", lwd = 2, lty = 2)
legend("topright",
       legend = c("Density", "Feb 2021 Min"),
       col    = c("red", "blue"),
       lty    = c("solid", "dashed"),
       lwd    = c(1, 1))

# --- Panel 2: QQ Normal --- is the full distribution normal? ---
prob_all   <- (1:length(tmin_winter)) / (length(tmin_winter) + 1)
quant_norm <- qnorm(prob_all,
                    mean = mean(tmin_winter, na.rm = TRUE),
                    sd   = sd(tmin_winter,   na.rm = TRUE))
plot(sort(tmin_winter), quant_norm,
     main = "QQ Normal -- Full Jan-Feb Tmin",
     xlab = "Observations (C)",
     ylab = "Theoretical Normal Quantiles (C)",
     pch  = 19, cex = 0.5)
abline(0, 1, col = "red", lwd = 2)

# --- Panel 3: Histogram of clean block minima ---
hist(block_min_clean,
     col    = "gray",
     prob   = TRUE,
     breaks = 15,
     main   = "Annual Jan-Feb Block Minima -- DFW",
     xlab   = "Annual Minimum Temperature (C)",
     ylab   = "Probability Density")
lines(density(block_min_clean), col = "red", lwd = 3)
abline(v = feb2021_tmin, col = "blue", lwd = 2, lty = 2)
legend("topright",
       legend = c("Density", "Feb 2021"),
       col    = c("red", "blue"),
       lty    = c("solid", "dashed"),
       lwd    = c(1, 1))

# --- Panel 4: QQ GEV --- does the GEV fit the block minima? ---
prob_bm   <- (1:length(block_max_neg)) / (length(block_max_neg) + 1)
quant_gev <- qevd(prob_bm, loc = mu, scale = sigma, shape = xi, type = "GEV")

plot(quant_gev, sort(block_max_neg),
     main = "QQ GEV Fit -- Block Minima (Negated)",
     xlab = "Theoretical GEV Quantiles",
     ylab = "Observed (Negated Block Minima)",
     pch  = 19)
abline(0, 1, col = "red", lwd = 2)


########################################
# 5. RISK CURVE WITH GEV FIT
########################################

par(mfrow = c(1, 1))

# Empirical return periods for cold extremes
# P(Tmin <= t) from empirical CDF; return period = 1 / P(Tmin <= t)
cdf_data           <- ecdf(block_min_clean)
p_data             <- cdf_data(sort(block_min_clean))
return_periods_emp <- 1 / p_data

# Temperature sequence for GEV curve (actual C)
t_seq     <- seq(min(block_min_clean) - 10, max(block_min_clean) + 5, by = 0.1)
t_seq_neg <- -t_seq   # negated for GEV CDF evaluation

# GEV exceedance probability: P(Tmin <= t) = 1 - F_GEV(-t)
p_cold <- 1 - pevd(t_seq_neg, loc = mu, scale = sigma, shape = xi, type = "GEV")
rp_gev <- 1 / p_cold

# Plot empirical data points
plot(sort(block_min_clean), return_periods_emp,
     pch  = 19,
     log  = "y",
     main = "DFW Jan-Feb Cold Extreme Risk Curve",
     xlab = "Annual Minimum Temperature (C)",
     ylab = "Return Period (years)",
     xlim = c(min(block_min_clean) - 10, max(block_min_clean) + 2),
     ylim = c(1, 500),
     col  = "black")

# GEV fit line
lines(t_seq, rp_gev, col = "red", lwd = 3)

# Mark Feb 2021 event
abline(v = feb2021_tmin, col = "darkorange", lwd = 2, lty = 2)
abline(h = feb2021_rp,   col = "darkorange", lwd = 2, lty = 2)

# Mark 1/50 and 1/100 year reference lines
abline(h = 50,  col = "darkgreen", lwd = 2, lty = 2)
abline(h = 100, col = "darkgreen", lwd = 2, lty = 2)

legend("topright",
       legend = c("Empirical Data", "GEV Fit (L-moments)","90% Bootstrap CI",
                  "Feb 2021", "50 / 100 yr"),
       col    = c("black", "red", "blue","darkorange", "darkgreen"),
       lty    = c(NA, "solid", "dashed","dashed", "dotted"),
       pch    = c(19, NA, NA, NA, NA),
       lwd    = c(NA, 3, 2, 2, 1))


########################################
# 6. BOOTSTRAP UNCERTAINTY
########################################

N <- 500   # bootstrap samples
set.seed(42)

mu_array    <- rep(NA, N)
sigma_array <- rep(NA, N)
xi_array    <- rep(NA, N)

for (i in 1:N) {
  # Resample negated clean block minima with replacement
  dummy_data <- sample(block_max_neg, size = length(block_max_neg), replace = TRUE)

  # Fit GEV using L-moments (same method as best estimate)
  dummy_fit <- tryCatch(
    fevd(dummy_data, type = "GEV", method = "Lmoments"),
    error = function(e) NULL
  )

  if (!is.null(dummy_fit)) {
    mu_array[i]    <- dummy_fit$results["location"]
    sigma_array[i] <- dummy_fit$results["scale"]
    xi_array[i]    <- dummy_fit$results["shape"]
  }
}

# Remove failed fits
valid_idx   <- !is.na(xi_array)
mu_array    <- mu_array[valid_idx]
sigma_array <- sigma_array[valid_idx]
xi_array    <- xi_array[valid_idx]
cat(sprintf("\nSuccessful bootstrap fits: %d / %d\n", sum(valid_idx), N))

# 90% confidence intervals on parameters
mu_q    <- quantile(mu_array,    probs = c(0.05, 0.95))
sigma_q <- quantile(sigma_array, probs = c(0.05, 0.95))
xi_q    <- quantile(xi_array,    probs = c(0.05, 0.95))

cat("\n--- Bootstrap 90% Confidence Intervals ---\n")
cat(sprintf("  mu    (location): [%.4f, %.4f]\n", mu_q[1],    mu_q[2]))
cat(sprintf("  sigma (scale):    [%.4f, %.4f]\n", sigma_q[1], sigma_q[2]))
cat(sprintf("  xi    (shape):    [%.4f, %.4f]\n", xi_q[1],    xi_q[2]))

# Add bootstrap uncertainty envelope to risk curve using
# 5th and 95th percentile parameter sets
p_ci_lo  <- 1 - pevd(t_seq_neg, loc = mu_q[1], scale = sigma_q[1],
                     shape = xi_q[1], type = "GEV")
p_ci_hi  <- 1 - pevd(t_seq_neg, loc = mu_q[2], scale = sigma_q[2],
                     shape = xi_q[2], type = "GEV")
rp_ci_lo <- 1 / p_ci_lo
rp_ci_hi <- 1 / p_ci_hi

lines(t_seq, rp_ci_lo, col = "blue", lwd = 2, lty = 2)
lines(t_seq, rp_ci_hi, col = "blue", lwd = 2, lty = 2)

# Update legend to include CI
legend("topright",
       legend = c("Empirical Data", "GEV Fit (L-moments)",
                  "90% Bootstrap CI", "Feb 2021", "50 / 100 yr"),
       col    = c("black", "red", "blue", "darkorange", "darkgreen"),
       lty    = c(NA, "solid", "dashed", "dashed", "dotted"),
       pch    = c(19, NA, NA, NA, NA),
       lwd    = c(NA, 3, 2, 2, 1))


########################################
# 7. BOOTSTRAP PARAMETER HISTOGRAMS
########################################

par(mfrow = c(1, 3))

hist(mu_array,
     col    = "lightblue",
     prob   = TRUE,
     breaks = 20,
     main   = "Bootstrap: Location (mu)",
     xlab   = "mu (negated scale)")
abline(v = mean(mu_array), col = "red", lwd = 3)

hist(sigma_array,
     col    = "lightgreen",
     prob   = TRUE,
     breaks = 20,
     main   = "Bootstrap: Scale (sigma)",
     xlab   = "sigma")
abline(v = mean(sigma_array), col = "red", lwd = 3)

hist(xi_array,
     col    = "lightyellow",
     prob   = TRUE,
     breaks = 20,
     main   = "Bootstrap: Shape (xi)",
     xlab   = "xi")
abline(v = mean(xi_array), col = "red", lwd = 3)


########################################
# 8. BOOTSTRAP RETURN LEVEL CONFIDENCE INTERVALS
########################################

# Compute return levels for each bootstrap sample then take quantiles
boot_rl_50 <- sapply(seq_along(xi_array), function(i) {
  -qevd(1 - 1/50,  loc = mu_array[i], scale = sigma_array[i],
        shape = xi_array[i], type = "GEV")
})

boot_rl_100 <- sapply(seq_along(xi_array), function(i) {
  -qevd(1 - 1/100, loc = mu_array[i], scale = sigma_array[i],
        shape = xi_array[i], type = "GEV")
})

rl50_ci  <- quantile(boot_rl_50,  probs = c(0.05, 0.95), na.rm = TRUE)
rl100_ci <- quantile(boot_rl_100, probs = c(0.05, 0.95), na.rm = TRUE)


########################################
# 9. SUMMARY TABLES
########################################

cat("\n====== RETURN LEVEL SUMMARY TABLE ======\n")
cat(sprintf("  1-in-50  year event:  %.2f C  [90%% CI: %.2f to %.2f C]\n",
            rl_50, rl50_ci[1], rl50_ci[2]))
cat(sprintf("  1-in-100 year event:  %.2f C  [90%% CI: %.2f to %.2f C]\n",
            rl_100, rl100_ci[1], rl100_ci[2]))
cat(sprintf("  Feb 2021 return period (best estimate): %.1f years\n", feb2021_rp))
cat("=========================================\n")

cat("\n====== GEV PARAMETER SUMMARY TABLE ======\n")
cat(sprintf("  %-14s  %-16s  %-24s\n",
            "Parameter", "Best Estimate", "90% CI (Bootstrap)"))
cat(sprintf("  %s\n", paste(rep("-", 58), collapse = "")))
cat(sprintf("  %-14s  %-16.4f  [%.4f, %.4f]\n", "mu (location)", mu,    mu_q[1],    mu_q[2]))
cat(sprintf("  %-14s  %-16.4f  [%.4f, %.4f]\n", "sigma (scale)", sigma, sigma_q[1], sigma_q[2]))
cat(sprintf("  %-14s  %-16.4f  [%.4f, %.4f]\n", "xi (shape)",    xi,    xi_q[1],    xi_q[2]))
cat("==========================================\n")
