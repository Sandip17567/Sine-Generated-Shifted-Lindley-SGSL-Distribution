# ---- Define Functions ----

F_shifted_lindley <- function(x, theta, mu) {
  1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

f_shifted_lindley <- function(x, theta, mu) {
  (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

G_SGSL <- function(x, theta, mu) {
  sin((pi / 2) * F_shifted_lindley(x, theta, mu))
}

g_SGSL <- function(x, theta, mu) {
  (pi / 2) * cos((pi / 2) * F_shifted_lindley(x, theta, mu)) * f_shifted_lindley(x, theta, mu)
}


# ---- Parameter and Data Setup ----

theta <- 2
mu <- 1
x_values <- seq(mu, 10, by = 0.01)  # Finer resolution for smoother plots


# ---- Calculate & Plot PDF, CDF, and Moments ----

# Calculate
pdf_values <- g_SGSL(x_values, theta, mu)
cdf_values <- G_SGSL(x_values, theta, mu)

# Numerical Integration for Mean & Variance
mean_SGSL <- integrate(function(x) x * g_SGSL(x, theta, mu), mu, Inf)$value
second_moment_SGSL <- integrate(function(x) x^2 * g_SGSL(x, theta, mu), mu, Inf)$value
var_SGSL <- second_moment_SGSL - mean_SGSL^2

# Plotting
par(mfrow = c(1, 3))  # Arrange plots in a row

# PDF Plot
plot(x_values, pdf_values, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Density", main = "SGSL PDF")
abline(v = mean_SGSL, col = "red", lty = 2, lwd = 1.5)  # Add mean line
grid()

# CDF Plot
plot(x_values, cdf_values, type = "l", col = "green", lwd = 2,
     xlab = "x", ylab = "Cumulative Probability", main = "SGSL CDF")
grid()

# Moments Plot (Mean and Variance)
barplot(c(mean_SGSL, var_SGSL), names.arg = c("Mean", "Variance"),
        col = c("orange", "purple"), main = "SGSL Moments")
###############################################################


# ---- Define Functions ----

F_shifted_lindley <- function(x, theta, mu) {
  1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

f_shifted_lindley <- function(x, theta, mu) {
  (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

G_SGSL <- function(x, theta, mu) {
  sin((pi / 2) * F_shifted_lindley(x, theta, mu))
}

g_SGSL <- function(x, theta, mu) {
  (pi / 2) * cos((pi / 2) * F_shifted_lindley(x, theta, mu)) * f_shifted_lindley(x, theta, mu)
}


# ---- Parameter and Data Setup ----

mu <- 1  # Fixed mu
theta_values <- c(0.5, 2, 5)  # Varying theta
x_values <- seq(mu, 10, by = 0.01)  # Fine resolution for smoother plots

# Colors for different theta values
colors <- c("blue", "green", "red")


# ---- Plotting ----

par(mfrow = c(1, 2))  # Arrange plots in a row

# PDF Plot
plot(x_values, g_SGSL(x_values, theta_values[1], mu), type = "l", col = colors[1], lwd = 2,
     xlab = "x", ylab = "Density", main = "SGSL PDF for Different Theta Values")
for (i in 2:length(theta_values)) {
  lines(x_values, g_SGSL(x_values, theta_values[i], mu), col = colors[i], lwd = 2)
}
legend("topright", legend = paste("theta =", theta_values), col = colors, lwd = 2)

# CDF Plot
plot(x_values, G_SGSL(x_values, theta_values[1], mu), type = "l", col = colors[1], lwd = 2,
     xlab = "x", ylab = "Cumulative Probability", main = "SGSL CDF for Different Theta Values")
for (i in 2:length(theta_values)) {
  lines(x_values, G_SGSL(x_values, theta_values[i], mu), col = colors[i], lwd = 2)
}
legend("bottomright", legend = paste("theta =", theta_values), col = colors, lwd = 2)

#################################


# ---- Define Functions ----

# ... (F_shifted_lindley, f_shifted_lindley, G_SGSL, g_SGSL remain the same)


# ---- Parameter and Data Setup ----

theta_values <- c(2, 5)  # Two fixed theta values
mu_values <- c(0.5, 1, 1.5)  # Varying mu values for each theta
x_values <- seq(0, 10, by = 0.05)  # Range of x values (wider for better visualization)

# Colors for different mu values (same for both theta plots)
colors <- c("blue", "green", "red")


# ---- Plotting ----

par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid

for (theta in theta_values) {
  # PDF Plot for the current theta
  plot(x_values, g_SGSL(x_values, theta, mu_values[1]), type = "l", col = colors[1], lwd = 1,
       xlab = "x", ylab = "Density", main = paste("SGSL PDF (theta =", theta, ")"))
  for (i in 2:length(mu_values)) {
    lines(x_values, g_SGSL(x_values, theta, mu_values[i]), col = colors[i], lwd = 1)
  }
  legend("topright", legend = paste("mu =", mu_values), col = colors, lwd = 1)
  
  # CDF Plot for the current theta
  plot(x_values, G_SGSL(x_values, theta, mu_values[1]), type = "l", col = colors[1], lwd = 1,
       xlab = "x", ylab = "Cumulative Probability", main = paste("SGSL CDF (theta =", theta, ")"))
  for (i in 2:length(mu_values)) {
    lines(x_values, G_SGSL(x_values, theta, mu_values[i]), col = colors[i], lwd = 1)
  }
  legend("bottomright", legend = paste("mu =", mu_values), col = colors, lwd = 1)
}

###################################


library(tidyverse)
library(moments)  # For calculating skewness and kurtosis

# ---- Functions ----

# ... (F_shifted_lindley, f_shifted_lindley, G_SGSL, g_SGSL remain the same)

# Function to calculate moments of SGSL
calc_moments <- function(theta, mu) {
  fn_mean <- function(x) x * g_SGSL(x, theta, mu)
  fn_2nd_moment <- function(x) x^2 * g_SGSL(x, theta, mu)
  
  mean <- integrate(fn_mean, mu, Inf)$value
  second_moment <- integrate(fn_2nd_moment, mu, Inf)$value
  variance <- second_moment - mean^2
  
  # Sample from the SGSL distribution to calculate skewness and kurtosis
  samples <- rSGSL(100000, theta, mu) 
  skewness <- skewness(samples)
  kurtosis <- kurtosis(samples)
  
  return(data.frame(mean = mean, var = variance, skewness = skewness, kurtosis = kurtosis))
}

# ---- Create a grid of parameter values ----

theta_values <- seq(0.1, 5, length.out = 20)
mu_values <- seq(0.1, 5, length.out = 20)
params_grid <- expand.grid(theta = theta_values, mu = mu_values)

# ---- Calculate moments for each parameter combination ----

moment_df <- params_grid %>%
  rowwise() %>%
  mutate(calc_moments(theta, mu)) %>%
  ungroup() %>%
  pivot_longer(-c(theta, mu), names_to = "moment", values_to = "value")

# ---- Create Heatmaps ----

ggplot(moment_df, aes(mu, theta, fill = value)) +
  geom_tile() +
  facet_wrap(~ moment, ncol = 2, scales = "free") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(moment_df$value)) +
  labs(x = "mu", y = "theta", fill = "Value") +
  ggtitle("Heatmap of SGSL Moments (Mean, Variance, Skewness, Kurtosis)")

