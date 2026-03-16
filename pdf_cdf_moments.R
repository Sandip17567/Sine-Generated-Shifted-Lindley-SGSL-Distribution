# ---- Define Functions (with support restriction) ----
F_shifted_lindley <- function(x, theta, mu) {
  ifelse(x > mu, 1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu)), 0)
}

f_shifted_lindley <- function(x, theta, mu) {
  ifelse(x > mu, (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu)), 0)
}

G_SGSL <- function(x, theta, mu) {
  sin((pi / 2) * F_shifted_lindley(x, theta, mu))
}

g_SGSL<- function(x, theta, mu) {
  ifelse(x > mu, (pi / 2) * abs(cos((pi / 2) * F_shifted_lindley(x, theta, mu))) * f_shifted_lindley(x, theta, mu), 0)
}

# ---- Parameter and Data Setup ----
mu <- 1  # Fixed theta
theta_values <- c(0.5,2, 2.5)  # Varying mu
x_values <- seq(0, 10, by = 0.01)  # Range of x, starting from 0 to show full support
colors <- c("blue", "green", "red")  # Colors for different mu values

# ---- Plotting ----
par(mfrow = c(1,2))
# PDF Plot (all mu values on one plot)
plot(x_values, g_SGSL(x_values, mu, theta_values[1]), type = "l", col = colors[1], lwd = 2,
     xlab = "x", ylab = "Density", main = "SGSL PDF (mu = 1)")
for (i in 2:length(theta_values)) {
  lines(x_values, g_SGSL(x_values, mu, theta_values[i]), col = colors[i], lwd = 2)
}
legend("topright", legend = paste("theta=", theta_values), col = colors, lwd = 2)
grid()
  # Reset to single plot layout
# CDF Plot (all mu values on one plot)
plot(x_values, G_SGSL(x_values, mu, theta_values[1]), type = "l", col = colors[1], lwd = 1,
     xlab = "x", ylab = "Cumulative Probability", main = "SGSL CDF (mu = 1)")
for (i in 2:length(theta_values)) {
  lines(x_values, G_SGSL(x_values, mu, theta_values[i]), col = colors[i], lwd = 1,add=T)
}
legend("bottomright", legend = paste("mu =", mu_values), col = colors, lwd = 2,)
grid()
# Moments Calculation & Plot (mean and variance)
mean_values <- sapply(mu_values, function(mu) integrate(function(x) x * g_SGSL(x, theta, mu), mu, Inf)$value)
var_values <- sapply(mu_values, function(mu) {
  second_moment <- integrate(function(x) x^2 * g_SGSL(x, theta, mu), mu, Inf)$value
  second_moment - mean_values[match(mu, mu_values)]^2
})

barplot(rbind(mean_values, var_values), beside = TRUE, col = c("orange", "purple"),
        names.arg = mu_values, legend.text = c("Mean", "Variance"), args.legend = list(x = "topleft"),
        xlab = "mu", ylab = "Value", main = "SGSL Moments (theta = 2)")
grid()

#generate sample observation from lindley and shifted lindley distribution and plot 