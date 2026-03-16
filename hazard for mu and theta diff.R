# Define the shifted Lindley CDF
F_shifted_lindley <- function(x, theta, mu) {
  1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

# Define the shifted Lindley PDF
f_shifted_lindley <- function(x, theta, mu) {
  (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

# Define the SGSL CDF
G_SGSL <- function(x, theta, mu) {
  sin((pi / 2) * F_shifted_lindley(x, theta, mu))
}

# Define the SGSL PDF
g_SGSL <- function(x, theta, mu) {
  (pi / 2) * cos((pi / 2) * F_shifted_lindley(x, theta, mu)) * f_shifted_lindley(x, theta, mu)
}

# Define the SGSL hazard function
h_SGSL <- function(x, theta, mu) {
  g_SGSL(x, theta, mu) / (1 - G_SGSL(x, theta, mu))
}

# Example parameters
theta <- 2
mu <- 1

# Example x values
x_values <- seq(mu, mu + 5, by = 0.1)

# Calculate the hazard function values
hazard_values <- h_SGSL(x_values, theta, mu)

# Plot the hazard function
plot(x_values, hazard_values, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Hazard function h_SGSL(x; theta, mu)",
     main = "Hazard Function of the SGSL Distribution")

#################################

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

h_SGSL <- function(x, theta, mu) {
  g_SGSL(x, theta, mu) / (1 - G_SGSL(x, theta, mu))
}


# ---- Parameter Exploration & Plots ----

# Parameters
theta_values <- c(0.5, 1, 2)  # Varying theta
mu_values <- c(0.5, 1, 2)   # Varying mu
x_values <- seq(0, 10, by = 0.1)  # Wider range of x values

# Plot Setup
par(mfrow = c(length(theta_values), length(mu_values)))

for (theta in theta_values) {
  for (mu in mu_values) {
    
    # Calculate Hazard
    hazard_values <- h_SGSL(x_values, theta, mu)
    
    # Enhanced Plot
    plot(x_values, hazard_values, type = "l", col = "blue", lwd = 2,
         xlab = "x", ylab = "Hazard function h_SGSL(x)",
         main = paste("SGSL Hazard: theta =", theta, ", mu =", mu),
         ylim = c(0, max(hazard_values) * 1.2))  # Adjust y-axis for clarity
    grid(lty = 1, col = "lightgray")  # Add gridlines
  }
}


# ---- Statistical Interpretation of Hazard Function ----

# Hazard function's interpretation:
# - Instantaneous risk of an event (e.g., failure, death) occurring at time x,
#   given that the individual or item has survived up to time x.
# - High hazard indicates a greater risk of the event happening soon.
# - Low hazard suggests a lower risk of the event in the near future.

# Examples of interpretations based on your plots:
# - If the hazard is increasing, the risk of the event increases over time.
# - If the hazard is decreasing, the risk of the event decreases over time.
# - A constant hazard indicates a constant risk of the event throughout.