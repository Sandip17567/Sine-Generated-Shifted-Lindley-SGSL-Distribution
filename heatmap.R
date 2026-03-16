# Load necessary libraries
library(ggplot2)
library(reshape2)

# Define the range for theta and mu
theta_vals <- seq(0.5, 5, by = 0.5)
mu_vals <- seq(0.5, 5, by = 0.5)

# Initialize matrices to store moments
mean_matrix <- matrix(0, nrow = length(theta_vals), ncol = length(mu_vals))
variance_matrix <- matrix(0, nrow = length(theta_vals), ncol = length(mu_vals))
skewness_matrix <- matrix(0, nrow = length(theta_vals), ncol = length(mu_vals))
kurtosis_matrix <- matrix(0, nrow = length(theta_vals), ncol = length(mu_vals))

# Placeholder functions for moments (replace with actual implementations)
moment_mean <- function(theta, mu) {
  # Placeholder calculation
  return(theta + mu)
}

moment_variance <- function(theta, mu) {
  # Placeholder calculation
  return(theta * mu)
}

moment_skewness <- function(theta, mu) {
  # Placeholder calculation
  return(theta / mu)
}

moment_kurtosis <- function(theta, mu) {
  # Placeholder calculation
  return(theta^2 + mu^2)
}

# Calculate moments for each combination of theta and mu
for (i in 1:length(theta_vals)) {
  for (j in 1:length(mu_vals)) {
    theta <- theta_vals[i]
    mu <- mu_vals[j]
    mean_matrix[i, j] <- moment_mean(theta, mu)
    variance_matrix[i, j] <- moment_variance(theta, mu)
    skewness_matrix[i, j] <- moment_skewness(theta, mu)
    kurtosis_matrix[i, j] <- moment_kurtosis(theta, mu)
  }
}

# Convert matrices to data frames for ggplot
mean_df <- melt(mean_matrix)
mean_df
variance_df <- melt(variance_matrix)
skewness_df <- melt(skewness_matrix)
kurtosis_df <- melt(kurtosis_matrix)

# Add column names
colnames(mean_df) <- colnames(variance_df) <- colnames(skewness_df) <- colnames(kurtosis_df) <- c("ThetaIndex", "MuIndex", "Value")
mean_df$Theta <- rep(theta_vals, each = length(mu_vals))
mean_df$Mu <- rep(mu_vals, length(theta_vals))
variance_df$Theta <- mean_df$Theta
variance_df$Mu <- mean_df$Mu
skewness_df$Theta <- mean_df$Theta
skewness_df$Mu <- mean_df$Mu
kurtosis_df$Theta <- mean_df$Theta
kurtosis_df$Mu <- mean_df$Mu

# Plot heatmaps
plot_heatmap <- function(data, title) {
  ggplot(data, aes(x = Theta, y = Mu, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(data$Value)) +
    labs(title = title, x = expression(theta), y = expression(mu)) +
    theme_minimal()
}

mean_plot <- plot_heatmap(mean_df, "Mean of SGSL Distribution")
variance_plot <- plot_heatmap(variance_df, "Variance of SGSL Distribution")
skewness_plot <- plot_heatmap(skewness_df, "Skewness of SGSL Distribution")
kurtosis_plot <- plot_heatmap(kurtosis_df, "Kurtosis of SGSL Distribution")

# Print plots
print(mean_plot)
print(variance_plot)
print(skewness_plot)
print(kurtosis_plot)
