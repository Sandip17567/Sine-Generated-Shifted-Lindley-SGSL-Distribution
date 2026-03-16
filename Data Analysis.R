# SGSL Complete Toolbox
# Purpose: Random number generation, density/CDF, MGF, diagnostics, and model fitting for SGSL
# Date: May 20, 2025

## 0 ──────────────────────────────────────────────────────────────────────
# Check and install required packages
suns
pkgs<-c("fitdistrplus", "ggplot2", "stats4", "dplyr", "patchwork", "knitr", "tibble")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load dataset
data <- c(0.4, 0.6, 0.7, 0.9, 1.0, 1.2, 1.4, 1.6, 1.9, 2.2, 2.5, 2.9, 3.3, 3.8, 4.5, 
          5.2, 6.0, 7.0, 8.2, 9.5, 11.0, 13.0, 15.5, 18.0, 21.0, 25.0, 30.0, 36.0, 
          43.0, 51.0, 60.0, 70.0, 82.0, 95.0, 110.0, 130.0, 150.0, 175.0, 200.0, 
          230.0, 270.0, 320.0, 380.0, 450.0, 530.0, 620.0, 720.0, 850.0, 1000.0, 1200.0)
samples<-data/1000
## 1  Shifted-Lindley building blocks ─────────────────────────────────────
F_shifted <- function(x, theta, mu) {
  if (theta <= 0 || mu < 0) stop("theta must be positive, mu non-negative")
  ifelse(x > mu,
         1 - (1 + theta * (1 + x)) / (1 + theta * (1 + mu)) * exp(-theta * (x - mu)),
         0)
}

f_shifted <- function(x, theta, mu) {
  if (theta <= 0 || mu < 0) stop("theta must be positive, mu non-negative")
  ifelse(x > mu,
         theta^2 / (1 + theta * (1 + mu)) * exp(-theta * (x - mu)),
         0)
}

## 2  SGSL CDF / PDF -------------------------------------------------------
G_SGSL <- function(x, theta, mu) {
  sin(pi / 2 * F_shifted(x, theta, mu))
}

g_SGSL <- function(x, theta, mu) {
  cos_term <- cos(pi / 2 * F_shifted(x, theta, mu))
  dens <- ifelse(x > mu, (pi / 2) * cos_term * f_shifted(x, theta, mu), 0)
  pmax(dens, 0)  # Ensure non-negative density
}

## 3  Quantile of shifted-Lindley (numeric inversion) ----------------------
q_shifted <- function(p, theta, mu) {
  if (p < 0 || p > 1) stop("p must be between 0 and 1")
  uniroot(function(x) F_shifted(x, theta, mu) - p,
          interval = c(mu, mu + 1000), tol = 1e-10, extendInt = "yes")$root
}

## 4  SGSL random generator ------------------------------------------------
rSGSL <- function(n, theta, mu) {
  if (n <= 0) stop("n must be positive")
  u <- runif(n)
  p <- (2 / pi) * asin(u)  # Map U~Unif → shifted-Lindley prob.
  vapply(p, q_shifted, numeric(1), theta = theta, mu = mu)
}

## 5  Numerical MGF --------------------------------------------------------
MGF_SGSL <- function(t, theta, mu) {
  sapply(t, function(tt)
    integrate(function(x) exp(tt * x) * g_SGSL(x, theta, mu),
              lower = mu, upper = mu + 100, rel.tol = 1e-6)$value)
}

## 6  Densities for fitdistrplus  (names must start with d…) -------------
dSGSL <- function(x, theta, mu, log = FALSE) {
  z <- g_SGSL(x, theta, mu)
  z <- pmax(z, 1e-10)  # Prevent log(0)
  if (log) log(z) else z
}

dShL <- function(x, theta, mu, log = FALSE) {
  z <- f_shifted(x, theta, mu)
  z <- pmax(z, 1e-10)  # Prevent log(0)
  if (log) log(z) else z
}

## 7  Fitting helper -------------------------------------------------------
fit_all <- function(data, start = c(theta = 1 / (mean(data) - min(data)), mu = min(data) - 0.05)) {
  list(
    SGSL = tryCatch({
      fit <- fitdist(data, "SGSL", start = as.list(start), lower = c(1e-3, 0),
                     method = "mle", optim.method = "Nelder-Mead",
                     control = list(maxit = 2000, trace = 1))
      message("SGSL fit: theta = ", round(fit$estimate[1], 4), ", mu = ", round(fit$estimate[2], 4),
              ", loglik = ", round(fit$loglik, 2))
      fit
    }, error = function(e) {
      message("SGSL fit failed: ", e$message)
      NULL
    }),
    ShL = tryCatch({
      fit <- fitdist(data, "ShL", start = as.list(start), lower = c(1e-3, 0),
                     method = "mle", optim.method = "Nelder-Mead",
                     control = list(maxit = 2000))
      message("ShL fit: theta = ", round(fit$estimate[1], 4), ", mu = ", round(fit$estimate[2], 4))
      fit
    }, error = function(e) {
      message("ShL fit failed: ", e$message)
      NULL
    }),
    Lindley = tryCatch({
      fit <- fitdist(data, "gamma", start = list(shape = 1, scale = mean(data)))
      message("Lindley fit: shape = ", round(fit$estimate[1], 4), ", scale = ", round(fit$estimate[2], 4))
      fit
    }, error = function(e) {
      message("Lindley fit failed: ", e$message)
      NULL
    }),
    Exp = tryCatch({
      fit <- fitdist(data, "exp")
      message("Exp fit: rate = ", round(fit$estimate[1], 4))
      fit
    }, error = function(e) {
      message("Exp fit failed: ", e$message)
      NULL
    })
  )
}

## 8  Diagnostic plot (hist+PDF, ECDF+CDF) ---------------------------------
check_plots <- function(theta, mu, n = 1e4, file = "SGSL_check.pdf") {
  samp <- rSGSL(n, theta, mu)
  xgrid <- seq(mu, quantile(samp, 0.999, na.rm = TRUE), length = 400)
  pdfdf <- tibble(x = xgrid, y = g_SGSL(xgrid, theta, mu))
  cdfdf <- tibble(x = xgrid, y = G_SGSL(xgrid, theta, mu))
  
  g1 <- ggplot() +
    geom_histogram(aes(samp, ..density..), bins = 40,
                   fill = "grey80", colour = "grey30") +
    geom_line(data = pdfdf, aes(x, y),
              colour = "firebrick", linewidth = 1) +
    labs(title = bquote("SGSL PDF (theta = " * .(theta) * ", mu = " * .(mu) * ")"),
         x = NULL, y = "Density") +
    theme_minimal(14)
  
  g2 <- ggplot() +
    stat_ecdf(aes(samp), colour = "steelblue", linewidth = 0.9) +
    geom_line(data = cdfdf, aes(x, y),
              colour = "firebrick", linewidth = 1) +
    labs(title = "ECDF vs Theoretical CDF",
         x = NULL, y = "Probability") +
    theme_minimal(14)
  
  tryCatch(
    ggsave(file, g1 + g2, width = 8, height = 4.5, device = cairo_pdf),
    error = function(e) message("Failed to save plot: ", e$message)
  )
  message("Diagnostic figure saved as: ", file)
}

## 9  Overlay of fitted density ------------------------------------------
plot_fits <- function(data, fits, theta0, mu0, file = "SGSL_fits.pdf") {
  xmax <- max(data, na.rm = TRUE) * 1.05
  xg <- seq(mu0, xmax, length = 600)
  dens <- tibble(
    x = rep(xg, 4),
    density = c(
      if (!is.null(fits$SGSL)) dSGSL(xg, fits$SGSL$estimate[1], fits$SGSL$estimate[2]) else rep(NA, length(xg)),
      if (!is.null(fits$ShL)) dShL(xg, fits$ShL$estimate[1], fits$ShL$estimate[2]) else rep(NA, length(xg)),
      if (!is.null(fits$Lindley)) dgamma(xg, shape = fits$Lindley$estimate[1], scale = fits$Lindley$estimate[2]) else rep(NA, length(xg)),
      if (!is.null(fits$Exp)) dexp(xg, rate = 1 / fits$Exp$estimate) else rep(NA, length(xg))
    ),
    model = factor(rep(c("SGSL", "Shifted-Lindley", "Lindley≈Gamma", "Exponential"),
                       each = length(xg)),
                   levels = c("SGSL", "Shifted-Lindley", "Lindley≈Gamma", "Exponential"))
  )
  
  p <- ggplot() +
    geom_histogram(aes(data, ..density..), bins = 40,
                   fill = "grey85", colour = "grey25") +
    geom_line(data = dens, aes(x, density, colour = model),
              linewidth = 1, na.rm = TRUE) +
    scale_colour_manual(values = c("midnightblue", "purple", "tomato", "darkgreen")) +
    labs(title = bquote("Fitted Density (theta = " * .(theta0) * ", mu = " * .(mu0) * ")"),
         colour = NULL, x = NULL, y = "Density") +
    theme_minimal(15)
  
  tryCatch(
    ggsave(file, p, width = 6.5, height = 4.5, device = cairo_pdf),
    error = function(e) message("Failed to save plot: ", e$message)
  )
  message("Overlay figure saved as: ", file)
}

## 10  FULL DEMO – reproducible --------------------------------------------
set.seed(2025)
theta0 <- 1.5
mu0 <- 0.5
n <- 1000
#samples<-rSGSL(n,2.24,2.15)
# 10a  Diagnostic PDF/CDF plot
check_plots(theta0, mu0, n = 8000)  # → SGSL_check.pdf

# 10b  Fit four models
fits <- fit_all(samples, start = c(theta = 1 / (mean(samples) - min(samples)), mu = 0.5))
if (!is.null(fits$SGSL)) {
  message("SGSL fit successful: theta = ", round(fits$SGSL$estimate[1], 4),
          ", mu = ", round(fits$SGSL$estimate[2], 4),
          ", AIC = ", round(AIC(fits$SGSL), 1),
          ", BIC = ", round(-2 * fits$SGSL$loglik + length(fits$SGSL$estimate) * log(length(samples)), 1))
}

# 10c  AIC table (and save LaTeX)
AICtab <- data.frame(
  Model = c("SGSL", "Shifted-Lindley", "Lindley≈Gamma", "Exponential"),
  AIC = c(
    if (!is.null(fits$SGSL)) AIC(fits$SGSL) else NA,
    if (!is.null(fits$ShL)) AIC(fits$ShL) else NA,
    if (!is.null(fits$Lindley)) fits$Lindley$aic else NA,
    if (!is.null(fits$Exp)) fits$Exp$aic else NA
  )
)
tryCatch(
  writeLines(
    knitr::kable(AICtab, format = "latex", booktabs = TRUE,
                 caption = "AIC comparison on provided data",
                 digits = 1),
    "AIC_table.tex"
  ),
  error = function(e) message("Failed to save AIC table: ", e$message)
)
message("LaTeX AIC table written to: AIC_table.tex")

# 10d  Overlay density plot
plot_fits(samples, fits, theta0, mu0)  # → SGSL_fits.pdf

# 10c.1  BIC table
logL <- function(fit) if (!is.null(fit)) fit$loglik else NA
n <- length(samples)

k_SGSL <- if (!is.null(fits$SGSL)) length(fits$SGSL$estimate) else 0
k_ShL <- if (!is.null(fits$ShL)) length(fits$ShL$estimate) else 0
k_Lindley <- if (!is.null(fits$Lindley)) length(fits$Lindley$estimate) else 0
k_Exp <- if (!is.null(fits$Exp)) length(fits$Exp$estimate) else 0

BICtab <- data.frame(
  Model = c("SGSL", "Shifted-Lindley", "Lindley≈Gamma", "Exponential"),
  BIC = c(
    if (!is.na(logL(fits$SGSL))) -2 * logL(fits$SGSL) + k_SGSL * log(n) else NA,
    if (!is.na(logL(fits$ShL))) -2 * logL(fits$ShL) + k_ShL * log(n) else NA,
    if (!is.na(logL(fits$Lindley))) -2 * logL(fits$Lindley) + k_Lindley * log(n) else NA,
    if (!is.na(logL(fits$Exp))) -2 * logL(fits$Exp) + k_Exp * log(n) else NA
  )
)

tryCatch(
  writeLines(
    knitr::kable(BICtab, format = "latex", booktabs = TRUE,
                 caption = "BIC comparison on provided data",
                 digits = 1),
    "BIC_table.tex"
  ),
  error = function(e) message("Failed to save BIC table: ", e$message)
)
message("LaTeX BIC table written to: BIC_table.tex")

