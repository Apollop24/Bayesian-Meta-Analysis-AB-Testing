set.seed(123)  # Set seed for reproducibility

# ==============================================================================
# AUTOMATICALLY INSTALL AND LOAD REQUIRED LIBRARIES
# ==============================================================================

install_and_load_packages <- function() {
  required_libraries <- c("runjags", "coda", "dplyr", "tidyr", "ggplot2", "parallel", "future", "future.apply")
  
  # Identify packages that are not installed
  new_packages <- required_libraries[!sapply(required_libraries, requireNamespace, quietly = TRUE)]
  
  # Install missing packages
  if (length(new_packages) > 0) {
    install.packages(new_packages, dependencies = TRUE)
  }
  
  # Load all required libraries
  invisible(sapply(required_libraries, function(lib) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }))
}

# Run the package installation and loading function
install_and_load_packages()

# Set up future plan for parallel processing (multisession for Windows)
future::plan(multisession, workers = parallel::detectCores())

# ==============================================================================
# DEFINE GROUPED DATA FOR ANALYSIS AS DATA FRAMES
# ==============================================================================

# Create a new environment to store analysis datasets
analysis_env <- new.env()
with(analysis_env, {
  
  
  mattress_low_conversion_1 <- data.frame(
    revenue_control = c(0, 1406, 1105, 0, 0, 0, 1347, 0, 648, 0),
    revenue_variation = c(405, 0, 897, 1980, 357, 1650, 0, 134, 0, 974),
    conversions_control = c(0, 1, 1, 0, 0, 0, 3, 0, 1, 0),
    conversions_variation = c(2, 0, 2, 2, 1, 1, 0, 1, 0, 1),
    recipients_control = c(4091, 5758, 14231, 5142, 6774, 4594, 10147, 3589, 5472, 13072),
    recipients_variation = c(4074, 5749, 14192, 5121, 6790, 4590, 10154, 3586, 5440, 13055)
  )
  
  
  
 
  
})

# ==============================================================================
# FUNCTION: PREPROCESS DATA FOR ANALYSIS
# Purpose: Handle missing values, zeros and perform validation.
# ==============================================================================

preprocess_data <- function(data) {
  library(dplyr)
  
  # Convert data to data frame format
  processed <- as.data.frame(data)
  
  # Identify numeric columns
  numeric_cols <- names(processed)[sapply(processed, is.numeric)]
  
  # Initialize log messages for tracking preprocessing changes
  log_messages <- c()
  
  # ============================================================================
  # FUNCTION: VALIDATE DATASET VALUES
  # Purpose: Ensure no invalid values in the data
  # ============================================================================
  validation_checks <- function(data) {
    if (any(data$revenue_control < 0) || any(data$revenue_variation < 0)) {
      stop("Invalid negative revenue detected")
    }
    if (any(data$recipients_control <= 0) || any(data$recipients_variation <= 0)) {
      stop("Invalid recipient counts - must be positive")
    }
    for (group in c("control", "variation")) {
      conv_rate <- data[[paste0("conversions_", group)]] / 
        data[[paste0("recipients_", group)]]
      if (any(conv_rate > 1, na.rm = TRUE)) {
        stop(paste("Invalid conversion rate >100% detected in", group))
      }
    }
  }
  
  # Run validation checks with error handling
  tryCatch({
    validation_checks(processed)
  }, error = function(e) {
    stop(paste("Validation failed:", e$message))
  })
  
  # ============================================================================
  # FUNCTION: HANDLE ZERO VALUES
  # Purpose: Handle zeros to avoid computational issues in Bayesian modeling.
  # ============================================================================
  handle_zeros <- function(data) {
    min_nonzero_revenue <- min(c(data$revenue_control[data$revenue_control > 0], 
                                 data$revenue_variation[data$revenue_variation > 0]), na.rm = TRUE)
    small_revenue <- min_nonzero_revenue * 0.001
    
    if (any(data$revenue_control == 0)) {
      data$revenue_control[data$revenue_control == 0] <- small_revenue
      log_messages <<- c(log_messages, paste("Replaced zeros in control revenue with", small_revenue))
    }
    
    if (any(data$revenue_variation == 0)) {
      data$revenue_variation[data$revenue_variation == 0] <- small_revenue
      log_messages <<- c(log_messages, paste("Replaced zeros in variation revenue with", small_revenue))
    }
    
    data$conversions_control <- round(data$conversions_control)
    data$conversions_variation <- round(data$conversions_variation)
    
    return(data)
  }
  
  # Run zero handling function
  processed <- handle_zeros(processed)
  
  # ============================================================================
  # FUNCTION: HANDLE MISSING VALUES
  # Purpose: Fill missing values in revenue and conversion columns
  # ============================================================================
  handle_missing_values <- function(data, col) {
    if (any(is.na(data[[col]]))) {
      na_count <- sum(is.na(data[[col]]))
      
      if (grepl("revenue|conversion", col, ignore.case = TRUE)) {
        data[[col]][is.na(data[[col]])] <- 0
        log_messages <<- c(log_messages, paste("Replaced", na_count, "missing values with 0 in", col))
      }
    }
    return(data)
  }
  
  # Apply missing value handling to all numeric columns
  for (col in numeric_cols) {
    processed <- handle_missing_values(processed, col)
  }
  
  # Store preprocessing logs and data quality attributes
  attr(processed, "preprocessing_log") <- log_messages
  
  return(processed)
}

# ==============================================================================
# JAGS MODEL FOR REVENUE PER RECIPIENT (RPR)
# ==============================================================================
model_rpr <- "
model {
  # Priors for global parameters
  mu_global ~ dnorm(0, 0.001)    # Prior for global mean (log scale)
  
  # Use a half-Cauchy prior for the between-study standard deviation.
  # Here, sigma_between is given a truncated t-distribution (df=1, scale=10)
  # which is equivalent to a half-Cauchy. Do not supply initial values for tau_between.
  sigma_between ~ dt(0, 0.01, 1) T(0,)
  tau_between <- 1/(sigma_between^2)  # Between-study precision computed from sigma_between
  
  psi ~ dbeta(1, 1)  # Prior for probability of binary indicator (Bernoulli)
  
  for (i in 1:n) {
    # Binary latent indicators for control and variation groups
    z_control[i] ~ dbern(psi)
    z_variation[i] ~ dbern(psi)

    # Study-level mean effect (theta)
    theta[i] ~ dnorm(mu_global, tau_between)

    # Revenue distributions (Gamma distributed)
    revenue_control[i] ~ dgamma(shape_c[i], rate_c[i])
    revenue_variation[i] ~ dgamma(shape_v[i], rate_v[i])

    # Shape parameter transformations for gamma distribution
    log(shape_c[i]) <- alpha_c[i]
    log(shape_v[i]) <- alpha_v[i]

    # Priors for alpha parameters (log-scale shape parameters)
    alpha_c[i] ~ dnorm(0, 0.001)
    alpha_v[i] ~ dnorm(0, 0.001)

    # Compute rate parameters based on mean revenue per recipient
    rate_c[i] <- shape_c[i] / (mu_c[i] * recipients_control[i])
    rate_v[i] <- shape_v[i] / (mu_v[i] * recipients_variation[i])

    # Log-transformed mean revenue per recipient
    log(mu_c[i]) <- theta[i]
    log(mu_v[i]) <- theta[i] + delta[i]

    # Study-level treatment effect (delta)
    delta[i] ~ dnorm(0, tau_between)

    # Compute Revenue Per Recipient (RPR)
    rpr_control[i] <- revenue_control[i] / recipients_control[i]
    rpr_variation[i] <- revenue_variation[i] / recipients_variation[i]

    # Compute lift in RPR percentage
    lift_rpr[i] <- ((rpr_variation[i] + epsilon) - (rpr_control[i] + epsilon)) / 
                   (rpr_control[i] + epsilon) * 100

    # Compute standardized effect size
    effect_size[i] <- abs(delta[i]) / sqrt(1/shape_c[i] + 1/shape_v[i])

    # Compute precision
    precision[i] <- 1 / sqrt(1/shape_c[i] + 1/shape_v[i])

    # Compute probability of a positive effect
    prob_positive[i] <- mean(delta[i] > 0)
  }

  # Small constant to prevent division errors
  epsilon <- 0.0001

  # Compute global lift estimate (percentage)
  global_lift_rpr <- (exp(mu_global + mean(delta[1:n])) - exp(mu_global)) / 
                     exp(mu_global) * 100

  # Compute probability of a positive effect
  prob_positive_rpr <- mean(delta[1:n] > 0)

  # Compute heterogeneity measures
  I2 <- tau_between / (tau_between + mean(precision[1:n]))  # Proportion of variance due to heterogeneity
  H2 <- 1 / (1 - I2)  # Heterogeneity ratio

  # Compute predictive intervals
  pred_sigma <- sqrt(1/tau_between + mean(1/precision[1:n]))
  pi_lower <- mu_global - 1.96 * pred_sigma
  pi_upper <- mu_global + 1.96 * pred_sigma
}
"



# ==============================================================================
# JAGS MODEL FOR CONVERSION RATE (CR)
# ==============================================================================
# ==============================================================================
# JAGS MODEL FOR CONVERSION RATE (CR)
# ==============================================================================
model_cr <- "
model {
  # Priors for global parameters
  mu_global ~ dnorm(0, 0.001)   # Prior for global mean effect
  
  # Use a half-Cauchy prior for the between-study standard deviation.
  sigma_between ~ dt(0, 0.01, 1) T(0,)
  tau_between <- 1/(sigma_between^2)  # Derived between-study precision

  for (i in 1:n) {
    # Binomial likelihood for conversion counts
    conversions_control[i] ~ dbin(p_control[i], recipients_control[i])
    conversions_variation[i] ~ dbin(p_variation[i], recipients_variation[i])

    # Logit-transformed conversion probabilities
    p_control[i] <- ilogit(mu_c[i])
    p_variation[i] <- ilogit(mu_v[i])

    # Study-level mean conversion rates
    mu_c[i] ~ dnorm(base[i], tau_between)
    mu_v[i] ~ dnorm(base[i] + delta[i], tau_between)

    # Baseline conversion rate
    base[i] ~ dnorm(mu_global, tau_between)

    # Study-level treatment effect (delta)
    delta[i] ~ dnorm(0, tau_between)

    # Compute lift in conversion rate percentage
    lift_cr[i] <- ((p_variation[i] + epsilon) - (p_control[i] + epsilon)) / 
                  (p_control[i] + epsilon) * 100

    # Compute standardized effect size
    effect_size[i] <- abs(delta[i]) * sqrt(recipients_control[i] * p_control[i] * (1 - p_control[i]))

    # Compute precision
    precision[i] <- sqrt(recipients_control[i] * p_control[i] * (1 - p_control[i]))

    # Compute probability of a positive effect
    prob_positive[i] <- step(delta[i])
  }

  # Small constant to prevent division errors
  epsilon <- 0.0001

  # Compute global lift estimate (percentage)
  global_lift_cr <- (ilogit(mu_global + mean(delta[1:n])) - ilogit(mu_global)) / 
                    ilogit(mu_global) * 100

  # Compute probability of positive effect
  prob_positive_cr <- step(mean(delta[1:n]))

  # Compute heterogeneity measures
  I2 <- tau_between / (tau_between + mean(precision[1:n]))  # Proportion of variance due to heterogeneity
  H2 <- 1 / (1 - I2)  # Heterogeneity ratio

  # Compute predictive intervals
  pred_sigma <- sqrt(1/tau_between + mean(1/precision[1:n]))
  pi_lower <- mu_global - 1.96 * pred_sigma
  pi_upper <- mu_global + 1.96 * pred_sigma
}
"


# ==============================================================================
# FUNCTION: GENERATE INITIAL VALUES FOR MCMC CHAINS
# ==============================================================================
# ==============================================================================
# FUNCTION: GENERATE INITIAL VALUES FOR MCMC CHAINS
# ==============================================================================
generate_inits <- function(n_chains, data, monitor_params) {
  inits <- list()
  
  # Generate reasonable starting values for the model
  get_starting_values <- function() {
    # Compute mean conversion and revenue rates
    mean_conv_control <- mean(data$conversions_control / data$recipients_control, na.rm = TRUE)
    mean_rev_control <- mean(data$revenue_control / data$recipients_control, na.rm = TRUE)
    
    # Initialize key parameters with random but reasonable values.
    # Note: We no longer include tau_between as it is derived from sigma_between.
    starting_values <- list(
      mu_global = rnorm(1, log(mean_rev_control + 0.0001), 0.1),  # Log-scale mean revenue
      sigma_between = runif(1, 0.5, 2),  # Initial value for between-study SD (must be > 0)
      delta = rnorm(nrow(data), 0, 0.1)  # Study-level treatment effects
    )
    
    # If model includes baseline conversion rates, initialize them
    if ("base" %in% monitor_params) {
      starting_values$base <- rep(log(mean_conv_control / (1 - mean_conv_control + 0.0001)), nrow(data))
    }
    
    return(starting_values)
  }
  
  # Generate initial values for each MCMC chain
  for (i in 1:n_chains) {
    inits[[i]] <- c(
      get_starting_values(),
      list(
        .RNG.name = switch(i,
                           "base::Wichmann-Hill",
                           "base::Marsaglia-Multicarry",
                           "base::Super-Duper",
                           "base::Mersenne-Twister"),  # Different RNGs for each chain
        .RNG.seed = sample.int(1e6, 1)  # Random seed for reproducibility
      )
    )
  }
  
  return(inits)
}


# ==============================================================================
# FUNCTION: RUN BAYESIAN META-ANALYSIS
# Purpose: This function runs Bayesian meta-analysis using JAGS for two 
# key metrics: Revenue Per Recipient (RPR) and Conversion Rate (CR). 
# It preprocesses data, sets up JAGS input, runs models in parallel, 
# and returns results with diagnostics.
# ==============================================================================
run_bayesian_meta_analysis <- function(data, group_name) {
  # Preprocess the input dataset 
  processed_data <- preprocess_data(data)
  
  # Prepare data for RPR analysis (Revenue Per Recipient)
  jags_data_rpr <- list(
    revenue_control = processed_data$revenue_control,
    revenue_variation = processed_data$revenue_variation,
    recipients_control = processed_data$recipients_control,
    recipients_variation = processed_data$recipients_variation,
    n = nrow(processed_data)
  )
  
  # Prepare data for CR analysis (Conversion Rate)
  jags_data_cr <- list(
    conversions_control = processed_data$conversions_control,
    conversions_variation = processed_data$conversions_variation,
    recipients_control = processed_data$recipients_control,
    recipients_variation = processed_data$recipients_variation,
    n = nrow(processed_data)
  )
  
  # Define monitored parameters for RPR analysis
  monitor_params_rpr <- c(
    "mu_global", "tau_between", "sigma_between", "psi",
    "global_lift_rpr", "prob_positive_rpr",
    "lift_rpr", "rpr_control", "rpr_variation",
    "effect_size", "precision",
    "I2", "H2", "pi_lower", "pi_upper"
  )
  
  # Define monitored parameters for CR analysis
  monitor_params_cr <- c(
    "mu_global", "tau_between", "sigma_between",
    "global_lift_cr", "prob_positive_cr",
    "lift_cr", "p_control", "p_variation",
    "effect_size", "precision",
    "I2", "H2", "pi_lower", "pi_upper"
  )
  
  # Determine the number of chains based on available CPU cores (max 4)
  n_chains_used <- ifelse(parallel::detectCores() < 4, parallel::detectCores(), 4)
  
  # Generate initial values for the RPR model
  inits_rpr <- generate_inits(n_chains = n_chains_used, data = processed_data, monitor_params = monitor_params_rpr)
  
  # ----------------------------------------------------------------------------
  # RUN JAGS MODEL FOR RPR ANALYSIS
  # ----------------------------------------------------------------------------
  results_rpr <- tryCatch({
    run.jags(
      model = model_rpr,
      data = jags_data_rpr,
      monitor = monitor_params_rpr,
      adapt = 5000,
      burnin = 20000,
      sample = 20000,
      n.chains = n_chains_used,
      inits = inits_rpr,
      method = "parallel",
      summarise = FALSE   # Disable automatic summary computation
    )
  }, error = function(e) {
    cat("\nError in RPR analysis:", e$message, "\n")
    return(NULL)
  })
  
  # ----------------------------------------------------------------------------
  # RUN JAGS MODEL FOR CR ANALYSIS
  # ----------------------------------------------------------------------------
  results_cr <- tryCatch({
    run.jags(
      model = model_cr,
      data = jags_data_cr,
      monitor = monitor_params_cr,
      adapt = 5000,
      burnin = 20000,
      sample = 20000,
      n.chains = n_chains_used,
      inits = generate_inits(n_chains_used, processed_data, monitor_params_cr),
      method = "parallel",
      summarise = FALSE   # Disable automatic summary computation
    )
  }, error = function(e) {
    cat("\nError in CR analysis:", e$message, "\n")
    return(NULL)
  })
  
  results_combined <- list(
    rpr = if (!is.null(results_rpr)) as.mcmc.list(results_rpr) else NULL,
    cr = if (!is.null(results_cr)) as.mcmc.list(results_cr) else NULL,
    diagnostics = list(
      rpr_convergence = !is.null(results_rpr),
      cr_convergence = !is.null(results_cr),
      data_quality = attr(processed_data, "data_quality"),
      sample_sizes = attr(processed_data, "sample_sizes")
    )
  )
  
  return(results_combined)
}



# ==============================================================================
# FUNCTION: CHECK CONVERGENCE 
# Purpose: This function evaluates the convergence of MCMC chains by:
# 1. Filtering out parameters with zero variance.
# 2. Computing Geweke, Gelman-Rubin, and effective sample size diagnostics.
# 3. Returning convergence status for each diagnostic.
# ==============================================================================
check_convergence <- function(mcmc_object, var_threshold = 1e-6) {
  tryCatch({
    # Extract parameter names from the first chain
    params <- colnames(as.matrix(mcmc_object[[1]]))
    
    # Filter parameters: keep only those with variance greater than var_threshold
    nonzero_params <- sapply(params, function(param) {
      values <- unlist(lapply(mcmc_object, function(chain) chain[, param]))
      return(var(values, na.rm = TRUE) > var_threshold)
    })
    params_to_use <- params[nonzero_params]
    
    # If no parameter varies, return NA for all diagnostics.
    if (length(params_to_use) == 0) {
      return(list(geweke = NA, gelman = NA, n_eff = NA))
    }
    
    # Compute Geweke diagnostic for each parameter individually.
    geweke_results <- sapply(params_to_use, function(p) {
      tryCatch({
        diag <- coda::geweke.diag(mcmc_object[, p], autoburnin = FALSE)
        return(diag$z[1])
      }, error = function(e) { NA })
    })
    
    # Compute Gelman-Rubin diagnostic (PSRF) for each parameter individually.
    gelman_results <- sapply(params_to_use, function(p) {
      tryCatch({
        diag <- coda::gelman.diag(mcmc_object[, p], autoburnin = FALSE)
        return(diag$psrf[1])
      }, error = function(e) { NA })
    })
    
    # Compute effective sample size for each parameter individually.
    n_eff_values <- sapply(params_to_use, function(p) {
      tryCatch({
        return(coda::effectiveSize(mcmc_object[, p]))
      }, error = function(e) { NA })
    })
    
    # Remove NA values from each diagnostic vector.
    valid_idx <- !is.na(geweke_results) & !is.na(gelman_results) & !is.na(n_eff_values)
    
    # If no valid diagnostic values remain, return NA.
    if (sum(valid_idx) == 0) {
      return(list(geweke = NA, gelman = NA, n_eff = NA))
    }
    
    # Use only valid (non-NA) diagnostic values.
    geweke_vals <- geweke_results[valid_idx]
    gelman_vals <- gelman_results[valid_idx]
    n_eff_vals <- n_eff_values[valid_idx]
    
    # Overall convergence criteria:
    # - All Geweke z-scores must be within ±1.96.
    # - All Gelman-Rubin PSRF values must be below 1.1.
    # - All effective sample sizes must exceed 100.
    geweke_result <- all(abs(geweke_vals) < 1.96)
    gelman_result <- all(gelman_vals < 1.1)
    n_eff_result <- all(n_eff_vals > 100)
    
    return(list(geweke = geweke_result, gelman = gelman_result, n_eff = n_eff_result))
    
  }, error = function(e) {
    # In case of any unexpected error, return NA for all diagnostics.
    return(list(geweke = NA, gelman = NA, n_eff = NA))
  })
}

# ==============================================================================
# FUNCTION: LONG-TERM PREDICTION ANALYSIS 
# Purpose: Forecasts cumulative impact using a state-space model with:
# - Parallel computation for efficiency.
# - Bayesian posterior samples for uncertainty estimation.
# - Stability and sustainability assessment for long-term effects.
# ==============================================================================
generate_long_term_prediction <- function(posterior_matrix, analysis_type) {
  vars <- get_analysis_vars(analysis_type)  # Retrieve relevant variable names for the given analysis type
  
  tryCatch({
    # Extract key parameters from the posterior samples
    global_lift <- posterior_matrix[, vars$global_lift]
    sigma <- posterior_matrix[, vars$sigma]
    
    # Ensure the extracted samples are valid (non-NA, finite values)
    valid_indices <- !is.na(global_lift) & !is.na(sigma) & is.finite(global_lift) & is.finite(sigma)
    
    # If there are too few valid samples, stop execution
    if (sum(valid_indices) < 100) {
      stop("Insufficient valid posterior samples for reliable prediction")
    }
    
    # Filter valid samples
    global_lift <- global_lift[valid_indices]
    sigma <- sigma[valid_indices]
    
    # Define forecast parameters
    n_posterior_samples <- length(global_lift)  # Number of valid posterior samples
    forecast_horizon <- 12  # Forecast for the next 12 months
    
    # Generate noise components for the state-space model
    system_noise <- sigma / 4  # System noise component (scaled by sigma)
    state_noise <- matrix(rnorm(n_posterior_samples * forecast_horizon), 
                          nrow = n_posterior_samples, ncol = forecast_horizon) *
      matrix(rep(system_noise, each = forecast_horizon), nrow = n_posterior_samples)
    
    observation_noise <- sapply(1:forecast_horizon, function(t) {
      rnorm(n_posterior_samples, mean = 0, sd = sigma / sqrt(t))
    })
    
    # Generate decay rates from a Beta distribution
    decay_rates <- rbeta(n_posterior_samples, 1, 1)
    
    # ==========================================================================
    # VECTORIZED PREDICTION GENERATION (State-space model)
    # ==========================================================================
    predictions <- matrix(0, nrow = n_posterior_samples, ncol = forecast_horizon)
    
    for (t in 1:forecast_horizon) {
      multiplier0 <- (1 - decay_rates)^t  # Compute decay for the current time step
      
      if (t == 1) {
        noise_component <- state_noise[, 1]
      } else {
        multipliers <- sapply(1:t, function(j) (1 - decay_rates)^(t - j))
        noise_component <- rowSums(multipliers * state_noise[, 1:t, drop = FALSE])
      }
      
      # Compute state-space estimates for the current time step
      state_t <- multiplier0 * global_lift + noise_component
      pred_t <- state_t + observation_noise[, t]
      
      # Ensure predictions stay within reasonable bounds (-100% to +200%)
      predictions[, t] <- pmin(pmax(pred_t, -100), 200)
    }
    
    # ==========================================================================
    # COMPUTE HIGH POSTERIOR DENSITY (HPD) INTERVALS
    # ==========================================================================
    compute_hpd <- function(x, level = 0.95) {
      hpd <- tryCatch(coda::HPDinterval(as.mcmc(x), prob = level), error = function(e) NULL)
      if (!is.null(hpd)) {
        return(c(lower = hpd[1, 1], upper = hpd[1, 2]))
      } else {
        return(c(lower = quantile(x, 0.025), upper = quantile(x, 0.975)))
      }
    }
    
    hpd_results <- future.apply::future_apply(predictions, 2, compute_hpd)
    
    # Create a forecast dataframe
    forecast_df <- data.frame(
      Month = 1:forecast_horizon,
      Expected_Lift = colMeans(predictions),
      Lower_CI = hpd_results[1, ],
      Upper_CI = hpd_results[2, ]
    )
    
    # ==========================================================================
    # COMPUTE STABILITY AND SUSTAINABILITY METRICS
    # ==========================================================================
    
    # Compute trend stability based on weighted trajectory
    weights <- seq(0.5, 1, length.out = forecast_horizon)
    weighted_traj <- predictions * matrix(weights, nrow = n_posterior_samples, ncol = forecast_horizon, byrow = TRUE)
    
    # Assess long-term effect sustainability
    sign_consistency <- rowSums(weighted_traj > 0) == forecast_horizon
    magnitude_check <- rowSums(abs(weighted_traj) < 200) == forecast_horizon
    
    trend_stability <- future.apply::future_apply(weighted_traj, 1, function(x) {
      trend <- coef(lm(x ~ seq_along(x)))[2]  # Extract trend slope from linear model
      abs(trend) < mean(abs(x)) / 4  # Check if trend slope is small
    })
    
    prob_sustained <- mean(sign_consistency & magnitude_check & trend_stability)
    
    # Compute eigenvalues for stability score
    cov_matrix <- cov(predictions)
    eigen_values <- eigen(cov_matrix, symmetric = TRUE)$values
    stability_score <- 1 / (1 + sqrt(sum(eigen_values^2)))
    stability_score <- pmin(pmax(stability_score, 0), 1) * 100  # Normalize to [0, 100] range
    
    # Return computed forecast and stability metrics
    return(list(
      forecast = forecast_df,
      prob_sustained = prob_sustained,
      stability_score = stability_score,
      predictions = predictions
    ))
    
  }, error = function(e) {
    warning(sprintf("Error in long-term prediction for %s: %s", analysis_type, e$message))
    return(NULL)
  })
}


# ==============================================================================
# FUNCTION: BAYESIAN POWER ANALYSIS (
# Purpose: Compute statistical power for different sample sizes using Bayesian
# posterior distributions. The function is fully vectorized and leverages
# parallel computation via `future.apply` for efficiency.
# ==============================================================================

compute_bayesian_power <- function(posterior_matrix, analysis_type) {
  
  # Extract the global lift parameter from posterior samples
  mu <- posterior_matrix[, ifelse(analysis_type == "RPR", "global_lift_rpr", "global_lift_cr")]
  
  # Extract the between-study standard deviation
  sigma <- posterior_matrix[, "sigma_between"]
  
  # Identify valid samples (remove NA or infinite values)
  valid_idx <- !is.na(mu) & !is.na(sigma) & is.finite(mu) & is.finite(sigma)
  mu <- mu[valid_idx]
  sigma <- sigma[valid_idx]
  
  # Define different sample sizes to evaluate power
  sample_sizes <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  
  # ============================================================================
  # FUNCTION: Compute Power for a Given Sample Size
  # Purpose: Uses effect size and chi-squared test to determine power.
  # ============================================================================
  compute_power_for_n <- function(n) {
    
    # Compute effect size based on the analysis type
    if (analysis_type == "RPR") {
      effect_sizes <- abs(mu) / (sigma * sqrt(2 / n))
    } else { 
      # Compute probability values for conversion rates (CR analysis)
      p1 <- plogis(mu)
      p2 <- plogis(mu + sigma)
      effect_sizes <- 2 * abs(asin(sqrt(p1)) - asin(sqrt(p2)))  # Compute odds ratio transformation
    }
    
    # Compute non-centrality parameter for chi-squared distribution
    ncp <- effect_sizes^2 * n
    
    # Compute the power by comparing with chi-square threshold
    q_val <- qchisq(0.95, df = 1)  # 95% confidence threshold
    power_samples <- 1 - pchisq(q_val, df = 1, ncp = ncp)
    
    # Return mean power and confidence intervals
    list(
      mean_power = mean(power_samples),
      ci_lower = as.numeric(quantile(power_samples, 0.025)),  # 2.5% CI
      ci_upper = as.numeric(quantile(power_samples, 0.975))   # 97.5% CI
    )
  }
  
  # Compute power for all sample sizes using parallel processing
  power_results <- future.apply::future_lapply(sample_sizes, compute_power_for_n)
  
  # Create a dataframe to store power results
  power_df <- data.frame(
    Sample_Size = sample_sizes,
    Mean_Power = sapply(power_results, function(x) x$mean_power),
    CI_Lower = sapply(power_results, function(x) x$ci_lower),
    CI_Upper = sapply(power_results, function(x) x$ci_upper)
  )
  
  # Identify the minimum sample size required to reach 80% power
  min_n_80 <- min(sample_sizes[power_df$Mean_Power >= 0.8], Inf)
  
  # Compute power for the average sample size
  current_power <- compute_power_for_n(mean(sample_sizes))
  
  # Return power analysis results
  return(list(
    power_df = power_df,
    min_n_80 = min_n_80,
    current_power = current_power$mean_power,
    power_ci = c(current_power$ci_lower, current_power$ci_upper)
  ))
}


# ==============================================================================
# Helper function to get analysis-specific variables
# Purpose: Returns variable names for "RPR" or "CR" analyses.
# ==============================================================================
get_analysis_vars <- function(analysis_type) {
  vars <- list(
    RPR = list(
      global_lift = "global_lift_rpr",
      prob_positive = "prob_positive_rpr",
      sigma = "sigma_between",
      name = "REVENUE PER RECIPIENT (RPR)"
    ),
    CR = list(
      global_lift = "global_lift_cr",
      prob_positive = "prob_positive_cr",
      sigma = "sigma_between",
      name = "CONVERSION RATE (CR)"
    )
  )
  return(vars[[analysis_type]])  # Return the corresponding variable set
}


# ==============================================================================
# FUNCTION: CREATE AND PRINT FORMATTED RESULTS
# Purpose: Processes Bayesian meta-analysis results and prints them in 
# a structured format. Includes posterior summaries, credible intervals,
# heterogeneity analysis, and power calculations.
# ==============================================================================
create_and_print_results <- function(results, group_name) {
  cat("\n===========================================================")
  cat(sprintf("\nResults for: %s", group_name))
  cat("\n===========================================================\n")
  
  # Helper function to process and print analysis for each type (RPR or CR)
  process_analysis_type <- function(results, analysis_type) {
    result_key <- tolower(analysis_type)
    metrics <- list(
      RPR = list(
        name = "REVENUE PER RECIPIENT (RPR)",
        global_lift = "global_lift_rpr",
        prob_positive = "prob_positive_rpr",
        sigma = "sigma_between"
      ),
      CR = list(
        name = "CONVERSION RATE (CR)",
        global_lift = "global_lift_cr",
        prob_positive = "prob_positive_cr",
        sigma = "sigma_between"
      )
    )
    
    if (!is.null(results[[result_key]])) {
      metrics_vars <- metrics[[analysis_type]]
      cat(sprintf("\n%s ANALYSIS:\n", metrics_vars$name))
      
      cat("\nFULL POSTERIOR SUMMARY FOR", analysis_type, ":\n")
      
      # Convert the MCMC results to a matrix
      posterior_matrix <- as.matrix(results[[result_key]])
      
      # Compute the variance of each parameter using matrixStats if available for speed.
      if (requireNamespace("matrixStats", quietly = TRUE)) {
        param_variances <- matrixStats::colVars(posterior_matrix, na.rm = TRUE)
        names(param_variances) <- colnames(posterior_matrix)
      } else {
        param_variances <- apply(posterior_matrix, 2, function(x) var(x, na.rm = TRUE))
      }
      
      # Identify parameters that are nonconstant (variance greater than a small threshold)
      nonzero_params <- names(param_variances)[param_variances > 1e-10]
      
      # Warn if any parameters are omitted
      if (length(nonzero_params) < ncol(posterior_matrix)) {
        warning(sprintf("The following parameters have zero variance and will be omitted from the summary: %s",
                        paste(setdiff(colnames(posterior_matrix), nonzero_params), collapse = ", ")))
      }
      
      # If at least one parameter has variation, compute and print the summary; otherwise, print a message.
      if (length(nonzero_params) > 0) {
        filtered_matrix <- posterior_matrix[, nonzero_params, drop = FALSE]
        print(summary(filtered_matrix), digits = 4)
      } else {
        cat("All monitored parameters have zero variance. Summary statistics cannot be computed.\n")
      }
      
      # Extract statistics
      global_lift <- posterior_matrix[, metrics_vars$global_lift]
      sigma <- posterior_matrix[, metrics_vars$sigma]
      prob_positive <- mean(posterior_matrix[, metrics_vars$prob_positive])
      
      # Print main results dataframe
      results_df <- data.frame(
        Category = c(
          rep("Pooled Effect Metrics", 4),
          rep("Credible Intervals", 4),
          rep("Heterogeneity Analysis", 3)
        ),
        Metric = c(
          "Global Posterior Mean Lift (%)",
          "Median Lift (%)",
          "Standard Deviation",
          "Probability of Positive Effect (%)",
          "95% CI",
          "90% CI",
          "85% CI",
          "80% CI",
          "Between-study SD (σ)",
          "Heterogeneity CI",
          "Heterogeneity Magnitude"
        ),
        Value = c(
          sprintf("%.2f", mean(global_lift)),
          sprintf("%.2f", median(global_lift)),
          sprintf("%.2f", sd(global_lift)),
          sprintf("%.1f", prob_positive * 100),
          sprintf("[%.2f, %.2f]", quantile(global_lift, 0.025), quantile(global_lift, 0.975)),
          sprintf("[%.2f, %.2f]", quantile(global_lift, 0.05), quantile(global_lift, 0.95)),
          sprintf("[%.2f, %.2f]", quantile(global_lift, 0.075), quantile(global_lift, 0.925)),
          sprintf("[%.2f, %.2f]", quantile(global_lift, 0.1), quantile(global_lift, 0.9)),
          sprintf("%.2f", mean(sigma)),
          sprintf("[%.2f, %.2f]", quantile(sigma, 0.025), quantile(sigma, 0.975)),
          ifelse(mean(sigma) < 10, "Low",
                 ifelse(mean(sigma) < 30, "Moderate", "High"))
        )
      )
      print(results_df)
      
      # Generate and print long-term prediction analysis
      cat(sprintf("\n📈 %s LONG-TERM PREDICTION ANALYSIS:\n", analysis_type))
      pred_results <- generate_long_term_prediction(posterior_matrix, analysis_type)
      stability_score <- pred_results$stability_score
      
      # Print forecast
      if (!is.null(pred_results$forecast) && !all(is.na(pred_results$forecast))) {
        print(pred_results$forecast)
        
        cat("\nLong-term Effect Analysis:\n")
        if (!is.na(pred_results$stability_score)) {
          cat(sprintf("• Stability Score: %.1f%%\n", pred_results$stability_score))
          cat(sprintf("• Probability of Sustained Effect: %.1f%%\n", 
                      pred_results$prob_sustained * 100))
          
          cat("\nStability Assessment:\n")
          if (pred_results$stability_score > 80) {
            cat("✅ High stability - predictions are reliable\n")
          } else if (pred_results$stability_score > 50) {
            cat("🟡 Moderate stability - interpret with caution\n")
          } else {
            cat("⚠️ Low stability - predictions may be unreliable\n")
          }
        }
      } else {
        cat("\n⚠️ Unable to generate reliable long-term predictions\n")
      }
      
      # Add convergence diagnostics
      conv_checks <- check_convergence(results[[result_key]])
      
      cat("\n📊 Convergence Summary:\n")
      if (!is.na(conv_checks$geweke)) {
        cat(sprintf("  - Geweke Test: %s\n", 
                    ifelse(conv_checks$geweke, "✓ Passed", "✗ Failed")))
      }
      if (!is.na(conv_checks$gelman)) {
        cat(sprintf("  - Gelman-Rubin Test: %s\n", 
                    ifelse(conv_checks$gelman, "✓ Passed", "✗ Failed")))
      }
      if (!is.na(conv_checks$n_eff)) {
        cat(sprintf("  - Effective Adaptation Sample Size: %s\n",
                    ifelse(conv_checks$n_eff, "✓ Adequate", "✗ Inadequate")))
      }
      if (all(is.na(c(conv_checks$geweke, conv_checks$gelman, conv_checks$n_eff)))) {
        cat("\n⚠️ Warning: Convergence checks failed. Results may be unreliable.\n")
        cat("  Consider increasing adaptation and sample size or checking model specification.\n")
      }
      
      # Data sufficiency assessment
      cat(sprintf("\n📊 %s DATA SUFFICIENCY ASSESSMENT:\n", analysis_type))
      ci_width <- quantile(global_lift, 0.975) - quantile(global_lift, 0.025)
      if (ci_width > 50) {
        cat(sprintf("\n⚠️ Warning: Wide %s credible interval (%.2f) indicates high uncertainty\n", 
                    analysis_type, ci_width))
        cat("• Consider collecting more data or running additional tests\n")
      }
      
      # Power Analysis Section
      cat(sprintf("\n📊 %s POWER ANALYSIS:\n", analysis_type))
      power_results <- compute_bayesian_power(posterior_matrix, analysis_type)
      
      # Print power analysis table
      power_df <- power_results$power_df
      power_df$Power <- sprintf("%.1f%% [%.1f%%, %.1f%%]",
                                power_df$Mean_Power * 100,
                                power_df$CI_Lower * 100,
                                power_df$CI_Upper * 100)
      power_df$CI_Lower <- NULL
      power_df$CI_Upper <- NULL
      power_df$Mean_Power <- NULL
      
      cat("\nPower Analysis Results:\n")
      print(power_df)
      
      # Print power analysis conclusions
      cat("\nPower Analysis Summary:\n")
      if (power_results$min_n_80 < Inf) {
        cat(sprintf("• Minimum sample size for 80%% power: %d\n", 
                    power_results$min_n_80))
      } else {
        cat("• Unable to achieve 80% power with evaluated sample sizes\n")
      }
      
      cat(sprintf("• Current achieved power: %.1f%% [%.1f%%, %.1f%%]\n",
                  power_results$current_power * 100,
                  power_results$power_ci[1] * 100,
                  power_results$power_ci[2] * 100))
      
      # Power assessment
      cat("\nPower Assessment:\n")
      if (power_results$current_power >= 0.8) {
        cat("✅ Test is well-powered\n")
      } else if (power_results$current_power >= 0.5) {
        cat("🟡 Test is moderately powered\n")
        cat("   Consider increasing sample size if possible\n")
      } else {
        cat("⚠️ Test is underpowered\n")
        cat("   Recommendations:\n")
        cat("   - Increase sample size\n")
        cat("   - Consider pooling data from similar tests\n")
        cat("   - Review test design for potential improvements\n")
      }
      
      # Print final recommendations
      cat(sprintf("\n📋 %s FINAL RECOMMENDATIONS:\n", analysis_type))
      if (prob_positive > 0.95 && pred_results$stability_score > 70) {
        cat("✅ Strong evidence for positive effect - Implement change\n")
      } else if (prob_positive > 0.8 && pred_results$stability_score > 50) {
        cat("🟡 Moderate evidence for positive effect - Consider implementation\n")
      } else if (prob_positive < 0.2 || pred_results$stability_score < 30) {
        cat("❌ Evidence suggests negative effect or unstable results - Do not implement\n")
      } else {
        cat("⚠️ Results inconclusive - Consider additional testing\n")
      }
    }
  }
  
  # Helper function for printing diagnostics remains unchanged
  print_diagnostics <- function(diagnostics) {
    cat("\nConvergence Summary:\n")
    cat(sprintf("• Total Parameters: %d\n", diagnostics$summary$total_params))
    cat(sprintf("• Converged Parameters: %d\n", diagnostics$summary$converged_params))
    
    if (length(diagnostics$zero_var_params) > 0) {
      cat("\n⚠️ Parameters with zero variance:\n")
      cat(paste("  -", diagnostics$zero_var_params, collapse = "\n"), "\n")
    }
    
    if (length(diagnostics$non_converged_params) > 0) {
      cat("\n⚠️ Parameters with convergence issues:\n")
      cat(paste("  -", diagnostics$non_converged_params, collapse = "\n"), "\n")
    }
    
    convergence_ratio <- diagnostics$summary$converged_params / 
      diagnostics$summary$total_params
    
    cat("\nOverall Model Assessment:\n")
    if (convergence_ratio >= 0.9) {
      cat("✅ Model shows good convergence\n")
    } else if (convergence_ratio >= 0.7) {
      cat("🟡 Model shows acceptable convergence\n")
    } else {
      cat("❌ Model shows significant convergence issues\n")
    }
  }
  
  # Process both RPR and CR analyses
  process_analysis_type(results, "RPR")
  process_analysis_type(results, "CR")
  
  cat("\n=================================================================\n")
}

# ==============================================================================
# RUN ANALYSIS
# Purpose: This section scans the environment for valid datasets, filters them,
#          runs Bayesian meta-analysis, and prints results.
# ==============================================================================

# Initialize an empty list to store valid datasets
datasets <- list()

# Retrieve all object names stored in the `analysis_env` environment
section_objects <- ls(envir = analysis_env)

# Define the required columns that each dataset must contain for analysis
required_cols <- c("revenue_control", "revenue_variation", 
                   "conversions_control", "conversions_variation",
                   "recipients_control", "recipients_variation")

# ==============================================================================
# Identify and Validate Datasets
# ==============================================================================
for (obj_name in section_objects) {
  # Retrieve the object from the environment
  obj <- get(obj_name, envir = analysis_env)
  
  # Check if the object is a valid data frame and contains all required columns
  if (is.data.frame(obj) && all(required_cols %in% colnames(obj))) {
    datasets[[obj_name]] <- obj  # Store the valid dataset
  }
}

# ==============================================================================
# Handle Case Where No Valid Datasets Are Found
# ==============================================================================
if (length(datasets) == 0) {
  cat("\nNo valid datasets were defined in this section.\n")
  cat("Datasets must be data frames with the following required columns:\n")
  cat(paste("-", required_cols), sep = "\n")  # Print missing column names
  
} else {
  # ============================================================================
  # Print the Number of Valid Datasets Found and Their Names
  # ============================================================================
  cat("\nAnalyzing", length(datasets), "dataset(s) from current section:\n")
  cat(paste("-", names(datasets)), sep = "\n")  # Print dataset names
  
  # ============================================================================
  # Run Bayesian Meta-Analysis on Each Dataset
  # ============================================================================
  for (dataset_name in names(datasets)) {
    results <- run_bayesian_meta_analysis(datasets[[dataset_name]], dataset_name)  # Run analysis
    create_and_print_results(results, dataset_name)  # Generate and print results
  }
}

# ==============================================================================
# Clean Up: Remove Analysis Environment to Free Up Memory
# ==============================================================================
rm(analysis_env)

