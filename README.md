# 🎯 Bayesian Meta-Analysis for A/B Testing

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![JAGS](https://img.shields.io/badge/JAGS-4.3+-green.svg)](https://mcmc-jags.sourceforge.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/Status-Production%20Ready-brightgreen.svg)]()
[![Client Rating](https://img.shields.io/badge/Client%20Rating-5.0%2F5.0-gold.svg)]()

> **A production-ready Bayesian meta-analysis framework for email marketing A/B testing with hierarchical modeling, uncertainty quantification, and long-term predictions.**

## 🚀 Overview

This repository contains a sophisticated Bayesian meta-analysis system designed specifically for email marketing A/B testing. The framework enables data scientists and marketing analysts to combine results from multiple A/B tests to estimate pooled effects with proper uncertainty quantification.

### 🎯 Key Features

- **📊 Hierarchical Bayesian Modeling**: Random-effects meta-analysis using JAGS
- **🔍 Dual Metric Analysis**: Revenue Per Recipient (RPR) and Conversion Rate (CR)
- **📈 Uncertainty Quantification**: Multiple credible intervals (95%, 90%, 85%, 80%)
- **🔮 Long-term Predictions**: State-space modeling for future performance
- **⚡ Parallel Processing**: Optimized for multi-core computation
- **📋 Comprehensive Diagnostics**: Convergence checks and power analysis
- **🎨 Production-Ready**: Robust error handling and edge case management

## 📋 Table of Contents

- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Features](#-features)
- [Data Format](#-data-format)
- [Usage Examples](#-usage-examples)
- [Model Specifications](#-model-specifications)
- [Output Interpretation](#-output-interpretation)
- [Performance](#-performance)
- [Contributing](#-contributing)
- [License](#-license)

## 🛠️ Installation

### Prerequisites

- R version 4.0 or higher
- JAGS 4.3.0 or higher

### Quick Installation

```r
# Install JAGS first (system-specific)
# Windows: Download from https://mcmc-jags.sourceforge.io/
# macOS: brew install jags
# Ubuntu: sudo apt-get install jags

# Clone and run the script - it will auto-install required R packages
git clone https://github.com/yourusername/bayesian-meta-analysis-ab-testing.git
cd bayesian-meta-analysis-ab-testing
Rscript Bayesian_Meta_Analysis_Optimized-cauchy.R
```

### Manual Package Installation

```r
required_packages <- c("runjags", "coda", "dplyr", "tidyr", "ggplot2", 
                      "parallel", "future", "future.apply")
install.packages(required_packages, dependencies = TRUE)
```

## 🚀 Quick Start

1. In GitHub, go to `code/Bayesian_Meta_Analysis_Optimized-cauchy.R`  
2. Click the green **Code** button → **Download ZIP** to grab the script  
   (or click **Raw** to copy-paste).  
3. In R, run:
   ```r
   source("Bayesian_Meta_Analysis_Optimized-cauchy.R")
   results <- run_bayesian_meta_analysis(my_data, "My Campaign")


```r
# Load the script
source("Bayesian_Meta_Analysis_Optimized-cauchy.R")

# Prepare your data
my_data <- data.frame(
  revenue_control = c(1406, 1105, 0, 1347, 648),
  revenue_variation = c(0, 897, 1980, 0, 0),
  conversions_control = c(1, 1, 0, 3, 1),
  conversions_variation = c(0, 2, 2, 0, 0),
  recipients_control = c(5758, 14231, 5142, 10147, 5472),
  recipients_variation = c(5749, 14192, 5121, 10154, 5440)
)

# Run analysis
results <- run_bayesian_meta_analysis(my_data, "My Campaign")
create_and_print_results(results, "My Campaign")
```

## 🎯 Features

### 📊 Advanced Statistical Modeling

- **Random Effects Meta-Analysis**: Accounts for between-study heterogeneity
- **Non-informative Priors**: Unbiased parameter estimation
- **Robust Likelihood Functions**: Handles zero values and edge cases
- **Hierarchical Structure**: Properly models study-level and global effects

### 🔍 Comprehensive Output

```r
# What you get from each analysis:
✅ Global posterior mean lift (%)
✅ Multiple credible intervals (95%, 90%, 85%, 80%)
✅ Probability of positive effect
✅ Between-study heterogeneity (σ)
✅ Convergence diagnostics
✅ Power analysis
✅ Long-term predictions
✅ Stability assessment
```

### ⚡ Performance Optimizations

- **Parallel MCMC Chains**: Leverages all available CPU cores
- **Efficient Memory Management**: Optimized for large datasets
- **Vectorized Computations**: Fast matrix operations
- **Smart Initialization**: Improves convergence speed

## 📊 Data Format

Your data must include these columns:

| Column | Description | Example |
|--------|-------------|---------|
| `revenue_control` | Revenue from control group | 1406 |
| `revenue_variation` | Revenue from variation group | 897 |
| `conversions_control` | Conversions in control | 1 |
| `conversions_variation` | Conversions in variation | 2 |
| `recipients_control` | Recipients in control | 5758 |
| `recipients_variation` | Recipients in variation | 5749 |

### 📁 Sample Datasets Included

- **Low Conversion Client**: 10 email tests with sparse conversion data
- **Fishing Gear Company**: 7 tests with moderate conversion rates
- **SMS Tests**: 4 tests including major outliers

## 💡 Usage Examples

### Basic Analysis

```r
# Run analysis on built-in dataset
results <- run_bayesian_meta_analysis(
  analysis_env$mattress_low_conversion_1, 
  "Mattress Campaign"
)
```

### Custom Analysis with Validation

```r
# Your custom dataset
custom_data <- data.frame(
  revenue_control = c(1000, 1500, 2000),
  revenue_variation = c(1200, 1400, 2300),
  conversions_control = c(10, 15, 20),
  conversions_variation = c(12, 14, 23),
  recipients_control = c(5000, 6000, 7000),
  recipients_variation = c(5000, 6000, 7000)
)

# Preprocess and validate
processed_data <- preprocess_data(custom_data)
results <- run_bayesian_meta_analysis(processed_data, "Custom Analysis")
```

## 🔬 Model Specifications

### Revenue Per Recipient (RPR) Model

```r
# Hierarchical structure:
# Level 1: Revenue ~ Gamma(shape, rate)
# Level 2: log(μ) ~ Normal(θ, τ)
# Level 3: θ ~ Normal(μ_global, τ_between)
```

### Conversion Rate (CR) Model

```r
# Hierarchical structure:
# Level 1: Conversions ~ Binomial(n, p)
# Level 2: logit(p) ~ Normal(μ, τ)
# Level 3: μ ~ Normal(μ_global, τ_between)
```

### Prior Specifications

- **Global mean**: `μ_global ~ Normal(0, 0.001)`
- **Between-study SD**: `σ_between ~ HalfCauchy(0, 10)`
- **Study effects**: `θ_i ~ Normal(μ_global, τ_between)`

## 📈 Output Interpretation

### Key Metrics

| Metric | Interpretation |
|--------|----------------|
| **Global Posterior Mean** | Average effect across all studies |
| **95% Credible Interval** | Range containing true effect with 95% probability |
| **Probability of Positive Effect** | Likelihood that the intervention is beneficial |
| **Between-study SD (σ)** | Heterogeneity between studies |
| **Stability Score** | Reliability of long-term predictions |

### Decision Framework

| Condition | Recommendation |
|-----------|----------------|
| P(positive) > 95% & Stability > 70% | ✅ **Implement change** |
| P(positive) > 80% & Stability > 50% | 🟡 **Consider implementation** |
| P(positive) < 20% \| Stability < 30% | ❌ **Do not implement** |
| Otherwise | ⚠️ **Additional testing needed** |

## 🎨 Visualization Examples

The framework generates comprehensive diagnostic plots and summaries:

```r
# Convergence diagnostics
✓ Geweke Test: Passed
✓ Gelman-Rubin Test: Passed  
✓ Effective Sample Size: Adequate

# Power analysis for different sample sizes
Sample_Size | Mean_Power | CI_Lower | CI_Upper
100         | 23.5%      | 18.2%    | 29.1%
500         | 67.8%      | 61.4%    | 73.9%
1000        | 89.2%      | 85.1%    | 92.8%
```

## 🚀 Performance

- **Parallel Processing**: Utilizes all available CPU cores
- **Memory Efficient**: Optimized for large datasets (1000+ studies)
- **Fast Convergence**: Typically converges in <2 minutes for 10 studies
- **Robust**: Handles edge cases (zero conversions, missing data)

## 🎯 Real-World Applications

This framework has been successfully used for:

- **Email Marketing Optimization**: A/B testing subject lines, send times, content
- **E-commerce Conversion**: Product page layouts, checkout processes
- **Digital Advertising**: Ad creative performance, audience targeting
- **Product Development**: Feature rollouts, UI/UX changes

## 📊 Client Success Story

> *"Philip communicated well, and clearly, and performed the task in a timely manner. I needed a hierarchical Bayesian meta-analysis script for R, and Philip created one... He did really well within the scope of the project. Philip also strikes me as someone who cares about what he's doing."*
> 
> **⭐⭐⭐⭐⭐ 5.0/5.0** - Client Review

## 🛠️ Advanced Features

### Long-term Prediction Engine

```r
# Generates 12-month forecasts using state-space modeling
forecast_results <- generate_long_term_prediction(posterior_samples, "RPR")

# Includes:
- Expected lift trajectory
- Confidence intervals
- Stability assessment
- Sustainability metrics
```

### Power Analysis Suite

```r
# Computes statistical power for different sample sizes
power_results <- compute_bayesian_power(posterior_samples, "CR")

# Provides:
- Minimum sample size for 80% power
- Current achieved power
- Power curves for planning
```

## 🔧 Customization Options

### Model Configuration

```r
# Adjust MCMC parameters
adapt = 5000      # Adaptation iterations
burnin = 20000    # Burn-in iterations  
sample = 20000    # Posterior samples
n.chains = 4      # Number of parallel chains
```

### Prior Specifications

```r
# Modify priors for domain knowledge
mu_global ~ dnorm(domain_mean, domain_precision)
sigma_between ~ dt(0, scale, 1) T(0,)  # Half-t prior
```

## 📝 Contributing

We welcome contributions!
### Development Setup

```bash
git clone https://github.com/yourusername/bayesian-meta-analysis-ab-testing.git
cd bayesian-meta-analysis-ab-testing
```

### Running Tests

```r
# Run validation tests
source("tests/validation_tests.R")
run_all_tests()
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Professional Services

Looking for custom Bayesian analysis solutions? I specialize in:

- **Custom Meta-Analysis Frameworks**
- **Bayesian A/B Testing Systems**
- **Marketing Analytics Pipelines**
- **Statistical Consulting & Training**

### 📞 Contact

- **Upwork**: [https://www.upwork.com/freelancers/~01055a2b89788071d4?mp_source=share]
- **Portfolio**: [https://apollop24.github.io/]

## 🎯 Why Choose This Framework?

- ✅ **Production-Tested**: Successfully deployed in real marketing campaigns
- ✅ **Scientifically Rigorous**: Proper Bayesian methodology with peer-reviewed techniques
- ✅ **User-Friendly**: Comprehensive documentation and examples
- ✅ **Scalable**: Handles everything from small pilot tests to large-scale analyses
- ✅ **Maintained**: Regular updates and improvements based on client feedback

---

<div align="center">

### 🌟 Star this repository if you find it useful!

*Built with ❤️ by a data scientist who cares about statistical rigor and practical applications.*

</div>

## 🔄 Recent Updates

- **v1.2.0**: Added long-term prediction engine with state-space modeling
- **v1.1.0**: Implemented parallel processing for faster computation
- **v1.0.0**: Initial release with full Bayesian meta-analysis framework


---

*Last Updated: July 2025*
