############################################################################################ 3
# Example usage:
# Define dataset of usage (data.table required) and selected variables for coint. analysis. The dependant var should be listed first.
# Ensure the data.table package is loaded
library(data.table)
library(plm)
library(systemfit)
library(magrittr)

# Create the sample dataset
set.seed(123) # For reproducibility
countries <- c("Austria", "Germany", "Italy")
table_dt <- data.table(
    reporter = rep(countries, each = 28),
    year = rep(1992:2019, 3),
    tech_exports = rnorm(30, 5000, 1000), # Sample tech_exports data
    rprices = rnorm(30, 100, 20), # Sample rprices data
    fincome = rnorm(30, 40000, 5000) # Sample fincome data
)
sel_variables <- c("tech_exports", "rprices", "fincome")
lags <- 1
iterations <- 1

# Get an Unrestricted ECM using systemfit methods
pre_exp <- uecm_systemfit(
    dt = table_dt,
    col_names = sel_variables,
    nlags = lags,
    method = "SUR",
    iterations = iterations,
    method_solv = "EViews" # only 3sls
)

# Get a Restricted ECM using systemfit methods
pos_exp <- recm_systemfit(
    uecm_model = pre_exp,
    dt = table_dt,
    col_names = sel_variables,
    nlags = lags,
    method = "SUR",
    iterations = iterations,
    method_solv = "EViews" # only 3sls
) %>%
    summary() %>%
    print()

# Finally, apply and F Bound-Test for equations in systems following Pesaran (2001)
bounds_F_results <- systemfit_boundsF_test(
    system_ecm = pos_exp,
    units = countries
) %>% print()
