###############################################
# Define dataset of usage (data.table required) and selected variables for coint. analysis. The dependant var should be listed first.

library(data.table) # for simple and performant data manipulation
library(plm) # needed for systemfit to handle panel structure
library(systemfit) # for FGLS system linear models
library(magrittr) # For piping with %<% without dplyr dependencies
library(aod) # for performing F Bounds test

# install and import this library
library(systemfitECM)

# Create the sample dataset
set.seed(123) # For reproducibility
countries <- c("Austria", "Germany", "Italy")
period <- 1992:2019
table_dt <- data.table(
    reporter = rep(countries, each = length(period)),
    year = rep(period, length(countries)),
    tech_exports = rnorm(length(countries) * length(period), 5000, 1000), # Sample tech_exports data
    rprices = rnorm(length(countries) * length(period), 100, 20), # Sample rprices data
    fincome = rnorm(length(countries) * length(period), 40000, 5000) # Sample fincome data
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

# Apply and F Bound-Test for equations in systems following Pesaran (2001)
bounds_F_results <- systemfit_boundsF_test(
    system_ecm = pre_exp,
    units = countries
) %>% print()

# Finally, get a Restricted ECM using systemfit methods
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
