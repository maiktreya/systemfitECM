###############################################
# Define dataset of usage (data.table required) and selected variables for coint. analysis. The dependant var should be listed first.

library(data.table) # for simple and performant data manipulation
library(plm) # needed for systemfit to handle panel structure
library(systemfit) # for FGLS system linear models
library(magrittr) # For piping with %>% without dplyr dependencies
library(aod) # for performing F Bounds test

# install and import this library
library(systemfitECM)

# Create the sample dataset
set.seed(1234) # For reproducibility
countries <- c("Austria", "Germany", "Italy")
period <- 1992:2019
table_dt <- data.table(
    reporter = rep(countries, each = length(period)),
    year = rep(period, length(countries)),
    tech_exports = rnorm(length(countries) * length(period), 5000, 1000), # Sample tech_exports data
    rprices = rnorm(length(countries) * length(period), 100, 20), # Sample rprices data
    fincome = rnorm(length(countries) * length(period), 40000, 5000), # Sample fincome data
    investment = rnorm(length(countries) * length(period), 40000, 5000), # Sample rprices data
    consumption = rnorm(length(countries) * length(period), 40000, 5000) # Sample fincome data
)

# Set remaining control parameters
sel_variables <- c("tech_exports", "rprices", "fincome") # first is dependant variable in systemfit
instruments <- c("fincome", "investment", "consumption") # first is endogenous regressor and remaining their instruments.
method <- "3SLS"
estimation3SLS <- "EViews"
lags <- 1
iterations <- 1

# Get an Unrestricted ECM using systemfit methods
pre_exp <- uecm_systemfit(
    dt = table_dt,
    col_names = sel_variables,
    nlags = lags,
    grouping = "reporter",
    method = method,
    iterations = iterations,
    method_solv = estimation3SLS, # only 3sls,
    inst_list = instruments # endo first, then remaining
)
pre_exp %>%
    summary() %>%
    print()

# Apply and F Bound-Test for equations in systems following Pesaran (2001)
bounds_F_results <- systemfit_boundsF_test(
    system_ecm = pre_exp,
    units = countries,
    sel_variables = sel_variables
)
bounds_F_results %>% print()

# Finally, get a Restricted ECM using systemfit methods
pos_exp <- recm_systemfit(
    uecm_model = pre_exp,
    dt = table_dt,
    col_names = sel_variables,
    nlags = lags,
    grouping = "reporter",
    method = method,
    iterations = iterations,
    nunits = length(countries),
    nperiods = length(period),
    method_solv = estimation3SLS, # only 3sls,
    inst_list = instruments # endo first, then remaining
)
pos_exp %>%
    summary() %>%
    print()


aa <- get_ect_systemfit(
    systemfit_uecm_coefs = pre_exp,
    nperiods = length(period),
    nunits = length(countries),
    sel_variables = sel_variables,
    table_dt = table_dt
)
