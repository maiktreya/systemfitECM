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
table_dt <- fread(".data/eudata_final_nom.csv")[tech == "HIGH", ]

# Set remaining control parameters
sel_variables <- c("tech_imports", "rprices", "fincome") # first is dependant variable in systemfit
instruments <- c("fincome", "investment", "consump") # first is endogenous regressor and remaining their instruments.
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
    method_solv = estimation3SLS, # only 3sls,
    inst_list = instruments # endo first, then remaining
)
pos_exp %>%
    summary() %>%
    print()
