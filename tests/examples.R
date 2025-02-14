# systemfitECM/tests/testthat/test_systemfitECM.R #nolint

library(testthat)
library(systemfitECM)
library(data.table)
library(plm)
library(aod) # For wald.test

# Test Data (adapted from your script)
set.seed(1234)
countries <- c("Austria", "Germany", "Italy")
period <- 1992:2019
table_dt <- data.table(
    reporter = rep(countries, each = length(period)),
    year = rep(period, length(countries)),
    tech_exports = rnorm(length(countries) * length(period), 5000, 1000),
    rprices = rnorm(length(countries) * length(period), 100, 20),
    fincome = rnorm(length(countries) * length(period), 40000, 5000),
    investment = rnorm(length(countries) * length(period), 40000, 5000),
    consumption = rnorm(length(countries) * length(period), 40000, 5000)
)


sel_variables <- c("tech_exports", "rprices", "fincome")
instruments <- c("tech_exports", "fincome", "investment", "consumption") # tech_exports is now also instrument
method <- "SUR" # Or "2SLS", "3SLS" as needed
estimation3SLS <- "EViews" # If method is 3SLS
lags <- 2
iterations <- 1
nunits <- length(countries)
nperiods <- length(period)



test_that("uecm_systemfit works", {
    uecm_model <- uecm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        method_solv = estimation3SLS,
        inst_list = instruments
    )
    expect_s3_class(uecm_model, "systemfit")
    expect_named(
        coef(uecm_model),
        c(
            "(Intercept)",
            paste0(rep(sel_variables, each = lags + 1), "_lag", rep(0:lags, times = length(sel_variables))),
            paste0(sel_variables[-1], "_diff")
        )
    ) # Check if coefficients have the correct names. Adjust to your specific names.

    # More specific tests (e.g., check coefficient values if you have expected values)
    # Example:
    # expect_equal(coef(uecm_model)[["rprices_lag1"]], expected_value, tolerance = 0.01) #nolint
})


test_that("get_ect_systemfit works", {
    uecm_model <- uecm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        method_solv = estimation3SLS,
        inst_list = instruments
    )
    ect_test <- get_ect_systemfit(
        systemfit_uecm_coefs = uecm_model,
        nperiods = nperiods,
        nunits = nunits,
        sel_variables = sel_variables,
        table_dt = table_dt
    )
    expect_s3_class(ect_test, "data.table")
    expect_equal(nrow(ect_test), nunits * (nperiods + 1)) # Check correct number of rows
    expect_named(ect_test, c("key", "time", "ect_x")) # Check column names

    # Add more specific tests for the ECT values if you have expected values.
})


test_that("recm_systemfit works", {
    uecm_model <- uecm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        method_solv = estimation3SLS,
        inst_list = instruments
    )

    recm_model <- recm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        uecm_model = uecm_model,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        nunits = nunits,
        nperiods = nperiods,
        method_solv = estimation3SLS,
        inst_list = instruments
    )
    expect_s3_class(recm_model, "systemfit")
    # Add more specific tests (e.g., check coefficient values, compare with UECM if needed)
})


test_that("systemfit_boundsF_test works", {
    uecm_model <- uecm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        method_solv = estimation3SLS,
        inst_list = instruments
    )

    bounds_F_results <- systemfit_boundsF_test(
        system_ecm = uecm_model,
        units = countries,
        sel_variables = sel_variables
    )
    expect_equal(length(bounds_F_results), length(countries)) # One result per country
    expect_named(bounds_F_results, countries) # Results should be named by country

    # Add more specific tests if you have expected F-statistic values.
})


test_that("Error handling in uecm_systemfit", {
    expect_error(
        uecm_systemfit(data = table_dt, col_names = sel_variables, nlags = -1, grouping = "reporter", method = method, inst_list = instruments),
        "nlags must be non-negative"
    )
    expect_error(
        uecm_systemfit(data = table_dt, col_names = sel_variables, nlags = lags, grouping = "wrong_group", method = method, inst_list = instruments),
        "grouping variable wrong_group not found in data"
    )
    expect_error(
        uecm_systemfit(data = table_dt, col_names = sel_variables, nlags = lags, grouping = "reporter", method = "2SLS"), # Missing inst_list
        "inst_list must be provided for 2SLS/3SLS"
    )

    # Add more error condition tests as needed (e.g. non-existent variable names, etc.)
})

test_that("Error handling in recm_systemfit", {
    uecm_model <- uecm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        method_solv = estimation3SLS,
        inst_list = instruments
    )

    expect_error(
        recm_systemfit(
            data = table_dt,
            col_names = sel_variables,
            uecm_model = uecm_model,
            nlags = -1,
            grouping = "reporter",
            method = method,
            nunits = nunits,
            nperiods = nperiods,
            inst_list = instruments
        ),
        "nlags must be non-negative"
    )
    # ... Add more error condition tests as needed
})

test_that("Error handling in get_ect_systemfit", {
    uecm_model <- uecm_systemfit(
        data = table_dt,
        col_names = sel_variables,
        nlags = lags,
        grouping = "reporter",
        method = method,
        iterations = iterations,
        method_solv = estimation3SLS,
        inst_list = instruments
    )
    sel_variables_wrong <- c("wrong_variable", sel_variables[-1])

    expect_error(
        get_ect_systemfit(systemfit_uecm_coefs = uecm_model, sel_variables = sel_variables_wrong, table_dt = table_dt, nperiods = nperiods, nunits = nunits),
        "Missing coefficients for variables: wrong_variable"
    )
})
