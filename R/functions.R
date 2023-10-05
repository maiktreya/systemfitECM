#### MODULE FOR INTEGRATION OF ECM INTO SYSTEMFIT

#' Unrestricted Error Correction Model for Systemfit
#' @importFrom rlang :=
#' Computes the Unrestricted Error Correction Model for Systemfit.
#'
#' @param col_names A character vector of column names.
#' @param nlags An integer specifying the number of lags.
#' @param method Character string indicating the desired estimation method.
#' @param method_solv Character string indicating the solution method. Default is "EViews".
#' @param iterations An integer indicating the number of iterations.
#' @param dt A data.table object containing the data.
#' @param inst_list List of instruments for 2SLS and 3SLS.
#'
#' @return A model result from systemfit.
#' @export
uecm_systemfit <- function(
    col_names = c(),
    nlags = 1,
    grouping,
    method = "SUR",
    method_solv = "EViews",
    iterations = 1,
    dt = data.table::data.table(),
    inst_list = c()) {
    diff_cols <- c()
    all_lag_cols <- c()

    ifelse(method != "SUR", col_names_ext <- c(col_names, inst_list[-1]), col_names_ext <- col_names)

    for (col in col_names_ext) {
        # Add diff column
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := diff(c(NA, get(col))), by = get(grouping)]
        diff_cols <- c(diff_cols, diff_col) # Populate diff_cols vector
    }
    for (col in col_names) {
        # Add lag columns for each lag value
        for (lag in 1:nlags) {
            lag_col <- paste0(col, "_lag", lag)
            dt[, (lag_col) := data.table::shift(get(col), n = lag, type = "lag"), by = get(grouping)]
            all_lag_cols <- c(all_lag_cols, lag_col) # Populate all_lag_cols vector
        }
    }

    if (method != "SUR") {
        diff_inst <- diff_cols[!(diff_cols %like% inst_list[1])]
        inst_eq <- paste("~", paste(c(diff_inst[-1], all_lag_cols), collapse = " + "))
        diff_cols <- diff_cols[!(diff_cols %like% inst_list[-1])]
    }

    # Construct formula strings
    formula_str <- paste(diff_cols[1], "~", paste(c(diff_cols[-1], all_lag_cols), collapse = " + "))

    # Remove rows with NA values
    dt <- dt[complete.cases(dt), ]
    dt <- plm::pdata.frame(dt, index = c("reporter", "year"))

    # Run systemfit model
    if (method == "3SLS") {
        control_system <- systemfit::systemfit.control(
            methodResidCov = "noDfCor",
            residCovWeighted = FALSE,
            maxiter = iterations,
            tol = 1e-5,
            method3sls = "EViews" # GLS(default), IV, GMM, SCHMIDT, EVIEWS
        )
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system, inst = as.formula(inst_eq))
    }
    if (method == "2SLS") {
        control_system <- systemfit::systemfit.control(
            methodResidCov = "noDfCor",
            residCovWeighted = FALSE,
            tol = 1e-5,
        )
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system, inst = as.formula(inst_eq))
    }
    if (method == "SUR") {
        control_system <- systemfit::systemfit.control(
            methodResidCov = "noDfCor",
            residCovWeighted = FALSE,
            maxiter = iterations,
            tol = 1e-5,
        )
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system)
    }

    return(lm_result)
}


#' Extract Error Correction Term from Systemfit Model
#'
#' Computes the Error Correction Term from Unrestricted ECM coefficients.
#'
#' @param systemfit_uecm_coefs A list containing coefficients from UECM.
#' @param sel_variables A character vector of selected variable names.
#'
#' @return A data.table object containing the error correction term.
#' @export
get_ect_systemfit <- function(systemfit_uecm_coefs, sel_variables, table_dt) {
    coef_exp <- systemfit_uecm_coefs$coefficients

    lags_x <- coef_exp[names(systemfit_uecm_coefs$coefficients) %like% "lag"]

    # Initialize ect_x with the first term
    ect_x <- table_dt[, get(sel_variables[1])]

    # Loop through the rest of the variables in sel_variables
    for (i in 2:length(sel_variables)) {
        term <-
            table_dt[, get(sel_variables[i])] * lags_x[names(lags_x) %like% sel_variables[i]] /
                abs(lags_x[names(lags_x) %like% sel_variables[1]])
        ect_x <- ect_x - term
    }
    transf <- data.table::data.table(ect_x)
    return(transf)
}


#' Restricted Error Correction Model for Systemfit
#'
#' Computes the Restricted Error Correction Model for Systemfit.
#'
#' @param col_names A character vector of column names.
#' @param uecm_model An object of class systemfit, representing the UECM model.
#' @param method Character string indicating the desired estimation method.
#' @param method_solv Character string indicating the solution method. Default is "EViews".
#' @param iterations An integer indicating the number of iterations.
#' @param nlags An integer specifying the number of lags.
#' @param dt A data.table object containing the data.
#' @param inst_list List of instruments for 2SLS and 3SLS.
#'
#' @return A model result from systemfit.
#' @export
recm_systemfit <- function(
    col_names = c(),
    uecm_model,
    grouping,
    method = "SUR",
    method_solv = "EViews",
    iterations = 1,
    nlags = 1,
    dt = data.table::data.table(),
    inst_list = c()) {
    diff_cols <- c()
    all_lag_cols <- c()

    # get and incorporate ECT from UECM
    ect_test <- get_ect_systemfit(
        systemfit_uecm_coefs = uecm_model,
        sel_variables = col_names,
        table_dt = dt
    )
    dt <- cbind(dt, ect_test)
    ect <- dt$ect_test
    ifelse(method != "SUR", col_names_ext <- c(col_names, inst_list[-1]), col_names_ext <- col_names)


    # Add lag columns for each lag value
    for (col in col_names_ext) {
        # Add diff column
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := diff(c(NA, get(col))), by = get(grouping)]
        diff_cols <- c(diff_cols, diff_col)
    }
    for (col in col_names) {
        if (nlags >= 2) {
            # Add lag columns for each lag value
            for (lag in 2:nlags) {
                lag_col <- paste0(col, "_lag", lag)
                dt[, (lag_col) := data.table::shift(get(col), n = lag, type = "lag"), by = get(grouping)]
                all_lag_cols <- c(all_lag_cols, lag_col)
            }
        }
    }

    if (method != "SUR") {
        diff_inst <- diff_cols[!(diff_cols %like% inst_list[1])]
        inst_eq <- paste("~", paste(c(diff_inst[-1], all_lag_cols), collapse = " + "))
        diff_cols <- diff_cols[!(diff_cols %like% inst_list[-1])]
    }

    # Construct formula string
    ifelse(nlags >= 2,
        formula_str <- paste(diff_cols[1], "~", paste(c(diff_cols[-1], ect, all_lag_cols), collapse = " + ")),
        formula_str <- paste(diff_cols[1], "~", paste(c(diff_cols[-1], ect), collapse = " + "))
    )

    # Remove rows with NA values
    dt <- dt[complete.cases(dt), ]
    # Run systemfit model
    dt <- plm::pdata.frame(dt, index = c("reporter", "year"))

    # Run systemfit model
    if (method == "3SLS") {
        control_system <- systemfit::systemfit.control(
            methodResidCov = "noDfCor",
            residCovWeighted = FALSE,
            maxiter = iterations,
            tol = 1e-5,
            method3sls = "EViews" # GLS(default), IV, GMM, SCHMIDT, EVIEWS
        )
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system, inst = as.formula(inst_eq))
    }
    if (method == "2SLS") {
        control_system <- systemfit::systemfit.control(
            methodResidCov = "noDfCor",
            residCovWeighted = FALSE,
            tol = 1e-5,
        )
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system, inst = as.formula(inst_eq))
    }
    if (method == "SUR") {
        control_system <- systemfit::systemfit.control(
            methodResidCov = "noDfCor",
            residCovWeighted = FALSE,
            maxiter = iterations,
            tol = 1e-5,
        )
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system)
    }
    return(lm_result)
}


#' Bounds F-Test for Systemfit Error Correction Model
#'
#' Applies the Bounds F-Test to the system equations based on Pesaran (2001).
#'
#' @param system_ecm An object of class systemfit, representing the ECM.
#' @param units A character vector specifying the units or entities for the model.
#' @param sel_variables Variables included in the cointegration relationship
#'
#' @return A numeric vector with F-test results for each unit.
#' @export
systemfit_boundsF_test <- function(
    system_ecm,
    units,
    sel_variables) {
    bound_interx <- c()

    for (n in seq_along(units)) {
        ##### BOUND TEST ESTIMATION
        bound_interx[n] <-
            aod::wald.test(
                b = coef(system_ecm$eq[[n]]),
                Sigma = vcov(system_ecm$eq[[n]]),
                Terms = 2:(length(sel_variables) + 1)
            )$result$chi2[1] / length(sel_variables)
    }

    return(bound_interx)
}
