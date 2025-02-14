#' Unrestricted Error Correction Model for Systemfit
#'
#' Computes the Unrestricted Error Correction Model (UECM) for use with `systemfit`.
#'
#' @param col_names A character vector of column names to include in the model.
#' @param nlags An integer specifying the number of lags to include.
#' @param grouping A character vector specifying the column defining panel units.
#' @param method A character string indicating the desired estimation method (e.g., "SUR", "2SLS", "3SLS").
#' @param method_solv A character string indicating the solution method for 3SLS. Default is "EViews".
#' @param iterations An integer indicating the number of iterations for the estimation.
#' @param dt A `data.table` object containing the data.
#' @param inst_list A character vector of instruments for 2SLS and 3SLS methods.
#'
#' @return A model result from `systemfit`.
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

    dt <- copy(dt)
    diff_cols <- c()
    all_lag_cols <- c()

    col_names_ext <- if (method != "SUR") c(col_names, inst_list[-1]) else col_names

    for (col in col_names_ext) {
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := diff(c(NA, get(col))), by = get(grouping)]
        diff_cols <- c(diff_cols, diff_col)
    }

    for (col in col_names) {
        for (lag in 1:nlags) {
            lag_col <- paste0(col, "_lag", lag)
            dt[, (lag_col) := data.table::shift(get(col), n = lag, type = "lag"), by = get(grouping)]
            all_lag_cols <- c(all_lag_cols, lag_col)
        }
    }

    formula_str <- paste(diff_cols[1], "~", paste(c(diff_cols[-1], all_lag_cols), collapse = " + "))

    dt <- dt[complete.cases(dt), ]
    dt <- plm::pdata.frame(dt, index = c("reporter", "year"))

    control_system <- systemfit::systemfit.control(
        methodResidCov = "noDfCor",
        residCovWeighted = FALSE,
        maxiter = iterations,
        tol = 1e-5,
        method3sls = if (method == "3SLS") method_solv else NULL
    )

    if (method %in% c("2SLS", "3SLS")) {
        diff_inst <- diff_cols[!(diff_cols %like% inst_list[1])]
        inst_eq <- paste("~", paste(c(diff_inst[-1], all_lag_cols), collapse = " + "))
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system, inst = as.formula(inst_eq))
    } else {
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system)
    }

    return(lm_result)
}

#' Extract Error Correction Term from Systemfit Model
#'
#' Computes the Error Correction Term (ECT) from coefficients of an Unrestricted Error Correction Model (UECM).
#'
#' @param systemfit_uecm_coefs A list containing coefficients from a UECM.
#' @param sel_variables A character vector of selected variable names.
#' @param table_dt A `data.table` object containing the original data.
#' @param nperiods An integer specifying the number of time observations per unit.
#' @param nunits An integer reflecting the number of panel units.
#'
#' @return A `data.table` object containing the error correction term.
#' @export
get_ect_systemfit <- function(
    systemfit_uecm_coefs,
    sel_variables,
    table_dt,
    nperiods,
    nunits) {

    table_dt <- copy(table_dt)
    coef_exp <- systemfit_uecm_coefs$coefficients
    lags_x <- coef_exp[grepl("lag", names(coef_exp))]
    time <- rep(1:(nperiods + 1), nunits)
    key <- rep(1:nunits, each = nperiods + 1)

    ect_x <- table_dt[, get(sel_variables[1])]

    for (i in 2:length(sel_variables)) {
        term <- table_dt[, get(sel_variables[i])] * lags_x[grepl(sel_variables[i], names(lags_x))] / abs(lags_x[grepl(sel_variables[1], names(lags_x))])
        ect_x <- ect_x - term
    }

    transf <- data.table::data.table(ect_x = ect_x, key = key, time = time)
    transf[time == nperiods + 1, ect_x := NA]
    transf[time == nperiods + 1, time := 1]
    transf <- transf[order(key, time)]

    return(transf)
}

#' Restricted Error Correction Model for Systemfit
#'
#' Computes the Restricted Error Correction Model (RECM) for use with `systemfit` using a precomputed UECM model.
#'
#' @param col_names A character vector of column names to include in the model.
#' @param uecm_model An object of class `systemfit`, representing the precomputed UECM model.
#' @param grouping A character vector specifying the column defining panel units.
#' @param method A character string indicating the desired estimation method (e.g., "SUR", "2SLS", "3SLS").
#' @param method_solv A character string indicating the solution method for 3SLS. Default is "EViews".
#' @param iterations An integer indicating the number of iterations for the estimation.
#' @param nunits An integer reflecting the number of panel units.
#' @param nperiods An integer specifying the number of time observations per unit.
#' @param nlags An integer specifying the number of lags to include.
#' @param dt A `data.table` object containing the data.
#' @param inst_list A character vector of instruments for 2SLS and 3SLS methods.
#'
#' @return A model result from `systemfit`.
#' @export
recm_systemfit <- function(
    col_names = c(),
    uecm_model,
    grouping,
    method = "SUR",
    method_solv = "EViews",
    iterations = 1,
    nunits = 1,
    nperiods = 1,
    nlags = 1,
    dt = data.table::data.table(),
    inst_list = c()) {

    dt <- copy(dt)

    # Extract ECT from the precomputed UECM model
    ect_test <- get_ect_systemfit(
        systemfit_uecm_coefs = uecm_model,
        sel_variables = col_names,
        table_dt = dt,
        nperiods = nperiods,
        nunits = nunits
    )
    dt`${ect <- ect_test}$`ect_x

    col_names_ext <- if (method != "SUR") c(col_names, inst_list[-1]) else col_names

    diff_cols <- c()
    all_lag_cols <- c()

    for (col in col_names_ext) {
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := diff(c(NA, get(col))), by = get(grouping)]
        diff_cols <- c(diff_cols, diff_col)
    }

    for (col in col_names) {
        for (lag in 2:nlags) {
            lag_col <- paste0(col, "_diff", lag)
            dt[, (lag_col) := diff(c(rep(NA, lag), get(col)), differences = lag), by = get(grouping)]
            all_lag_cols <- c(all_lag_cols, lag_col)
        }
    }

    formula_str <- paste(diff_cols[1], "~", paste(c(diff_cols[-1], all_lag_cols, "ect"), collapse = " + "))

    dt <- dt[complete.cases(dt), ]
    dt <- plm::pdata.frame(dt, index = c("reporter", "year"))

    control_system <- systemfit::systemfit.control(
        methodResidCov = "noDfCor",
        residCovWeighted = FALSE,
        maxiter = iterations,
        tol = 1e-5,
        method3sls = if (method == "3SLS") method_solv else NULL
    )

    if (method %in% c("2SLS", "3SLS")) {
        diff_inst <- diff_cols[!(diff_cols %like% inst_list[1])]
        inst_eq <- paste("~", paste(c(diff_inst[-1], all_lag_cols, "ect"), collapse = " + "))
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system, inst = as.formula(inst_eq))
    } else {
        lm_result <- systemfit::systemfit(as.formula(formula_str), data = dt, method = method, control = control_system)
    }

    return(lm_result)
}

#' Bounds F-Test for Systemfit Error Correction Model
#'
#' Applies the Bounds F-Test to the system equations based on Pesaran (2001).
#'
#' @param system_ecm An object of class `systemfit`, representing the Error Correction Model (ECM).
#' @param units A character vector specifying the units or entities for the model.
#' @param sel_variables A character vector of variables included in the cointegration relationship.
#'
#' @return A numeric vector with F-test results for each unit.
#' @export
systemfit_boundsF_test <- function(
    system_ecm,
    units,
    sel_variables) {

    bound_interx <- sapply(seq_along(units), function(n) {
        aod::wald.test(
            b = coef(system_ecm$eq[[n]]),
            Sigma = vcov(system_ecm$eq[[n]]),
            Terms = 2:(length(sel_variables) + 1)
        )`${result}$`chi2[1] / length(sel_variables)
    })

    return(bound_interx)
}
