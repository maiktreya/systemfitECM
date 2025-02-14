#### MODULE FOR INTEGRATION OF ECM INTO SYSTEMFIT

#' Estimate an Unrestricted Error Correction Model (UECM) using systemfit
#'
#' Estimates a UECM for panel data using \code{systemfit}.
#'
#' @param data A data.frame or data.table.
#' @param col_names Character vector of variable names.
#' @param nlags Integer, number of lags.
#' @param grouping Character string, name of the grouping column.
#' @param method "SUR", "2SLS", or "3SLS".
#' @param method_solv Solution method for 3SLS (see \code{systemfit.control}).
#' @param iterations Integer, max iterations.
#' @param inst_list Character vector of instruments (for 2SLS/3SLS).
#'
#' @return A \code{systemfit} object.
#' @examples
#' # Example (replace with your data)
#' # data(your_data)
#' # uecm_model <- uecm_systemfit(data = your_data,...)
#' @importFrom data.table copy shift
#' @importFrom plm pdata.frame
#' @importFrom systemfit systemfit systemfit.control
#' @export
uecm_systemfit <- function(data, col_names, nlags = 1, grouping, method = "SUR",
                           method_solv = "EViews", iterations = 1, inst_list = NULL) {
    dt <- data.table::copy(data)
    stopifnot(is.data.frame(dt) || is.data.table(dt), "data must be a data.frame or data.table")
    stopifnot(is.character(col_names), "col_names must be a character vector")
    stopifnot(is.numeric(nlags), "nlags must be numeric")
    stopifnot(nlags >= 0, "nlags must be non-negative")
    stopifnot(is.character(grouping), "grouping must be a character string")
    stopifnot(grouping %in% names(dt), paste("grouping variable", grouping, "not found in data"))
    stopifnot(method %in% c("SUR", "2SLS", "3SLS"), "method must be one of 'SUR', '2SLS', '3SLS'")

    if (method %in% c("2SLS", "3SLS")) {
        stopifnot(!is.null(inst_list), "inst_list must be provided for 2SLS/3SLS")
        stopifnot(is.character(inst_list), "inst_list must be a character vector")
        stopifnot(inst_list %in% col_names, "First element of inst_list must be an endogenous variable")
    }

    diff_cols <- character()
    all_lag_cols <- character()
    all_inst_cols <- character() # for instruments

    for (col in col_names) {
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := diff(c(NA, get(col))), by = get(grouping)]
        diff_cols <- c(diff_cols, diff_col)
        for (lag in 1:nlags) {
            lag_col <- paste0(col, "_lag", lag)
            dt[, (lag_col) := data.table::shift(get(col), n = lag, type = "lag"), by = get(grouping)]
            all_lag_cols <- c(all_lag_cols, lag_col)
        }
    }

    if (method %in% c("2SLS", "3SLS")) {
        for (inst in inst_list[-1]) { # Instruments (excluding endogenous variable)
            for (lag in 1:nlags) {
                inst_col <- paste0(inst, "_lag", lag)
                dt[, (inst_col) := data.table::shift(get(inst), n = lag, type = "lag"), by = get(grouping)]
                all_inst_cols <- c(all_inst_cols, inst_col)
            }
        }
    }

    formula_str <- paste(diff_cols, "~", paste(c(diff_cols[-1], all_lag_cols), collapse = " + "))

    if (method %in% c("2SLS", "3SLS")) {
        inst_eq <- paste("~", paste(c(all_lag_cols, all_inst_cols), collapse = " + ")) # Instruments
    }

    dt <- dt[complete.cases(dt), ]
    dt <- plm::pdata.frame(dt, index = unique(grouping)) # Use unique to avoid error if grouping has many names
    control_system <- systemfit::systemfit.control(
        methodResidCov = "noDfCor", residCovWeighted = FALSE, maxiter = iterations, tol = 1e-5
    )

    if (method == "3SLS") {
        control_system$method3sls <- method_solv
        lm_result <- systemfit::systemfit(as.formula(formula_str),
            data = dt, method = method,
            control = control_system, inst = as.formula(inst_eq)
        )
    } else if (method == "2SLS") {
        lm_result <- systemfit::systemfit(as.formula(formula_str),
            data = dt, method = method,
            control = control_system, inst = as.formula(inst_eq)
        )
    } else {
        lm_result <- systemfit::systemfit(as.formula(formula_str),
            data = dt, method = method,
            control = control_system
        )
    }

    return(lm_result)
}


#' Extract Error Correction Term (ECT)
#'
#' Computes the ECT from UECM coefficients.
#'
#' @param systemfit_uecm_coefs \code{systemfit} object (UECM results).
#' @param sel_variables Character vector of variable names.
#' @param table_dt Data.table containing the original data.
#' @param nperiods Number of time periods.
#' @param nunits Number of panel units.
#'
#' @return A data.table with the ECT.
#' @export
get_ect_systemfit <- function(systemfit_uecm_coefs, sel_variables, table_dt, nperiods, nunits) {
    table_dt <- data.table::copy(table_dt)
    coef_exp <- systemfit_uecm_coefs$coefficients
    lags_x <- coef_exp[names(systemfit_uecm_coefs$coefficients) %like% "lag"]

    # Check if coefficients exist for all selected variables
    missing_coefs <- setdiff(sel_variables, gsub("_lag+", "", names(lags_x)))
    if (length(missing_coefs) > 0) {
        stop(paste("Missing coefficients for variables:", paste(missing_coefs, collapse = ", ")))
    }

    time <- rep(c(1:(nperiods + 1)), nunits)
    key <- c()

    # Initialize ect_x with the first term
    ect_x <- table_dt[, get(sel_variables)]

    # Loop through the rest of the variables in sel_variables
    for (i in 2:length(sel_variables)) {
        term <-
            table_dt[, get(sel_variables[i])] * lags_x[names(lags_x) %like% sel_variables[i]] /
                abs(lags_x[names(lags_x) %like% sel_variables])
        ect_x <- ect_x - term
    }

    for (i in 1:nunits) key <- c(key, rep(i, nperiods))
    transf <- data.table::data.table(ect_x)
    transf <- cbind(key, time, transf)
    transf[time == nperiods + 1, ect_x := NA]
    transf[time == nperiods + 1, time := 1]
    transf <- transf[order(key, time)]
    return(transf)
}

#' @param grouping Character string, name of the grouping column.
#' @param method "SUR", "2SLS", or "3SLS".
#' @param method_solv Solution method for 3SLS.
#' @param iterations Integer, max iterations.
#' @param nunits Number of panel units.
#' @param nperiods Number of time periods.
#' @param nlags Number of lags.
#' @param inst_list Character vector of instruments (for 2SLS/3SLS).
#'
#' @return A \code{systemfit} object.
#' @export
recm_systemfit <- function(data, col_names, uecm_model, grouping, method = "SUR",
                           method_solv = "EViews", iterations = 1, nunits, nperiods,
                           nlags = 1, inst_list = NULL) {
    dt <- data.table::copy(data)
    stopifnot(is.data.frame(dt) || is.data.table(dt), "data must be a data.frame or data.table")
    stopifnot(is.character(col_names), "col_names must be a character vector")
    stopifnot(is.numeric(nlags), "nlags must be numeric")
    stopifnot(nlags >= 0, "nlags must be non-negative")
    stopifnot(is.character(grouping), "grouping must be a character string")
    stopifnot(grouping %in% names(dt), paste("grouping variable", grouping, "not found in data"))
    stopifnot(method %in% c("SUR", "2SLS", "3SLS"), "method must be one of 'SUR', '2SLS', '3SLS'")

    if (method %in% c("2SLS", "3SLS")) {
        stopifnot(!is.null(inst_list), "inst_list must be provided for 2SLS/3SLS")
        stopifnot(is.character(inst_list), "inst_list must be a character vector")
        stopifnot(inst_list %in% col_names, "Instruments must be in col_names")
    }

    # Get and incorporate ECT from UECM
    ect_test <- get_ect_systemfit(
        systemfit_uecm_coefs = uecm_model,
        sel_variables = col_names,
        table_dt = dt,
        nperiods = nperiods,
        nunits = nunits
    )
    dt$ect <- ect_test$ect_x

    diff_cols <- character()
    all_lag_cols <- character()
    all_inst_cols <- character()

    for (col in col_names) {
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := diff(c(NA, get(col))), by = get(grouping)]
        diff_cols <- c(diff_cols, diff_col)

        for (lag in 2:nlags) { # Lags from 2 onwards for RECM
            lag_col <- paste0(col, "_diff", lag)
            dt[, (lag_col) := diff(c(rep(NA, lag), get(col)), differences = lag), by = get(grouping)]
            all_lag_cols <- c(all_lag_cols, lag_col)
        }
    }

    if (method %in% c("2SLS", "3SLS")) {
        for (inst in inst_list[-1]) { # Instruments (excluding endogenous variable)
            for (lag in 1:nlags) {
                inst_col <- paste0(inst, "_lag", lag)
                dt[, (inst_col) := data.table::shift(get(inst), n = lag, type = "lag"), by = get(grouping)]
                all_inst_cols <- c(all_inst_cols, inst_col)
            }
        }
    }

    formula_str <- paste(diff_cols, "~", paste(c(diff_cols[-1], all_lag_cols, "ect"), collapse = " + "))

    if (method %in% c("2SLS", "3SLS")) {
        inst_eq <- paste("~", paste(c(all_lag_cols, all_inst_cols, "ect"), collapse = " + ")) # Instruments including ect
    }

    dt <- dt[complete.cases(dt), ]
    dt <- plm::pdata.frame(dt, index = unique(grouping))
    control_system <- systemfit::systemfit.control(
        methodResidCov = "noDfCor", residCovWeighted = FALSE, maxiter = iterations, tol = 1e-5
    )

    if (method == "3SLS") {
        control_system$method3sls <- method_solv
        lm_result <- systemfit::systemfit(as.formula(formula_str),
            data = dt, method = method,
            control = control_system, inst = as.formula(inst_eq)
        )
    } else if (method == "2SLS") {
        lm_result <- systemfit::systemfit(as.formula(formula_str),
            data = dt, method = method,
            control = control_system, inst = as.formula(inst_eq)
        )
    } else {
        lm_result <- systemfit::systemfit(as.formula(formula_str),
            data = dt, method = method,
            control = control_system
        )
    }

    return(lm_result)
}


#' Perform Bounds F-Test
#'
#' Performs the Bounds F-test.
#'
#' @param system_ecm \code{systemfit} object (ECM results).
#' @param units Character vector of unit names.
#' @param sel_variables Character vector of variables in cointegration relationship.
#'
#' @return Numeric vector of F-test results.
#' @export
systemfit_boundsF_test <- function(system_ecm, units, sel_variables) {
    bound_interx <- c()

    for (n in seq_along(units)) {
        ##### Bound test computation applying F form of Wald Test
        bound_interx[n] <-
            aod::wald.test(
                b = coef(system_ecm$eq[[n]]),
                Sigma = vcov(system_ecm$eq[[n]]),
                Terms = 2:(length(sel_variables) + 1) # Terms related to the variables in the cointegration relationship
            )$result$chi2 / length(sel_variables)
    }

    return(bound_interx)
}
