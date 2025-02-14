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
#' @export
uecm_systemfit <- function(data, col_names, nlags = 1, grouping, method = "SUR",
                           method_solv = "EViews", iterations = 1, inst_list = NULL) {
    dt <- data.table::copy(data)

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
        for (inst in inst_list[-1]) {
            for (lag in 1:nlags) {
                inst_col <- paste0(inst, "_lag", lag)
                dt[, (inst_col) := data.table::shift(get(inst), n = lag, type = "lag"), by = get(grouping)]
                all_inst_cols <- c(all_inst_cols, inst_col)
            }
        }
    }

    formula_str <- paste(diff_cols, "~", paste(c(diff_cols[-1], all_lag_cols), collapse = " + "))

    if (method %in% c("2SLS", "3SLS")) {
        inst_eq <- paste("~", paste(c(all_lag_cols, all_inst_cols), collapse = " + "))
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
