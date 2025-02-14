#' Unrestricted Error Correction Model for Systemfit
#'
#' Computes the Unrestricted Error Correction Model (UECM) for use with `systemfit`.
#'  This function estimates a single equation from a system of equations.  It's designed
#'  for panel data.
#'
#' @param formula A formula of the form `y ~ x1 + x2 | z1 + z2`, where `y` is the
#'   dependent variable, `x1` and `x2` are independent variables, and `z1` and
#'   `z2` are instrumental variables (optional).  The formula *must* be
#'   explicitly specified; using a character string for the formula is *not*
#'   supported by `systemfit`.
#' @param data A `data.table` object containing the data.
#' @param grouping A character vector of length two, specifying the columns
#'   defining panel units and time periods (e.g., `c("reporter", "year")`).
#' @param nlags An integer specifying the number of lags to include.  Defaults to 1.
#' @param method A character string indicating the desired estimation method.
#'   Must be one of "SUR" (default), "OLS", "2SLS", or "3SLS".
#' @param method_solv A character string indicating the solution method for 3SLS.
#'   Default is "EViews".  See `systemfit.control` for other options.
#' @param iterations An integer indicating the number of iterations for the
#'   estimation. Defaults to 1.
#' @param inst_list A character vector of instruments for 2SLS and 3SLS methods.
#'   The first element should be the name of the endogenous variable being instrumented.
#'   The remaining elements are the instrument names. This is crucial for correct
#'   instrument specification.  Ignored if `method` is "SUR" or "OLS".
#'
#' @return A model result from `systemfit`.
#' @import data.table
#' @importFrom systemfit systemfit systemfit.control
#' @importFrom stats diff complete.cases as.formula coef vcov
#' @importFrom plm pdata.frame
#' @export
uecm_systemfit <- function(
    formula,
    data,
    grouping,
    nlags = 1,
    method = "SUR",
    method_solv = "EViews",
    iterations = 1,
    inst_list = NULL) {

    # Input validation
    if (!inherits(formula, "formula")) {
        stop("`formula` must be a formula object.")
    }
    if (!inherits(data, "data.table")) {
        stop("`data` must be a data.table.")
    }
    if (length(grouping) != 2) {
        stop("`grouping` must be a character vector of length 2 (unit, time).")
    }
    if (!method %in% c("SUR", "OLS", "2SLS", "3SLS")) {
      stop("`method` must be one of 'SUR', 'OLS', '2SLS', or '3SLS'.")
    }
    if (method %in% c("2SLS", "3SLS") && (is.null(inst_list) || length(inst_list) < 2)) {
      stop("For 2SLS and 3SLS, `inst_list` must be provided and have at least two elements (endogenous variable and instruments).")
    }


    dt <- data.table::copy(data)  # Avoid modifying the original data.table

    # --- Construct the formula dynamically (but we receive it ready) ---
    # All variables from the formula (much safer than extracting from col_names)
    all_vars <- all.vars(formula)
    response_var <- all_vars[1]  # The dependent variable (LHS of formula)
    explanatory_vars <- all_vars[-1]  # All RHS variables


    # --- Create differenced and lagged variables ---
    diff_cols <- character(0)  # Use character(0) for empty vectors
    all_lag_cols <- character(0)

    # Differenced variables (for all explanatory variables, including instruments)
    for (col in explanatory_vars) {
        diff_col <- paste0(col, "_diff")
        dt[, (diff_col) := c(NA, diff(get(col))), by = get(grouping[1])] #Efficient differencing
        diff_cols <- c(diff_cols, diff_col)
    }

    # Lagged variables (only for original, *non-instrument* variables in col_names)
    # We derive col_names from the formula to ensure we lag the correct things.
    # Get col_names (non-instrument variables) from the formula
    if(!is.null(inst_list)){
        col_names <- explanatory_vars[!(explanatory_vars %in% inst_list[-1])]
    } else {
        col_names <- explanatory_vars
    }

    for (col in col_names) {
        for (lag in 1:nlags) {
            lag_col <- paste0(col, "_lag", lag)
            dt[, (lag_col) := data.table::shift(get(col), n = lag, type = "lag"), by = get(grouping[1])]
            all_lag_cols <- c(all_lag_cols, lag_col)
        }
    }


    # --- Prepare for systemfit ---

    dt <- dt[complete.cases(dt), ]  # systemfit requires complete cases
    dt <- plm::pdata.frame(dt, index = grouping)

    control_system <- systemfit::systemfit.control(
        methodResidCov = "noDfCor",
        residCovWeighted = FALSE,
        maxiter = iterations,
        tol = 1e-5,
        method3sls = if (method == "3SLS") method_solv else NULL
    )

    # --- Instrumental Variables Handling ---
    if (method %in% c("2SLS", "3SLS")) {
        # Build the instrument formula dynamically
        endogenous_var <- inst_list[1]
        instruments <- inst_list[-1]

        # Non-differenced instruments, and lags of all *original* variables (from col_names)
        inst_formula_rhs <- paste(c(instruments, all_lag_cols), collapse = " + ")
        inst_formula_str <- paste("~", inst_formula_rhs)
        inst_formula <- as.formula(inst_formula_str)

        lm_result <- systemfit::systemfit(
            formula,
            data = dt,
            method = method,
            control = control_system,
            inst = inst_formula
        )
    } else {
        # For SUR and OLS, no instruments are needed
        lm_result <- systemfit::systemfit(
            formula,
            data = dt,
            method = method,
            control = control_system
        )
    }

    return(lm_result)
}



#' Extract Error Correction Term from Systemfit Model
#'
#' Computes the Error Correction Term (ECT) from coefficients of an
#' Unrestricted Error Correction Model (UECM) estimated using `systemfit`.
#'
#' @param systemfit_uecm A fitted model object of class `systemfit` from
#'   `uecm_systemfit`.
#' @param sel_variables A character vector of selected variable names *in
#'   levels* (not differenced) that form the cointegrating relationship. The
#'   first variable *must* be the dependent variable in levels.
#' @param data A `data.table` object containing the original data (must include
#'    the `grouping` columns).
#' @param grouping A character vector of length two, specifying the columns
#'   defining panel units and time periods (e.g., `c("reporter", "year")`).
#'
#' @return A `data.table` object containing the error correction term.
#' @import data.table
#' @importFrom stats coef
#' @export
get_ect_systemfit <- function(
    systemfit_uecm,
    sel_variables,
    data,
    grouping) {

    # Input validation
    if (!inherits(systemfit_uecm, "systemfit")) {
        stop("`systemfit_uecm` must be a systemfit object.")
    }
    if (!inherits(data, "data.table")) {
        stop("`data` must be a data.table.")
    }
    if(length(sel_variables) < 2){
        stop("`sel_variables` needs at least the dependent and one independent variable.")
    }
    if (length(grouping) != 2) {
        stop("`grouping` must be a character vector of length 2 (unit, time).")
    }


    table_dt <- data.table::copy(data)
    coef_exp <- coef(systemfit_uecm)

    # Get the names of the lagged variables from the coefficients
    lag_names <- grep("_lag[0-9]+$", names(coef_exp), value = TRUE)

    # Extract the base variable names from the lagged variable names
    base_var_names <- sub("_lag[0-9]+$", "", lag_names)

    # Ensure that all required sel_variables have corresponding lags
    if (!all(sel_variables %in% unique(base_var_names))) {
       stop("Not all `sel_variables` have corresponding lagged coefficients in the `systemfit_uecm` model.")
    }

    # --- Calculate the ECT ---
    # Initialize with the dependent variable (first sel_variable)
    ect_x <- table_dt[, get(sel_variables[1])]
    dep_var_lag_coef <- coef_exp[grep(paste0("^", sel_variables[1], "_lag[0-9]+$"), names(coef_exp))]

     if (length(dep_var_lag_coef) == 0) {
        stop(paste("No lag coefficient found for the dependent variable:", sel_variables[1]))
    }


    # Loop through the remaining sel_variables (independent variables)
    for (i in 2:length(sel_variables)) {
        var_name <- sel_variables[i]
        lag_coef_name <- grep(paste0("^",var_name, "_lag[0-9]+$"), names(coef_exp), value = TRUE)

        if(length(lag_coef_name) > 1) {
          warning(paste("Multiple lag coefficients found for", var_name,
                        ". Using the first one:", lag_coef_name[1]))
          lag_coef_name <- lag_coef_name[1] # Take only the first if there are multiples.
        }

        if(length(lag_coef_name) == 0){
            stop(paste("No lag coefficient found for:", var_name))
        }

        lag_coef <- coef_exp[lag_coef_name]
        term <- table_dt[, get(var_name)] * lag_coef / abs(dep_var_lag_coef[1])
        ect_x <- ect_x - term
    }

    # Create the ECT data.table, handling the initial NA correctly
    transf <- data.table::data.table(
        ect_x = ect_x,
        key = table_dt[, get(grouping[1])],
        time = table_dt[, get(grouping[2])]
    )

    #Shift time to account for differencing and lagging
    transf[, time := time - (nlags +1) ]
    transf <- transf[time >= min(data[,get(grouping[2])])] # Remove any leading NA rows


    return(transf)
}




#' Restricted Error Correction Model for Systemfit
#'
#' Computes the Restricted Error Correction Model (RECM) using a
#' precomputed UECM model and an extracted Error Correction Term (ECT).
#'
#' @param formula A formula of the form `y ~ x1 + x2 | z1 + z2`, where `y` is the
#'   dependent variable, `x1` and `x2` are independent variables, and `z1` and
#'   `z2` are instrumental variables (optional).  The formula *must* be
#'   explicitly specified.
#' @param data A `data.table` containing the original data.
#' @param grouping A character vector of length two specifying the panel unit
#'    and time period columns (e.g., `c("reporter", "year")`).
#' @param uecm_model A fitted model object of class `systemfit` representing the
#'   precomputed UECM model. This is used to calculate the ECT.
#' @param nlags An integer specifying the number of lags included in the UECM.
#'   Defaults to 1. This should match the `nlags` used in `uecm_systemfit`.
#' @param method A character string indicating the desired estimation method.
#'   Must be one of "SUR" (default), "OLS", "2SLS", or "3SLS".
#' @param method_solv A character string indicating the solution method for 3SLS.
#'   Default is "EViews".
#' @param iterations An integer indicating the number of iterations. Defaults to 1.
#' @param inst_list A character vector of instruments for 2SLS and 3SLS methods.  See
#'    `uecm_systemfit` for details.
#'
#' @return A model result from `systemfit`.
#' @import data.table
#' @importFrom systemfit systemfit systemfit.control
#' @importFrom stats diff complete.cases as.formula
#' @importFrom plm pdata.frame
#' @export
recm_systemfit <- function(
    formula,
    data,
    grouping,
    uecm_model,
    nlags = 1,
    method = "SUR",
    method_solv = "EViews",
    iterations = 1,
    inst_list = NULL) {

    # Input validation (similar to uecm_systemfit)
    if (!inherits(formula, "formula")) {
        stop("`formula` must be a formula object.")
    }
    if (!inherits(data, "data.table")) {
        stop("`data` must be a data.table.")
    }
    if (!inherits(uecm_model, "systemfit")) {
      stop("`uecm_model` must be a systemfit object.")
    }
    if (length(grouping) != 2) {
        stop("`grouping` must be a character vector of length 2 (unit, time).")
    }
    if (!method %in% c("SUR", "OLS", "2SLS", "3SLS")) {
        stop("`method` must be one of 'SUR', 'OLS', '2SLS', or '3SLS'.")
    }
    if (method %in% c("2SLS", "3SLS") && (is.null(inst_list) || length(inst_list) < 2)) {
      stop("For 2SLS and 3SLS, `inst_list` must be provided and have at least two elements.")
    }

    dt <- data.table::copy(data)

    # --- Prepare variables from formula ---
    all_vars <- all.vars(formula)
    response_var <- all_vars[1]
    explanatory_vars <- all_vars[-1]

    # Get col_names (non-instrument variables) from the formula
    if(!is.null(inst_list)){
        col_names <- explanatory_vars[!(explanatory_vars %in% inst_list[-1])]
    } else {
        col_names <- explanatory_vars
    }


    # --- Calculate ECT ---
    ect_test <- get_ect_systemfit(
        systemfit_uecm = uecm_model,
        sel_variables = c(response_var, col_names),  # Include response variable!
        data = dt,
        grouping = grouping
    )
    dt[, ect := ect_test$ect_x]

    # --- Create differenced and lagged variables ---
    diff_cols <- character(0)
    all_lag_cols <- character(0) # For differences of lags

    # Differenced variables (all variables in the formula)
    for (col in c(response_var, explanatory_vars)) { # Include response
        diff_col <- paste0(col, "_diff")
         dt[, (diff_col) := c(NA, diff(get(col))), by = get(grouping[1])]
        diff_cols <- c(diff_cols, diff_col)
    }

    # Lagged *differenced* variables (only for original, non-instrument variables)
    for (col in col_names) {
        for (lag in 1:nlags) { # Start lag at 1, not 2
            lag_col <- paste0(col, "_diff_lag", lag)
            dt[, (lag_col) := data.table::shift(get(paste0(col,"_diff")), n = lag, type = "lag"), by = get(grouping[1])]
            all_lag_cols <- c(all_lag_cols, lag_col)
        }
    }
    # --- Prepare for systemfit ---
    dt <- dt[complete.cases(dt), ]
    dt <- plm::pdata.frame(dt, index = grouping)

    control_system <- systemfit::systemfit.control(
        methodResidCov = "noDfCor",
        residCovWeighted = FALSE,
        maxiter = iterations,
        tol = 1e-5,
        method3sls = if (method == "3SLS") method_solv else NULL
    )

    # --- Construct the RECM formula ---
    # We already have diff_cols and all_lag_cols
    recm_formula_rhs <- paste(c(diff_cols[-1], all_lag_cols, "ect"), collapse = " + ")
    recm_formula_str <- paste(diff_cols[1], "~", recm_formula_rhs)
    recm_formula <- as.formula(recm_formula_str)

    # --- Instrumental Variables (RECM) ---
    if (method %in% c("2SLS", "3SLS")) {
        endogenous_var <- inst_list[1]
        instruments <- inst_list[-1]

        # Build instrument formula, including the ECT
        inst_formula_rhs <- paste(c(instruments, all_lag_cols, "ect"), collapse = " + ")
        inst_formula_str <- paste("~", inst_formula_rhs)
        inst_formula <- as.formula(inst_formula_str)

        lm_result <- systemfit::systemfit(
            recm_formula,  # Use the dynamically created formula
            data = dt,
            method = method,
            control = control_system,
            inst = inst_formula
        )
    } else {
        # SUR and OLS
        lm_result <- systemfit::systemfit(
            recm_formula,  # Use the dynamically created formula
            data = dt,
            method = method,
            control = control_system
        )
    }

    return(lm_result)
}



#' Bounds F-Test for Systemfit Error Correction Model
#'
#' Applies the Bounds F-Test (Pesaran, Shin, and Smith, 2001) to a single equation of the
#' system equations estimated by `systemfit`. This tests for cointegration.
#'
#' @param system_ecm An object of class `systemfit`, representing the Error
#'   Correction Model (ECM), either UECM or RECM.  The test is applied to the
#'   first equation in the system.
#' @param sel_variables A character vector of variables *in levels* included in
#'   the cointegration relationship. The first variable must correspond to the
#'   dependent variable in the ECM.
#'
#' @return A single numeric value representing the F-statistic.
#' @importFrom aod wald.test
#' @importFrom stats coef vcov
#' @export
systemfit_boundsF_test <- function(
    system_ecm,
    sel_variables) {

    # Input validation
    if (!inherits(system_ecm, "systemfit")) {
        stop("`system_ecm` must be a systemfit object.")
    }
    if (length(sel_variables) < 2) {
      stop("`sel_variables` must contain at least the dependent and one independent variable.")
    }

    # Extract coefficients and variance-covariance matrix for the first equation
    eq1_coefs <- coef(system_ecm$eq[[1]])
    eq1_vcov <- vcov(system_ecm$eq[[1]])

    # Identify the indices of the lagged *level* variables in the coefficient vector
    lag_level_indices <- integer(0) # Initialize empty vector
    for(var in sel_variables){
        indices <- grep(paste0("^", var, "_lag[0-9]+$"), names(eq1_coefs))
        lag_level_indices <- c(lag_level_indices, indices)
    }

    if (length(lag_level_indices) == 0) {
      stop("No lagged level variables found matching `sel_variables` in the model.")
    }
    # Ensure that the selection is sorted, since combining indices via grep might
    #  return unsorted indices.
    lag_level_indices <- sort(lag_level_indices)


    # Perform the Wald test
    wald_result <- aod::wald.test(
        b = eq1_coefs,
        Sigma = eq1_vcov,
        Terms = lag_level_indices
    )

    # Calculate and return the F-statistic
    f_statistic <- wald_result$result$chi2[1] / length(sel_variables)
    return(f_statistic)
}

