# Package: systemfitECM

#### 04-10-2023

Description: Enabling the estimation of restricted and unrestricted Error Correction Models (ECM) and performing Pesaran(2001) Bounds-F test for cointegration within systemfit objects/models.

    - Type: Package
    - Title: Cointegration Techniques into systemfit
    - Version: 0.1.0
    - Author: Miguel Garcia Duch
    - License: GPL-3
    - Imports:
    - library(data.table) # for simple and performant data manipulation
    - library(plm) # needed for systemfit to handle panel structure
    - library(systemfit) # for FGLS system linear models
    - library(magrittr) # For piping with %<% without dplyr dependencies
    - library(aod) # for performing F Bounds test

-----------------------------------------------

## Installation

Install the latest version of this package with:

``` r
devtools::install_github("iliciuv/systemfitECM")
```

-----------------------------------------------

## Example usage

Run:

``` r
source('tests/examples.R')
```

Or run the code included in the script for a working application on simulated data

### Core functionality

Pacakage functions:

- **uecm_systemfit**: Estimates Unrestricted-ECM
  - first listed variable is used as dependant variable.
  - if instruments are included, the first listed variable is used as the endogenous variable and the remaining elements their instruments.
- **get_ect_systemfit**: Generates a Error Correction Term serie for a given UECM. Called internally by other functions.
- **recm_systemfit**: Estimates Restricted-ECM (requires a previous obj. class uecm_systemfit).
- **systemfit_boundsF_test**: Pesaran et al. (2001) F-Bounds Test (requires a previous obj. class uecm_systemfit).
  - corresponding p-values should be taken from the original article or Nayaran (2003). This version only works for case III (intercept, no time trend.)

------------------------------------------------

## PENDING TO IMPLEMENT

-> Check F Bound  test for ARDL distinct of (p=1,q=2) case III which is the only version already tested an implemented.

``` r
systemfit_boundsF_test <- function(
    system_ecm,
    units) {
    bound_interx <- c()
    for (n in seq_along(units)) {
        ##### BOUND TEST ESTIMATION
        bound_interx[n] <- aod::wald.test(b = coef(system_ecm$eq[[n]]), Sigma = vcov(system_ecm$eq[[n]]), Terms = 2:4)$result$chi2[1] / 3
    }

    return(bound_interx)
}
# possible refraction of code:
systemfit_boundsF_test <- function(
    system_ecm,
    units) {
    bound_interx <- c()
    for (n in seq_along(units)) {
        ##### BOUND TEST ESTIMATION
        bound_interx[n] <- aod::wald.test(b = coef(system_ecm$eq[[n]]), Sigma = vcov(system_ecm$eq[[n]]), Terms = 2:length(sel_variables)+ 1)$result$chi2[1] / length(sel_variables)
    }

    return(bound_interx)
}

```
