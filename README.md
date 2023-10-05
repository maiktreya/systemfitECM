# Package: systemfitECM

#### 04-10-2023

Description: Enabling the estimation of restricted and unrestricted Error Correction Models (ECM) and performing Pesaran(2001) Bounds-F test for cointegration within systemfit objects/models.

    - Type: Package
    - Title: Cointegration Techniques into systemfit
    - Version: 0.1.0
    - Author: Miguel Garcia-Duch
    - License: GPL-3
    - Imports:
    - library(data.table) # for simple and performant data manipulation
    - library(plm) # needed for systemfit to handle panel structure
    - library(systemfit) # for FGLS system linear models
    - library(magrittr) # For piping with %>% without dplyr dependencies
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

- Check F Bound  test for ARDL distinct of (p=1,q=2) implemented 04-10-23.
- Enable beyond case III which is the only version already tested an implemented.
- Check implications of distinct short-run parameters for UECM and RECM when:
  - applied a first step iterated + uniterated or iterated.
  - to keep reporting a "pure" reparametrization iterations should be set to 1 (not iterated).
