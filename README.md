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

### IMPLEMENTATION

Pacakage functions:

- **uecm_systemfit**: Estimates Unrestricted-ECM
- **get_ect_systemfit**: Generates a Error Correction Term serie for a given UECM. Called internally by other functions.
- **recm_systemfit**: Estimates Restricted-ECM (requires a previous obj. class uecm_systemfit).
- **systemfit_boundsF_test**: Pesaran et al. (2001) F-Bounds Test (requires a previous obj. class uecm_systemfit).

------------------------------------------------

## PENDING TO IMPLEMENT

-> Check F Bound  test for ARDL(p=1,q=2) which is the only version already tested an implemented.
