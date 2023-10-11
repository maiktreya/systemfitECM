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

--------------

**uecm_systemfit**: Estimates Unrestricted-ECM

- first listed variable is used as dependant variable.
- if instruments are included, the first listed variable is used as the endogenous variable and the remaining elements their instruments.

#### Unrestricted Error Correction Model (UECM)

The UECM can be represented in a more generalized format as:

$$\Delta y_t = c_0 + c_1 t + \pi_y y_{t-1} + \sum_{j=1}^{k} \pi_j x_{j,t-1} + \sum_{i=1}^{p-1} \psi_{y,i} \Delta y_{t-i} + \sum_{j=1}^{k} \sum_{l=1}^{q_j-1} \psi_{j,l} \Delta x_{j,t-l} + \sum_{j=1}^{k} \omega_j \Delta x_{j,t} + \epsilon_t$$

Where:

- $\Delta y_t$: Change in the dependent variable at time $t$.

- $c_0$: Constant term.

- $c_1$: Trend coefficient.

- $\pi_y$: Coefficient of the lagged level of the dependent variable, capturing the long-run adjustment to equilibrium.

- $y_{t-1}$: Lagged level of the dependent variable.

- $x_{j,t-1}$: Lagged levels of the independent variables.

- $\Delta y_{t-i}, \Delta x_{j,t-l}$: First differences of the dependent and independent variables, respectively.

- $\epsilon_t$: Error term at time$t$.

Restrictions:

- $\psi_{j,l} = 0$ for all $q_j \leq 1$

- $\psi_{y,i} = 0$ if $p = 1$

- --------------

**get_ect_systemfit**: Generates a Error Correction Term serie for a given UECM. Called internally by other functions.

--------------

**recm_systemfit**: Estimates Restricted-ECM (requires a previous obj. class uecm_systemfit).

#### Restricted Error Correction Model (RECM)

The RECM can be derived from the UECM by imposing restrictions. In the RECM, the ECT (Error Correction Term) captures the long-run relationship, which gets adjusted in the short run:

$$\Delta y_t = c_0 + \rho ECT_{t-1} + \sum_{i=1}^{p-1} \theta_{y,i} \Delta y_{t-i} + \sum_{j=1}^{m} \sum_{l=1}^{q_j-1} \theta_{j,l} \Delta x_{j,t-l} + \epsilon_t$$

Where:

- $\Delta y_t$: Change in the dependent variable at time $t$.
- $c_0$: Constant term.
- $ECT_{t-1}$: Error Correction Term, which is the lagged residual from the cointegration relationship, typically represented as $y_{t-1} - \phi x_{t-1}$.
- $\rho$: Coefficient capturing the speed of adjustment back to the long-run equilibrium.
- $\Delta y_{t-i}, \Delta x_{j,t-l}$: First differences of the dependent and independent variables, respectively.
- $\epsilon_t$: Error term at time $t$.

Restrictions:

- Certain lags or coefficients might be restricted or omitted based on statistical tests or theoretical considerations.

--------------

**systemfit_boundsF_test**: Pesaran et al. (2001) F-Bounds Test (requires a previous obj. class uecm_systemfit).

- corresponding p-values should be taken from the original article or Nayaran (2003). This version only works for case III (intercept, no time trend.)

------------------------------------------------

### PENDING TO IMPLEMENT

- Check F Bound  test for ARDL distinct of (p=1,q=2) implemented 04-10-23.
- Enable beyond case III which is the only version already tested an implemented.
- Check implications of distinct short-run parameters for UECM and RECM when:
  - applied a first step iterated + uniterated or iterated.
  - to keep reporting a "pure" reparametrization iterations should be set to 1 (not iterated).

### REVS

#### 10/10/2023

- Fixed computation of ECT. For that purpose two additional parameters are set for its calculation and the calculation of associated RECM (nunits and nperiods).
- Fixed RECM parameterization (lags were removed and only multi order diffs are included along ECT).
