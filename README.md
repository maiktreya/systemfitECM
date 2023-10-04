# Package: systemfitECM

    - Type: Package
    - Title: Cointegration Techniques into systemfit
    - Version: 0.1.0
    - Author: Miguel Garcia Duch
    - Description: Enabling estimating restricted and unrestricted Error Correction Models (ECM) and performing Pesaran(2001) Bounds-F test for cointegration within systemfit objects/models.
    - License: GPL-3
    - Imports:
library(data.table) # for simple and performant data manipulation 
library(plm) # needed for systemfit to handle panel structure 
library(systemfit # for FGLS system linear models 
library(magrittr) # For piping with %<% without dplyr dependencies 
library(aod) # for performing F Bounds test
source("R/functions.R")

### Package functions and example file inside R folder
