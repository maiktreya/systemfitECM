# pipeline for upgrading to last version in local dev

rm(list = ls())
gc(full = TRUE)

# devtools::document()
devtools::load_all(reset = TRUE)
devtools::install(dependencies = TRUE) # To ensure a clean install
devtools::build()

# install.packages("dist/systemfitECM_0.1.0.tar.gz",repos = NULL, type="source")
