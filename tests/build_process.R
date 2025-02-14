# pipeline for upgrading to last version in local dev

devtools::document()
devtools::load_all()
devtools::build()

install.packages("dist/systemfitECM_0.1.0.tar.gz",repos = NULL, type="source")
