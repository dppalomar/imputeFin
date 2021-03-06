##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type = "source")
# Installation from GitHub
devtools::install_github("dppalomar/imputeFin")
# Installation from CRAN
install.packages("imputeFin")
# Getting help
library(imputeFin)
help(package = "imputeFin")
package?imputeFin
?impute_AR1_Gaussian
vignette("ImputeFinancialTimeSeries", package = "imputeFin")
browseVignettes("imputeFin")
RShowDoc("ImputeFinancialTimeSeries", package = "imputeFin")
citation("imputeFin")
#tools::showNonASCIIfile("R/estimate_impute_AR1_Gaussian.R")


##
## Developer commands (https://r-pkgs.org/, https://style.tidyverse.org/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::document()  #to generate all documentation via roxygen
devtools::install()
devtools::install(build_vignettes = TRUE)
library(imputeFin)


# Code tests
devtools::test()


# Reverse dependencies
devtools::revdep(pkg = "imputeFin")
revdepcheck::revdep_check(num_workers = 4)  # https://github.com/r-lib/revdepcheck


# CRAN check and submission (https://r-pkgs.org/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()  # run_dont_test = TRUE
rcmdcheck::rcmdcheck()  # build_args = "--run-donttest"
devtools::build()
#devtools::check_win_release()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD check imputeFin_0.1.2.tar.gz --as-cran --run-donttest  # this is before submission to CRAN
#R CMD install imputeFin_0.1.2.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html
