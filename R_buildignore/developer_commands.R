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
citation("imputeFin")
vignette(package = "imputeFin")


##
## Developer commands (https://r-pkgs.org/, https://style.tidyverse.org/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::document()  #to generate all documentation via roxygen
devtools::install()
library(imputeFin)


# Code tests
devtools::test()
#covr::package_coverage()  #coverage of tests


# CRAN check and submission (https://r-pkgs.org/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()  # run_dont_test = TRUE
rcmdcheck::rcmdcheck()  # build_args = "--run-donttest"
devtools::build()
#devtools::revdep(pkg = "imputeFin")  # to check reverse dependencies
#devtools::check_win_release()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD check imputeFin_0.1.0.tar.gz --as-cran --run-donttest  # this is before submission to CRAN
#R CMD install imputeFin_0.1.0.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html
