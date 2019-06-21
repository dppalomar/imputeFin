##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/imputeFin")
# Installation from CRAN
install.packages("imputeFin")
# Getting help
library(imputeFin)
help(package = "imputeFin")
package?imputeFin
?imputeAR1Gaussian
citation("imputeFin")
vignette(package = "imputeFin")


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
library(imputeFin)

# Documentation
devtools::document()  #to generate all documentation via roxygen
?imputeAR1Gaussian


# Code tests
devtools::test()
#covr::package_coverage()  #coverage of tests


# CRAN check and submission (http://r-pkgs.had.co.nz/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()
rcmdcheck::rcmdcheck()
devtools::build()
#devtools::revdep(pkg = "imputeFin")  # to check reverse dependencies
#devtools::check_win_release()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD check imputeFin_0.1.0.tar.gz --as-cran # this is before submission to CRAN
#R CMD install imputeFin_0.1.0.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html
