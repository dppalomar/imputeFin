##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/imputeFin")
# Getting help
library(imputeFin)
help(package = "imputeFin")
package?imputeFin
?XXXX


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
#devtools::install(build_vignettes = TRUE)
library(imputeFin)


# Documentation
devtools::document()  #to generate all documentation via roxygen
?estimateAR1Gaussian


# Code tests
devtools::test()
