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
library(imputeFin)
#devtools::build()  # to generate the installation file


# Documentation
devtools::document()  #to generate all documentation via roxygen
?XXX
