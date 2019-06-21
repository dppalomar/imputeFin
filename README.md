---
output:
  html_document:
    variant: markdown_github
    keep_md: true
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# imputeFin
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/imputeFin)](https://CRAN.R-project.org/package=imputeFin)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/imputeFin)](https://CRAN.R-project.org/package=imputeFin)
[![CRAN Downloads Total](https://cranlogs.r-pkg.org/badges/grand-total/imputeFin?color=brightgreen)](https://CRAN.R-project.org/package=imputeFin)

Missing values often occur in financial data due to a variety 
of reasons (errors in the collection process or in the processing stage, 
lack of asset liquidity, lack of reporting of funds, etc.). However, 
most data analysis methods expect complete data and cannot be employed 
with missing values. One convenient way to deal with this issue without 
having to redesign the data analysis method is to impute the missing 
values. This package provides an efficient way to impute the missing 
values based on modeling the time series with an autoregressive (AR) model,
convenient to model log-prices and log-volume in financial data. In the
current version, the imputation is univariate-based (so no asset correlation 
is used).

The package is based on the paper:   
J. Liu, S. Kumar, and D. P. Palomar (2019). Parameter Estimation of 
Heavy-Tailed AR Model With Missing Data Via Stochastic EM. _IEEE Trans. on 
Signal Processing_, vol. 67, no. 8, pp. 2159-2172. <doi:10.1109/TSP.2019.2899816>.


## Installation
The package can be installed from [CRAN](https://CRAN.R-project.org/package=imputeFin) or [GitHub](https://github.com/dppalomar/imputeFin):

```r
# install stable version from CRAN
install.packages("imputeFin")

# install development version from GitHub
devtools::install_github("dppalomar/imputeFin")
```

To get help:

```r
library(imputeFin)
help(package = "imputeFin")
?imputeAR1Gaussian
```

To cite `imputeFin` in publications:

```r
citation("imputeFin")
```



## Usage of `imputeFin`
TBD



## Links
Package: [CRAN](https://CRAN.R-project.org/package=imputeFin) and [GitHub](https://github.com/dppalomar/imputeFin).

README file: [GitHub-readme](https://github.com/dppalomar/imputeFin/blob/master/README.md).

Vignette: [CRAN-vignette](https://CRAN.R-project.org/package=imputeFin/vignettes/ImputeFinancialTimeSeries.html).

