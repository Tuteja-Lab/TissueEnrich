# TissueEnrich: A R package to carry out tissue-specific gene enrichment.

**Requirement**

You need R version above 3.4 to run this application. In addition, this package requires the following packages.
*  `dplyr (>= 0.7.3)`
*  `ensurer (>= 1.1)`
*  `ggplot2 (>= 2.2.1)`
*  `utils (>= 3.4.1)`

**How to install the R package**

* Download or fork the bitbucket repository
* Open R terminal or RStudio terminal
* Set the current path to the Bitbucket repository `setwd(<Path of bitbucket repository>)`
* Run command `install.packages("TissueEnrich_1.0.0.tar.gz", repos = NULL, type="source")`
* The package can also be installed with dependencies using `devtools` package
* Set the current path to the Bitbucket repository R package `setwd(<Path of bitbucket repository>/TissueEnrich)`
* Run command `devtools::install(".",dependencies = TRUE)`
* Load Library `library(TissueEnrich)`

**More about the package**

* Check more details about the package in the vignette `vignette(“TissueEnrich”)`
