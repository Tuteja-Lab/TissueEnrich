# TissueEnrich: R package to carry out tissue-specific gene enrichment.

[![bioc](http://www.bioconductor.org/shields/years-in-bioc/TissueEnrich.svg)](http://bioconductor.org/packages/stats/bioc/TissueEnrich.html)
[![bioc](http://www.bioconductor.org/shields/build/devel/bioc/TissueEnrich.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/TissueEnrich/)
<!--[![bioc](http://www.bioconductor.org/shields/downloads/TissueEnrich.svg)](http://bioconductor.org/packages/stats/bioc/TissueEnrich/)-->

**Citation**

If you use TissueEnrich in published research, please cite:

> Jain, A, Tuteja, G. (2018)
> TissueEnrich: Tissue-specific gene enrichment analysis.
> *Bioinformatics*, **bty890**.
> [10.1093/bioinformatics/bty890](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty890/5140218)

**Requirement**

You need R version above 3.5 to run this application. Other dependencies are:

* `dplyr (>= 0.7.3)`
* `ensurer (>= 1.1)`
* `ggplot2 (>= 2.2.1)`
* `tidyr (>=0.8.0)`
* `SummarizedExperiment (>= 1.6.5)`
* `GSEABase (>= 1.38.2)`

**How to install the R package**

**From Bioconductor**
* source("https://bioconductor.org/biocLite.R") 
* biocLite("TissueEnrich") 

**From Github**
* Install Dependencies
* `install.packages(c("dplyr","ensurer","ggplot2","tidyr"))`
* `install.packages("BiocManager")`
* `BiocManager::install("SummarizedExperiment")`
* `BiocManager::install("GSEABase")`
* Now install the `devtools` package
* `install.packages("devtools")`
* `library(devtools)`
* Run command `install_github("Tuteja-Lab/TissueEnrich")`

**More about the package**

* Check more details about the package in the vignette `vignette("TissueEnrich")`
