# lisaClust

<img src="https://raw.githubusercontent.com/ellispatrick/lisaClust/master/inst/lisaClust_sticker_ai_upscaled.png" align="right" width="200" />

Clustering of Local Indicators of Spatial Association

## Overview


**lisaClust** provides a series of functions to identify and visualise 
    regions of tissue where spatial associations between cell-types is similar.
    This package can be used to provide a high-level summary of cell-type 
    colocalisation in multiplexed imaging data that has been segmented at a 
    single-cell resolution.

## Installation


For the Bioconductor release version, run the following.
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("lisaClust")
```

If you would like the most up-to-date features, install the most recent development version.
```r
# Install the development version from Bioconductor:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
# This will update all your Bioconductor packages to devel version
BiocManager::install(version='devel')

BiocManager::install("lisaClust")

# Otherwise install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("ellispatrick/lisaClust")
library(lisaClust)
```
## Submitting an issue or feature request

`listClust` is still under active development. We would greatly appreciate any and 
all feedback related to the package.

* R package related issues should be raised [here](https://github.com/ellispatrick/listClust/issues).
* For general questions and feedback, please contact us directly via [ellis.patrick@sydney.edu.au](mailto:ellis.patrick@sydney.edu.au).


## Authors

* **Nicolas Canete**
* **Ellis Patrick**  - [@TheEllisPatrick](https://twitter.com/TheEllisPatrick)

