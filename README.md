
# ibdRel

ibdRel is an R package for working with genetic relationship estimation from simulated identity by descent (IBD) segments. 
A Shiny application for estimating and ranking possible pairwise genetic relationships is associated with the package.

## Installation

ibdRel is currently under development. Missing functionality and issues are to be expected.

To install the latest version from GitHub, use:

```r
# install.packages("remotes")
remotes::install_github("hnordto/ibdrel", dependencies = TRUE)
```

Once the installation is complete, the Shiny application can be launched with:

```r
ibdrel::launchApp()
```
