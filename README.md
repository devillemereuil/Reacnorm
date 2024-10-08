# Reacnorm: an R package to decompose the variance in reaction norms

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Reacnorm)](https://cran.r-project.org/package=Reacnorm)

## What is this package?

Reacnorm allows for the computation of various quantitative genetics parameters from a statistical models used to study reaction norms.

This package is still in development and the article related to it is still in peer review. It has not been submitted to CRAN yet, but will be in the future.

## How to install this package

### Using CRAN (not for now, coming soon)
* ~~Simply use `install.packages("Reacnorm")` as for any package.~~

### From this GitHub

* Install the `remotes` package:
     `install.packages("remotes")`
* Install the packages on which Reacnorm depends: cubature, stringi and matrixStats.
    `install.packages(c("cubature","stringi", "matrixStats", "R.rsp"))`
* Using the `remotes` package, install Reacnorm directly from this Github:
    `remotes::install_github("devillemereuil/Reacnorm", build_vignette = TRUE)`

## Tutorial

You can access the tutorial of the package [here](https://github.com/devillemereuil/Reacnorm/blob/main/vignettes/TutoReacnorm.pdf), or by using `vignette("TutoReacnorm")` once the package is installed.

## Submit feedback

If you encounter any bug or usability issue, or if you have some suggestions or feature request, please use the [issue tracker](https://github.com/devillemereuil/Reacnorm/issues). Thank you!

