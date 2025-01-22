# Reacnorm: an R package to decompose the variance in reaction norms

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/Reacnorm)](https://cran.r-project.org/package=Reacnorm)
[![DOI](https://zenodo.org/badge/827684926.svg)](https://doi.org/10.5281/zenodo.14674905)

## What is this package?

Reacnorm provides various functions to partition the phenotypic variance of a plastic trait, studied through its reaction norm. The variance partition distinguishes between the variance arising from the average shape of the reaction norms (V_Plas) and the (additive) genetic variance . The latter is itself separated into an environment-blind component (V_G/V_A) and the component arising from plasticity (V_GxE/V_AxE). The package also provides a way to further partition V_Plas into aspects (slope/curvature) of the shape of the average reaction norm (π-decomposition) and partition V_Add (γ-decomposition) and V_AxE (ι-decomposition) into the impact of genetic variation in the reaction norm parameters.

This package is still in development and the article related to it ~~is still in peer review~~ has been recommended by PCI EvolBiol and is on its way to Peer Community Journal! It is being submitted to CRAN.

### Reference

The current reference is the recommended preprint of [de Villemereuil & Chevin (2025)](https://doi.org/10.32942/X2NC8B).

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

