###############################################################################
##                             ReacNorm R package                            ##
##      Functions to be used in the integration procedures using cubature    ##
##       ----------------------------------------------------------------    ##
##                           Pierre de Villemereuil                          ##
##       ----------------------------------------------------------------    ##
##                                     2024                                  ##
###############################################################################

## --------------------------------------------------------------- LICENCE ----

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>. 

##  ---------------------------------------------------- Backend functions ----

## Log-determinant of a VCV matrix
calc_logdet <- function(Sigma) {
    eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values |>
        log() |> 
        sum()
}

## Vectorised multivariate Gaussian density function
# Shamelessly stolen from cubature vignette (credit to Balasubramanian Narasimhan)
vec_mvnorm <- function(x, mean, Sigma, logdet = NULL) {
    # If logdet not provided, compute it
    if (is.null(logdet)) {
        logdet <- calc_logdet(Sigma)
    }
    # Compute Mahalanobis distance (corresponds to the exp. part of the density)
    distval <- stats::mahalanobis(t(x), center = mean, cov = Sigma)
    # Compute the vectorised MVN density
    out <- exp(matrix(-(nrow(x) * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
    return(out)
}

## Vectorised function to get upper-triangle squared matrix
vec_sq_uptri <- function(x) {
    apply(x, 2, \(col_) {
        mat <- (col_) %*% t(col_)
        return(mat[upper.tri(mat, diag = TRUE)])
    })
}
