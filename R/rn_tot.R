###############################################################################
##                             ReacNorm R package                            ##
##                   Functions to compute the total phenotypic               ##
##                          variance of the reaction norm                    ##
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

## Compute the reaction norm variance conditionnally to the environment
# Args: - e: the environmental value to condition to
#       - func: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - P_theta: the complete variance-covariance matrix of the parameters theta
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_g_e (numeric)
rn_vt_e <- function(e, func, theta, P_theta, fixed = NULL, width = 10) {
    
    # Handling when some terms are fixed
    if (!is.null(fixed)) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(P_theta) == length(theta)) {
            P_theta <- P_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }
    
    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(P_theta)) * width
    
    # Number of dimensions
    d <- length(w)
    
    # Computing the logdet of vcov
    logdet <- calc_logdet(P_theta)
    
    # Average
    avg <- rn_avg_e(e, func, theta, P_theta, fixed = fixed, width = 10)
    
    # Computing the integral
    cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!is.null(fixed)) { full_x[var, ] <- x } else { full_x <- x }
            func(e, full_x)^2 * vec_mvnorm(x, var_theta, P_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral - avg^2
}

##  --------------------------------------------------- Frontend functions ----
    
## Compute the total phenotypic variance in the reaction norm
# Args: - theta: the average parameters estimated by the model (must be named)
#       - P_theta: the complete variance-covariance matrix of the parameters theta
#       - V_res: residual variance or vector of residuals variances
#       - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - X: the design matrix if the model used was linear (incompatible with "shape")
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env or rows of X)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_gen (numeric)
rn_vtot <- function(theta,
                    P_theta,
                    V_res,
                    env = NULL,
                    shape = NULL,
                    X = NULL,
                    fixed = NULL,
                    wt_env = NULL,
                    average = TRUE,
                    width = 10) {
    # The parameter theta must be named
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }
    
    # Checking that V_res is of length 1 or length(env)
    if (!(length(V_res) %in% c(1, length(env), nrow(X)))) {
        stop("V_res should be either a single value, or same number of values as the environment")
    }

    # X is incompatible with shape
    if (is.null(X) & (is.null(shape) & is.null(env))) {
        stop("Either the shape and environment or the design matrix X of a reaction norm should be provided")
    } else if (!is.null(X) & !(is.null(shape) & is.null(env))) {
        stop("The arguments X and shape cannot be used together.\n If the design matrix is available, it is usually better to use the argument X.")
    }

    if (!is.null(X)) {
        # Check X and theta compatibility
        if (ncol(X) != length(theta)) {
            stop("The number of columns in X should be equal to the length of theta.")
        }
    }

    # Use weighted mean if wt_env is not NULL
    if (is.null(wt_env)) {
        func_mean <- mean
    } else {
        func_mean <- function(x) { weighted.mean(x, w = wt_env) }
    }

    if (!is.null(X)) {
        # Computing the row-by-row genetic variance
        out <- apply(X, 1, \(row_) { t(row_) %*% P_theta %*% row_ })
    } else {
        # Formatting the shape for the computation of the integral
        func <- rn_generate_shape(shape, names(theta))

        # Computing V_gen for each environment
        out <-
            sapply(env,
                   \(e) rn_vt_e(e       = e,
                                func    = func,
                                theta   = theta,
                                P_theta = P_theta,
                                fixed   = fixed,
                                width   = width)) +
            V_res
    }

    # Averaging if requested (default)
    if (average) { out <- func_mean(out) }

    return(out)
}
