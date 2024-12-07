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
#       - V_theta: the complete variance-covariance matrix of the parameters theta
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_g_e (numeric)
rn_vt_e <- function(e, func, theta, V_theta, fixed = NULL, width = 10) {
    
    # Handling when some terms are fixed
    if (!is.null(fixed)) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(V_theta) == length(theta)) {
            V_theta <- V_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }
    
    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(V_theta)) * width
    
    # Number of dimensions
    d <- length(w)
    
    # Computing the logdet of vcov
    logdet <- calc_logdet(V_theta)
    
    # Average
    avg <- rn_avg_e(e, func, theta, V_theta, fixed = fixed, width = 10)
    
    # Computing the integral
    cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!is.null(fixed)) { full_x[var, ] <- x } else { full_x <- x }
            func(e, full_x)^2 * vec_mvnorm(x, var_theta, V_theta, logdet)
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
    
## Compute the phenotypic variance in each environment
# Args: - theta: the average parameters estimated by the model (must be named)
#       - V_theta: the complete variance-covariance matrix of the parameters theta
#       - var_res: residual variance or vector of residuals variances
#       - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - X: the design matrix if the model used was linear (incompatible with "shape")
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env or rows of X)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_gen (numeric)
rn_vp_env <- function(theta,
                      V_theta,
                      var_res,
                      env = NULL,
                      shape = NULL,
                      X = NULL,
                      fixed = NULL,
                      wt_env = NULL,
                      width = 10) {
    # The parameter theta must be named
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }
    
    # Checking that var_res is of length 1 or length(env)
    if (!(length(var_res) %in% c(1, length(env), nrow(X)))) {
        stop("var_res should be either a single value, or same number of values as the environment")
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
        out <- apply(X, 1, \(row_) { t(row_) %*% V_theta %*% row_ })
    } else {
        # Formatting the shape for the computation of the integral
        func <- rn_generate_shape(shape, names(theta))

        # Computing V_gen for each environment
        out <-
            sapply(env,
                   \(e) rn_vt_e(e       = e,
                                func    = func,
                                theta   = theta,
                                V_theta = V_theta,
                                fixed   = fixed,
                                width   = width)) +
            var_res
    }

    return(out)
}

## Compute the total phenotypic variance
# Args: - theta: the average parameters estimated by the model (must be named)
#       - V_theta: the complete variance-covariance matrix of the parameters theta
#       - var_res: residual variance or vector of residuals variances
#       - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - X: the design matrix if the model used was linear (incompatible with "shape")
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env or rows of X)
#       - correction: should the Bessel correction be used or not?
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_gen (numeric)
rn_vtot <- function(theta,
                    V_theta,
                    var_res,
                    env = NULL,
                    shape = NULL,
                    X = NULL,
                    S = NULL,
                    fixed = NULL,
                    wt_env = NULL,
                    correction = FALSE,
                    width = 10) {
    # The parameter theta must be named
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }
    
    # Checking that var_res is of length 1 or length(env)
    if (!(length(var_res) %in% c(1, length(env), nrow(X)))) {
        stop("var_res should be either a single value, or same number of values as the environment")
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
    method <- ifelse(correction, "unbiased", "ML")
    if (is.null(wt_env)) {
        func_mean <- mean
    } else {
        func_mean <- function(x) { weighted.mean(x, w = wt_env) }
    }
    if (is.null(wt_env) & !is.null(env)) {
        wt_env <- rep(1/length(env), length(env))
    } else if (!is.null(X)) {
        wt_env <- rep(1/nrow(X), nrow(X))
    }
    func_var <- function(x) {
        cov.wt(cbind(x), wt = wt_env, method = method)[["cov"]]
    }
    
    # Compute V_tot by environment
    if (!is.null(X)) {
        # Computing the row-by-row genetic variance
        vtot_env <- apply(X, 1, \(row_) { t(row_) %*% V_theta %*% row_ })
    } else {
        # Formatting the shape for the computation of the integral
        func <- rn_generate_shape(shape, names(theta))

        # Computing V_gen for each environment
        vtot_env <-
            sapply(env,
                   \(e) rn_vt_e(e       = e,
                                func    = func,
                                theta   = theta,
                                V_theta = V_theta,
                                fixed   = fixed,
                                width   = width)) +
            var_res
    }
    
    # Compute V_Plas
    if (!is.null(X)) {
        # Compute the variance-covariance matrix of X
        cov_X <- func_var(X)
        # Compute the correcting factor due to the uncertainty
        var_uncert <- sum(cov_X * S)
        # Compute V_Plas
        vplas <- as.vector(t(theta) %*% cov_X %*% theta - var_uncert)
        
    } else {
        # Compute the average for each environment, then take the variance
        vplas <-
            rn_mean_by_env(theta    = theta,
                           V_theta  = V_theta,
                           env      = env,
                           shape    = shape,
                           fixed    = fixed,
                           width    = width) |>
            func_var() |>
            as.numeric()
    }
    
    # Averaging if requested (default)
    vtot <- func_mean(vtot_env) + vplas

    return(vtot)
}
