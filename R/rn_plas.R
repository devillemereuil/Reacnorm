###############################################################################
##                             ReacNorm R package                            ##
##             Functions to compute V_plas and related parameters            ##
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

## Compute the average phenotype conditional to the environment
# Args: - e: the environmental value to condition to
#       - func: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for E_g_e (numeric)
rn_avg_e <- function(e,
                     func,
                     theta,
                     G_theta,
                     fixed = NULL,
                     width = 10) {
    # Handling when some terms are fixed
    if (!is.null(fixed)) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

    # Average
    avg <- cubature::hcubature(
        f  = function(x) {
            full_x <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!is.null(fixed)) { full_x[var, ] <- x } else { full_x <- x }
            func(e, full_x) * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit      = var_theta - w,
        upperLimit      = var_theta + w,
        fDim            = 1,
        tol             = 0.001,
        absError        = 0.0001,
        vectorInterface = TRUE
    )$integral

    return(avg)
}

##  --------------------------------------------------- Frontend functions ----

## Compute the mean phenotype by environment
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model (must be named)
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_plas (numeric)
rn_mean_by_env <- function(env,
                           shape,
                           theta,
                           G_theta,
                           fixed = NULL,
                           width = 10) {
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }

    # Formatting the function of the shape for the computation of the integral
    func <- rn_generate_shape(shape, names(theta))

    sapply(env,
           \(e) rn_avg_e(e       = e,
                         func    = func,
                         theta   = theta,
                         G_theta = G_theta,
                         fixed   = fixed,
                         width   = width))
}

## Compute the plastic variance (V_plas)
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model (must be named)
#       - G_theta: the G matrix containing the genetic variance-covariances of
#                  the parameters
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env)
#       - correction: should the Bessel correction be used or not?
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_plas (numeric)
rn_vplas <- function(env,
                     shape,
                     theta,
                     G_theta,
                     fixed = NULL,
                     wt_env = NULL,
                     correction = FALSE,
                     width = 10) {
    # theta must be named to know the names of parameters to format the "shape"
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }

    # Configure variance function
    method <- ifelse(correction, "unbiased", "ML")
    if (is.null(wt_env)) {wt_env <- rep(1/length(env), length(env))}
    func_var <- function(x) {
        cov.wt(cbind(x), wt = wt_env, method = method)[["cov"]] |>
            as.numeric()
    }

    # Compute the average for each environment, then take the variance
    rn_mean_by_env(env      = env,
                   shape    = shape,
                   theta    = theta,
                   G_theta  = G_theta,
                   fixed    = fixed,
                   width    = width) |>
        func_var()
}

## Compute the π-decomposition of V_Plas
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model (must be named)
#       - G_theta: the G matrix containing the genetic variance-covariances of
#                  the parameters
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env)
#       - correction: should the Bessel correction be used or not?
#       - v_plas: optionnaly, provide the value for v_plas if already computed
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_plas (numeric)
rn_pi_decomp <- function(env,
                         shape,
                         theta,
                         G_theta,
                         fixed = NULL,
                         wt_env = NULL,
                         correction = FALSE,
                         v_plas = NA,
                         width = 10) {
    # theta must be named to know the names of parameters to format the "shape"
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }

    # Use weighted mean if wt_env is not NULL
    if (is.null(wt_env)) {
        func_mean    <- mean
    } else {
        func_mean    <- function(x) { weighted.mean(x, w = wt_env) }
    }

    # Configure variance function
    method <- ifelse(correction, "unbiased", "ML")
    if (is.null(wt_env)) {wt_env <- rep(1/length(env), length(env))}
    func_var <- function(x) {
        cov.wt(cbind(x), wt = wt_env, method = method)[["cov"]] |>
            as.numeric()
    }

    # Compute V_plas
    if (is.na(v_plas)) {
        v_plas <-
            rn_vplas(env        = env,
                     shape      = shape,
                     theta      = theta,
                     G_theta    = G_theta,
                     fixed      = fixed,
                     correction = correction,
                     width      = width)
    }

    # Compute the average slope
    d_func <- rn_generate_gradient(shape, "x", names(theta))

    mean_sl <-
        sapply(env,
               \(e) rn_avg_e(e       = e,
                             func    = d_func,
                             theta   = theta,
                             G_theta = G_theta,
                             fixed   = fixed,
                             width   = width)) |>
        func_mean()

    # Compute the variance due to the average slope
    var_sl <- (mean_sl)^2 * func_var(env)

    # Compute the average curvature
    d2_func <- rn_generate_2diff(shape, "x", names(theta))

    mean_cv <-
        sapply(env,
               \(e) rn_avg_e(e       = e,
                             func    = d2_func,
                             theta   = theta,
                             G_theta = G_theta,
                             fixed   = fixed,
                             width   = width)) |>
        func_mean()

    # Compute the variance due to the average curvature*
    var_cv <- 0.25 * (mean_cv)^2 * func_var(env^2)

    # Return the decomposition
    data.frame(V_Plas = v_plas,
               Pi_Sl  = var_sl / v_plas,
               Pi_Cv  = var_cv / v_plas)
}
