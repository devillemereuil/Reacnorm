###############################################################################
##                             ReacNorm R package                            ##
##              Functions to compute the (additive) genetic variances        ##
##                           and related parameters                          ##
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

## Compute the genetic variance conditionnally to the environment
# Args: - e: the environmental value to condition to
#       - func: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_g_e (numeric)
rn_vg_e <- function(e, func, theta, G_theta, fixed = NULL, width = 10) {

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
    avg <- rn_avg_e(e, func, theta, G_theta, fixed = fixed, width = 10)

    # Computing the integral
    cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!is.null(fixed)) { full_x[var, ] <- x } else { full_x <- x }
            func(e, full_x)^2 * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral - avg^2
}

## Compute the average gradient conditionnally to the environment
# Args: - e: the environmental value to condition to
#       - d_func: the differential of function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_g_e (numeric)
rn_psi_e <- function(e, d_func, theta, G_theta, fixed = NULL, width = 10) {

    # Handling when some terms are fixed
    if (!is.null(fixed)) {
        var        <- setdiff(1:length(theta), fixed)
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

    # Computing the integral
    cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!is.null(fixed)) { full_x[var, ] <- x } else { full_x <- x }
            d_func(e, full_x) * matrix(rep(vec_mvnorm(x, var_theta, G_theta, logdet), d),
                                       nrow = d,
                                       byrow = TRUE)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = d,
        tol        = 0.001,
        absError   = 0.0001,
#         maxEval    = 10^7,
        vectorInterface = TRUE
    )$integral
}

##  --------------------------------------------------- Frontend functions ----

## Compute the genetic variance (V_gen)
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model (must be named)
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_gen (numeric)
rn_vgen <- function(env,
                    shape,
                    theta,
                    G_theta,
                    fixed = NULL,
                    wt_env = NULL,
                    average = TRUE,
                    width = 10) {
    # The parameter theta must be named
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }

    # Use weighted mean if wt_env is not NULL
    if (is.null(wt_env)) {
        func_mean <- mean
    } else {
        func_mean <- function(x) { weighted.mean(x, w = wt_env) }
    }

    # Formatting the shape for the computation of the integral
    func <- rn_generate_shape(shape, names(theta))

    # Computing V_gen for each environment
    out <-
        sapply(env,
               \(e) rn_vg_e(e       = e,
                            func    = func,
                            theta   = theta,
                            G_theta = G_theta,
                            fixed   = fixed,
                            width   = width))

    # Averaging if requested (default)
    if (average) { out <- func_mean(out) }

    return(out)
}

## Compute the additive genetic variance decomposition
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model (must be named)
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - wt_env: weights to use for computing the average over env
#                 (must be the same length as env)
#       - compute_gamma: should the γ-decomposition of V_add be returned?
#       - compute_iota: should the ι-decomposition of V_AxE be returned?
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_gen (numeric)
rn_gen_decomp <- function(env,
                          shape,
                          theta,
                          G_theta,
                          fixed = NULL,
                          wt_env = NULL,
                          compute_gamma = TRUE,
                          compute_iota  = TRUE,
                          width = 10) {
    # The parameter theta must be named
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }

    # Getting the names of parameters
    if (is.null(fixed)) {
        all_names <- names(theta)
        var_names <- names(theta)
    } else {
        all_names <- names(theta)
        var_names <- names(theta)[-fixed]
    }

    # Use weighted mean if wt_env is not NULL
    if (is.null(wt_env)) {
        func_mean    <- mean
        func_colmean <- colMeans
        func_cov     <- function(x) { cov(x) * (nrow(x) - 1) / nrow(x)}
    } else {
        func_mean    <- function(x) { weighted.mean(x, w = wt_env) }
        func_colmean <- function(x) { colWeightedMeans(x, w = wt_env) }
        func_cov     <- function(x) { cov.wt(x, wt = wt_env, method = "ML")[["cov"]] }
    }

    # Formatting the shape for the computation of the integral
    func    <- rn_generate_shape(shape, all_names)
    d_func  <- rn_generate_gradient(shape, var_names, all_names)

    # Computing psi for each environment
    psi <-
        lapply(env,
               \(e) {rn_psi_e(e       = e,
                              d_func  = d_func,
                              theta   = theta,
                              G_theta = G_theta,
                              fixed   = fixed,
                              width   = width)})
    psi <- do.call("rbind", psi)
    colnames(psi) <- var_names

    # Computing the total additive genetic variance V_add
    v_add <-
        apply(psi, 1, \(psi_) t(psi_) %*% G_theta %*% psi_) |>
        func_mean()
    names(v_add) <- "V_Add"

    # Computing the marginal additive genetic variance V_A
    v_a <- (t(func_colmean(psi)) %*% G_theta %*% func_colmean(psi)) |>
           as.vector()
    names(v_a) <- "V_A"

    # Computing the interaction variance
    v_axe <- v_add - v_a
    names(v_axe) <- "V_AxE"

    # Computing the γ-decomposition if required (default)
    if (compute_gamma) {
        # Compute the direct γ_i for each parameters
        gamma_i <- t(apply(psi, 1, \(psi_) psi_^2 * diag(G_theta)))
        colnames(gamma_i) <- paste0("Gamma_",var_names)

        # Compute the γ_ij for each pair of parameters
        gamma_ij <- apply(psi, 1, simplify = FALSE, FUN = \(psi_) {
            out <- numeric(sum(upper.tri(G_theta)))
            names_out <- character(sum(upper.tri(G_theta)))
            k   <- 1
            for (i in 1:(length(psi_) - 1)) {
                for (j in (i+1):length(psi_)) {
                    out[k] <- 2 * psi_[i] * psi_[j] * G_theta[i, j]
                    names_out[k] <- paste(names(psi_)[i], names(psi_)[j], sep = "_")
                    k <- k + 1
                }
            }
            names(out) <- names_out
            return(out)
        })
        gamma_ij <- do.call("rbind", gamma_ij)
        colnames(gamma_ij) <- paste0("Gamma_",colnames(gamma_ij))

        # Returning the value for the γ-decomposition
        gamma <-
            cbind(
                gamma_i,
                gamma_ij
            ) |>
            func_colmean()
        gamma <- gamma / v_add
    } else {
        gamma <- NULL
    }

    # Computing the ι-decomposition if required (default)
    if (compute_iota) {
        # Computing the (weighted) variance-covariance of Psi
        # (needs to remove Bessel's correction for consistancy with the means above)
        Psi <- func_cov(psi)

        # Computing the pairwise multiplication of Psi and G_theta
        M <- Psi * G_theta
        colnames(M) <- rownames(M) <- var_names

        # Computing the direct ι_i estimates
        iota_i <- diag(M)
        names(iota_i) <- paste0("Iota_", names(iota_i))

        # Compute the ι_ij for each pair of parameters
        iota_ij <- 2 * M[upper.tri(M)]
        M_names <- matrix(paste(matrix(colnames(M),
                                       nrow = nrow(M),
                                       ncol = ncol(M)),
                                matrix(colnames(M),
                                       nrow = nrow(M),
                                       ncol = ncol(M),
                                       byrow = TRUE),
                                sep = "_"),
                          ncol = ncol(M),
                          nrow = nrow(M))
        names(iota_ij) <- M_names[upper.tri(M_names)]
        names(iota_ij) <- paste0("Iota_", names(iota_ij))

        # Returning the value for the ι-decomposition
        iota <- c(iota_i, iota_ij) / v_axe
    } else {
        iota <- NULL
    }

    # Formatting the output
    out <- cbind(V_Add = v_add,
                 V_A = v_a,
                 V_AxE = v_axe,
                 do.call("cbind", as.list(gamma)),
                 do.call("cbind", as.list(iota))) |>
           as.data.frame()
    rownames(out) <- NULL
    return(out)
}

## Performing the γ-decomposition for each environment
# Args: - env: the environmental values over which the model has been estimated
#       - shape: the shape of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#                (must be named)
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - width: the width over which the integral must be computed
#                (10 is a generally a good value)
# Value: The value for V_A and Gamma-decomposition as a data.frame for each environment
rn_gamma_env <- function(env,
                         shape,
                         theta,
                         G_theta,
                         fixed = NULL,
                         width = 10) {
    # The parameter theta must be named
    if (is.null(names(theta))) {
        stop("The vector theta must be named with the corresponding parameter names")
    }

    # Getting the names of parameters
    if (is.null(fixed)) {
        all_names <- names(theta)
        var_names <- names(theta)
    } else {
        all_names <- names(theta)
        var_names <- names(theta)[-fixed]
    }

    # Formatting the shape for the computation of the integral
    func    <- rn_generate_shape(shape, all_names)
    d_func  <- rn_generate_gradient(shape, var_names, all_names)

    # Computing psi for each environment
    psi <-
        lapply(env,
               \(e) {rn_psi_e(e       = e,
                              d_func  = d_func,
                              theta   = theta,
                              G_theta = G_theta,
                              fixed   = fixed,
                              width   = width)})
    psi <- do.call("rbind", psi)
    colnames(psi) <- var_names

    # Computing V_A_e for each environment
    v_a_e <- apply(psi, 1, \(psi_) t(psi_) %*% G_theta %*% psi_)

    # Computing the gamma-decomposition of V_A_e
    gamma_i <- t(apply(psi, 1, \(psi_) psi_^2 * diag(G_theta)))
    colnames(gamma_i) <- paste0("Gamma_",colnames(gamma_i))

    gamma_ij <- apply(psi, 1, simplify = FALSE, FUN = \(psi_) {
        out <- numeric(sum(upper.tri(G_theta)))
        names_out <- character(sum(upper.tri(G_theta)))
        k   <- 1
        for (i in 1:(length(psi_) - 1)) {
            for (j in (i+1):length(psi_)) {
                out[k] <- 2 * psi_[i] * psi_[2] * G_theta[i, j]
                names_out[k] <- paste(names(psi_)[i], names(psi_)[j], sep = "_")
                k <- k + 1
            }
        }
        names(out) <- names_out
        return(out)
    })
    gamma_ij <- do.call("rbind", gamma_ij)
    colnames(gamma_ij) <- paste0("Gamma_",colnames(gamma_ij))

    # Formatting the output
    out <-
        cbind(
            gamma_i,
            gamma_ij
        )
    out <- out / v_a_e
    out <- cbind(data.frame(Env = env, V_A = v_a_e), out)

    return(out)
}
