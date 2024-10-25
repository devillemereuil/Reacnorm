\name{rn_vgen}
\alias{rn_vgen}

\title{Compute the total genetic variance \eqn{V_{\text{Gen}}}}
\description{
This function calculates the total genetic variance \eqn{V_{\text{Gen}}} arising from genetic variation in the reaction norm.
}
\usage{
rn_vgen(theta, G_theta, env = NULL, shape = NULL, X = NULL, fixed = NULL, wt_env = NULL,
        average = TRUE, width = 10)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{X}{If the model used was linear in the parameters, the design matrix X of the model (numeric, incompatible with the arguments \code{env} and \code{shape}).}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{wt_env}{Weights to apply to the \code{env} vector values, providing an information regarding their relative probability in the biological context. The weights must non-negative, and at least one must non-zero. The vector \code{wt_env} must be the same length as \code{env}. By default, no weighting is applied. (numeric)}
  \item{average}{Should the values for each environment be given, or just their average (default)? (boolean)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
The variance \eqn{V_{\text{Gen}}} is the total genetic variance in the reaction norm. If the reaction norm is linear in its parameters, then is it also equal to the total additive genetic variance in the reaction norm (\eqn{V_{\text{Add}}}). Otherwise, see the function \code{\link{rn_gen_decomp}} to compute \eqn{V_{\text{Add}}}.

It is very important that the parameters are in the same order in the argument of the \code{shape} function, in \code{theta} and in \code{G_theta}.

}
\value{This function yields the genetic variance arising from all the genetic variation in the reaction norm (numeric).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{rn_gen_decomp}}, \code{\link{rn_vplas}}, \code{\link{rn_vtot}}
}
\examples{
# Some environment vector
vec_env <- seq(-2, 2)

# Shape function
expr <- expression(
     cmax * exp(
         - exp(rho * (x - xopt) - 6) -
             sigmagaus * (x - xopt)^2
     ))

# Theta
theta <- c(cmax = 1, xopt = 0.9, rho = 8, sigmagaus = 0.4)
# G, only for cmax and xopt
G     <- matrix(c(0.1,      0.01,
                  0.01,     0.05),
                ncol = 2)

# Computing V_gen
rn_vgen(theta   = theta,
        G_theta = G,
        env     = vec_env,
        shape   = expr,
        fixed   = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in G

# Applying some weighting according to a normal distribution
rn_vgen(theta   = theta,
        G_theta = G,
        env     = vec_env,
        shape   = expr,
        fixed   = c(3, 4),
        wt_env  = dnorm(vec_env))

# Setting average to FALSE allows to obtain the value for each environment
rn_vgen(theta   = theta,
        G_theta = G,
        env     = vec_env,
        shape   = expr,
        fixed   = c(3, 4),
        average = FALSE)

# If a polynomial was used, it is possible to use the linear modeling rather having
# to compute integrals
theta <- c(a = 1.5, b = 0.5, c = -0.5)
X     <- cbind(1, vec_env, (vec_env - mean(vec_env))^2)
G     <- 0.1 * diag(3)
rn_vgen(theta   = theta,
        G_theta = G,
        X       = X)
# Should be very close to the computation with integration
rn_vgen(theta   = theta,
        G_theta = G,
        env     = vec_env,
        shape   = expression(a + b * x + c * x^2))
}
