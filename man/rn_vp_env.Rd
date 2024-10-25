\name{rn_vp_env}
\alias{rn_vp_env}

\title{Compute the total phenotypic variance \eqn{V_{\text{Tot}}}}
\description{
This function calculates the total phenotypic variance \eqn{V_{\text{Tot}}} arising along and beyond the reaction norm, i.e. including the residual variance.
}
\usage{
rn_vp_env(theta, V_theta, var_res, env = NULL, shape = NULL, X = NULL,
          fixed = NULL, wt_env = NULL, width = 10)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{V_theta}{Full variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{var_res}{Residual variance beyond the reaction norm. It could be a scalar if this residual variance is assumed homogeneous or a vector the same length as \code{env} (or number of rows of \code{X}), if the residual variance is allowed to vary across the environment. (numeric scalar or vector)}
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{X}{If the model used was linear in the parameters, the design matrix X of the model (numeric, incompatible with the arguments \code{env} and \code{shape}).}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{wt_env}{Weights to apply to the \code{env} vector values, providing an information regarding their relative probability in the biological context. The weights must non-negative, and at least one must non-zero. The vector \code{wt_env} must be the same length as \code{env}. By default, no weighting is applied. (numeric)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
The phenotypic variance at a given environment depends on the variation in the \eqn{\theta} parameters and the residual variance. It thus account for everything but the variance due to the average shape of the reaction norm (i.e. \eqn{V_{\text{Plas}}}).
}
\value{This function yields the phenotypic variance in each environment (numerical vector).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{rn_vplas}}, \code{\link{rn_vgen}}
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
# Another covariance matrix for cmax and xopt (e.g. permanent environment, or maternal effet)
M     <- matrix(c(0.05,      0,
                  0,      0.01),
                ncol = 2)
# Full variance-covariance matrix
P     <- G + M
# Residual variance
vr    <- 0.1
                
# Computing V_tot
rn_vp_env(theta   = theta,
          V_theta = P,
          var_res = vr,
          env     = vec_env,
          shape   = expr,
          fixed   = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in P
}
