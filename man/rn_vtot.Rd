\name{rn_vtot}
\alias{rn_vtot}

\title{Compute the total phenotypic variance \eqn{V_{\text{Tot}}}}
\description{
This function calculates the total phenotypic variance \eqn{V_{\text{Tot}}} arising along and beyond the reaction norm, i.e. including the residual variance.
}
\usage{
rn_vtot(theta, V_theta, var_res, env = NULL, shape = NULL, X = NULL, S = NULL,
        fixed = NULL, wt_env = NULL, correction = FALSE, width = 10)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{V_theta}{Full variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{var_res}{Residual variance beyond the reaction norm. It could be a scalar if this residual variance is assumed homogeneous or a vector the same length as \code{env} (or number of rows of \code{X}), if the residual variance is allowed to vary across the environment. (numeric scalar or vector)}
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{X}{If the model used was linear in the parameters, the design matrix X of the model (numeric, incompatible with the arguments \code{env} and \code{shape}).}
  \item{S}{The error variance-covariance matrix of the estimated fixed effects \code{theta}. Used only if the \code{X} argument is used. (numeric, incompatible with the arguments \code{env} and \code{shape}).}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{wt_env}{Weights to apply to the \code{env} vector values, providing an information regarding their relative probability in the biological context. The weights must non-negative, and at least one must non-zero. The vector \code{wt_env} must be the same length as \code{env}. By default, no weighting is applied. (numeric)}
  \item{correction}{Should Bessel's correction (dividing the sum-of-squares by N-1 rather than N) be used (TRUE) or not (FALSE, default). The default is FALSE, because it is likely that other components such as the total phenotypic variance is computed over the number of individuals (generally large number) rather than the number of environments (generally small number). The best is to manually use Bessel's correction over the proper number of data points. (boolean)}  
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
The variance \eqn{V_{\text{Tot}}} is the total phenotypic variance of the trait along and beyond the reaction norm. In the simplest case, it is equal to the sum \eqn{V_{\text{Plas}} + V_{\text{Gen}} + V_{\text{Res}}}, in which case this function is rather useless. But in cases where the parameters \eqn{\theta} in the reaction norm curve is fitted with more than the (additive) genetic effect, then computing \eqn{V_{\text{Tot}}} becomes non-trivial. This function allows to do it, e.g. by providing (i) a way to specify the full variation assumed in \eqn{\theta}, using the \code{V_theta} full variance-covariance argument, and (ii) a way to specify a varying residual variance by environment.
}
\value{This function yields the total phenotypic variance along and beyond the reaction norm (numeric).}

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
rn_vtot(theta   = theta,
        V_theta = P,
        var_res = vr,
        env     = vec_env,
        shape   = expr,
        fixed   = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in P

# With varying residual variance along the environmental (here with environmental variance)
vr    <- 0.05 * (vec_env^2)
rn_vtot(theta   = theta,
        V_theta = P,
        var_res = vr,
        env     = vec_env,
        shape   = expr,
        fixed   = c(3, 4))
# Same result because the average of Vr is the same as the unique value above
}
