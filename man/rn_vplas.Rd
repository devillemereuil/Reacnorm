\name{rn_vplas}
\alias{rn_vplas}

\title{Compute the purely environmental plastic variance \eqn{V_{\text{Plas}}}}
\description{
This function calculates the purely environmental plastic variance arising from the average (across genotypes) shape of the reaction norm.
}
\usage{
rn_vplas(env, shape, theta, G_theta, fixed = NULL, wt_env = NULL,
         correction = FALSE, width = 10)
}
\arguments{
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{wt_env}{Weights to apply to the \code{env} vector values, providing an information regarding their relative probability in the biological context. The weights must non-negative, and at least one must non-zero. The vector \code{wt_env} must be the same length as \code{env}. By default, no weighting is applied. (numeric)}
  \item{correction}{Should Bessel's correction (dividing the sum-of-squares by N-1 rather than N) be used (TRUE) or not (FALSE, default). The default is FALSE, because it is likely that other components such as the total phenotypic variance is computed over the number of individuals (generally large number) rather than the number of environments (generally small number). The best is to manually use Bessel's correction over the proper number of data points. (boolean)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
The variance \eqn{V_{\text{Plas}}} is the variance arising from the shape of the reaction norm after avering over the genetic variance, i.e. the purely environmental part of the variance of plasticity in the reaction norm.

It is very important that the parameters are in the same order in \code{theta} (which, again, must be named) and in \code{G_theta}.
}
\value{This function yields \eqn{V_{\text{Plas}}} (numeric).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{rn_mean_by_env}}, \code{\link{rn_pi_decomp}}
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

# Computing V_plas
rn_vplas(env     = vec_env,
         shape   = expr,
         theta   = theta,
         G_theta = G,
         fixed   = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in G

# Note the quite large difference when Bessel's correction is used
rn_vplas(env     = vec_env,
         shape   = expr,
         theta   = theta,
         G_theta = G,
         fixed   = c(3, 4),
         correction = TRUE)
         
# It is possible to weight the environment, e.g. according to a normal distribution
rn_vplas(env     = vec_env,
         shape   = expr,
         theta   = theta,
         G_theta = G,
         fixed   = c(3, 4),
         wt_env  = dnorm(vec_env),
         correction = TRUE)
}
