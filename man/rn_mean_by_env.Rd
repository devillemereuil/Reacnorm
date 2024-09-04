\name{rn_mean_by_env}
\alias{rn_mean_by_env}

\title{Compute the phenotypic mean on the observed scale}
\description{
This function calculates the phenotypic mean on the observed scale from the latent mean and variance.
}
\usage{
rn_mean_by_env(theta, G_theta, env, shape, fixed = NULL, width = 10)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
This function yields the phenotypic average (across genotypes) for each value of the environment. If the reaction norm is not linear in its parameters (the parameters in \code{theta}), then this average will generaly differ from the value yielded by \code{shape} evaluated at values in \code{theta}.

It is very important that the parameters are in the same order in \code{theta} (which, again, must be named) and in \code{G_theta}.
}
\value{This function yields the phenotypic mean for each value of the environmental vector provided (numeric).}

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

# Computing mean by environment
rn_mean_by_env(env     = vec_env,
               shape   = expr,
               theta   = theta,
               G_theta = G,
               fixed = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in G

# This is (slightly) different from the function evaluated
# using the average values of the parameters in theta 
theta["cmax"] * exp(
         - exp(theta["rho"] * (vec_env - theta["xopt"]) - 6) -
             theta["sigmagaus"] * (vec_env - theta["xopt"])^2
     )
}
