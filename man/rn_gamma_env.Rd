\name{rn_gamma_env}
\alias{rn_gamma_env}

\title{Compute the \eqn{\gamma}-decomposition of additive genetic variance per environment}
\description{
This function computes the \eqn{\gamma}-decomposition of the additive genetic variance for each environment.
}
\usage{
rn_gamma_env(theta, G_theta, env = NULL, shape = NULL, X = NULL,
             fixed = NULL, width = 10)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{X}{If the model used was linear in the parameters, the design matrix X of the model (numeric, incompatible with the arguments \code{env} and \code{shape}).}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
This function provides the \eqn{\gamma}-decomposition (according to each parameter of the reaction norm and their covariances) of the additive genetic variance in each environment.

It is very important that the parameters are in the same order in the argument of the \code{shape} function, in \code{theta} and in \code{G_theta}.

If \code{X} is provided the second column is returned as the value for the "environment" (based on the assumption that a polynomial function was used).
}
\value{This function yields a data.frame containing the input values for the environment, the additive genetic variance in each of them, with its corresponding \eqn{\gamma}-decomposition (data.frame, all numeric).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{rn_vgen}}, \code{\link{rn_gen_decomp}}
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
rn_gamma_env(env     = vec_env,
              shape   = expr,
              theta   = theta,
              G_theta = G,
              fixed   = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in G

# If a polynomial was used, it is possible to use the linear modeling rather having
# to compute integrals
theta <- c(a = 1.5, b = 0.5, c = -0.5)
X     <- cbind(1, vec_env, (vec_env - mean(vec_env))^2)
G     <- 0.1 * diag(3)
rn_gamma_env(theta   = theta,
             G_theta = G,
             X       = X)
}
