\name{rn_gen_decomp}
\alias{rn_gen_decomp}

\title{Compute the decomposition of the additive genetic variance and their \eqn{\gamma}- and \eqn{\iota}-decomposition}
\description{
This function calculates the total additive genetic variance of the reaction norm (\eqn{V_{\text{Add}}}, and its \eqn{\gamma}-decomposition), the marginal additive genetic variance of the trait (\eqn{V_{\text{A}}}) and the additive genetic variance of plasticity (\eqn{V_{\text{A}\times\text{E}}}, and its \eqn{\iota}-decomposition).
}
\usage{
rn_gen_decomp(theta, G_theta, env = NULL, shape = NULL, X = NULL, fixed = NULL,
              wt_env = NULL, compute_gamma = TRUE, compute_iota  = TRUE, width = 10)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Expression providing the shape of the reaction where \code{x} is the environment. For example: \code{expression(a + b * x + c * x^2)}.}
  \item{X}{If the model used was linear in the parameters, the design matrix X of the model (numeric, incompatible with the arguments \code{env} and \code{shape}).}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{wt_env}{Weights to apply to the \code{env} vector values, providing an information regarding their relative probability in the biological context. The weights must non-negative, and at least one must non-zero. The vector \code{wt_env} must be the same length as \code{env}. By default, no weighting is applied. (numeric)}
  \item{compute_gamma}{Should the \eqn{\gamma}-decomposition be performed? Default is TRUE. (boolean)}
  \item{compute_iota}{Should the \eqn{\iota}-decomposition be performed? Default is TRUE. (boolean)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
} 

\details{
The variance \eqn{V_{\text{Add}}} is the total additive genetic variance in the reaction norm. It can be decomposed according to each of the parameters in \code{G_theta} (and their pairwise covariance), using the \eqn{\gamma}-decomposition, if the \code{compute_gamma} is activated. It can also be decomposed into the marginal additive genetic variance of the trait (\eqn{V_{\text{A}}}, the additive genetic variance of the trait that would be computed if plasticity were ignored), and the additive genetic variance of plasticity itself (\eqn{V_{\text{A}\times\text{E}}}, i.e. the genetic-by-environment interaction variance). This last component can also be decomposed according to each of the parameters in \code{G_theta} (and their pairwise covariance), using the \eqn{\iota}-decomposition (for interaction).

It is very important that the parameters are in the same order in the argument of the \code{shape} function, in \code{theta} and in \code{G_theta}.
}
\value{This function yields the values for \eqn{V_{\text{Add}}}, \eqn{V_{\text{A}}} and \eqn{V_{\text{A}\times\text{E}}}, as well as the \eqn{\gamma}- and \eqn{\iota}-decomposition if requested, formatted as a 1-row data.frame (data.frame, all numeric).}

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
out <- 
  rn_gen_decomp(theta   = theta,
                G_theta = G,
                env     = vec_env,
                shape   = expr,
                fixed   = c(3, 4))
# Note that fixed is set for the third and forth parameters than are not in G

# Checking that V_Add = V_A + V_AxE
out$V_Add == out$V_A + out$V_AxE
# Checking that the sum of \eqn{\gamma}'s is 1
out$Gamma_cmax + out$Gamma_xopt + out$Gamma_cmax_xopt
# Checking that the sum of \eqn{\iota}'s is 1
out$Iota_cmax + out$Iota_xopt + out$Iota_cmax_xopt

# Applying some weighting to the environmental values according to a normal distribution
rn_gen_decomp(theta   = theta,
              G_theta = G,
              env     = vec_env,
              shape   = expr,
              fixed   = c(3, 4),
              wt_env  = dnorm(vec_env))
              
# If a polynomial was used, it is possible to use the linear modeling rather having
# to compute integrals
theta <- c(a = 1.5, b = 0.5, c = -0.5)
X     <- cbind(1, vec_env, (vec_env - mean(vec_env))^2)
G     <- 0.1 * diag(3)
rn_gen_decomp(theta   = theta,
              G_theta = G,
              X        = X)
}
