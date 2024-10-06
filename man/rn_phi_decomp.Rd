\name{rn_phi_decomp}
\alias{rn_phi_decomp}

\title{Compute the \eqn{\varphi}-decomposition of \eqn{V_{\text{Plas}}}}
\description{
This function calculates, for a linear (in its parameters) model, the plastic variance of the average (over genotypes) reaction norm shape \eqn{V_{\text{Plas}}} and its \eqn{\varphi}-decomposition (depending on the parameters rather than slope/curvature). This is especially useful when the reaction norm model was a polynomial.
}
\usage{
rn_phi_decomp(theta, X, S = NULL, wt_env = NULL, correction = FALSE, v_plas = NA)
}
\arguments{
  \item{theta}{Average parameters of the shape function. It must be a named vector, with the names corresponding to the parameters in the \code{shape} expression. (numeric)}
  \item{X}{The design matrix X of the reaction norm model. For a polynomial, X contains all the exponents of the environment. (numerical matrix)}
  \item{S}{The error variance-covariance matrix of the estimated fixed effects \code{theta}. Not mandatory, but it is strongly advised to correct for the error variance in the parameter estimates. The function throws a warning if it is \code{NULL} (numerical matrix).}
  \item{wt_env}{Weights to apply to the \code{env} vector values, providing an information regarding their relative probability in the biological context. The weights must non-negative, and at least one must non-zero. The vector \code{wt_env} must be the same length as \code{env}. By default, no weighting is applied. (numeric)}
  \item{correction}{Should Bessel's correction (dividing the sum-of-squares by N-1 rather than N) be used (TRUE) or not (FALSE, default). The default is FALSE, because it is likely that other components such as the total phenotypic variance is computed over the number of individuals (generally large number) rather than the number of environments (generally small number). The best is to manually use Bessel's correction over the proper number of data points. (boolean)}
  \item{v_plas}{If already computed, an estimate for \eqn{V_{\text{Plas}}}. If NA, then \eqn{V_{\text{Plas}}} will be computed by the function. (numeric)}
} 

\details{
The variance \eqn{V_{\text{Plas}}} is the variance arising from the shape of the reaction norm after avering over the genetic variance, i.e. the purely environmental part of the variance of plasticity in the reaction norm. Its \eqn{\phi}-decomposition is based on the parameters of the linear model \eqn{\mathbf{X}\theta}, rather than on slope and curvature, as for the \eqn{\pi}-decomposition.

It is very important that the parameters are in the same order in \code{theta} (which, again, must be named) and in \code{G_theta}.
}
\value{This function yields the \eqn{\pi}-decomposition of \eqn{V_{\text{Plas}}} as a one-row data.frame. If \eqn{V_{\text{Plas}}} was provided as an argument (e.g. as computed from a separate character-state model), it also outputs \eqn{M^2}, the ratio between \eqn{V_{\text{Plas}}} as computed from the linear model parameters and the provided value for \eqn{V_{\text{Plas}}} (data.frame, all numeric).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{rn_pi_decomp}}, \code{\link{rn_vplas}}
}
\examples{
# Parameters
vec_env <- rnorm(20)
theta <- c(a = 1.5, b = 0.5, c = -0.5)
X     <- cbind(1, vec_env, (vec_env - mean(vec_env))^2)
S     <- (10^-4) * diag(3)

# Computing the phi-decomposition
rn_phi_decomp(theta = theta, X = X, S = S)
}
