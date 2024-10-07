\name{rn_cs_gen}
\alias{rn_cs_gen}

\title{Decompose the \eqn{\mathbf{G}}-matrix from a character-state model}
\description{
This function takes the \eqn{\mathbf{G}}-matrix estimated from a character-state model and returns the values for \eqn{V_{\text{Add}}},  \eqn{V_{\text{A}}}, \eqn{V_{\text{A}\times\text{E}}} and \eqn{n_{\text{eff}}}
}
\usage{
rn_cs_gen(G_cs, wt = NULL)
}
\arguments{
  \item{G_cs}{(Additive) genetic variance-covariance matrix estimated from a character-state model (i.e. where environments are treated as a categorical variable). (numerical matrix)}
  \item{wt}{Weights to apply to the different environments, e.g. reflecting their frequencies in the wild. The weights must non-negative, and at least one must be non-zero. The vector \code{wt} must be the same length as the rows and columns of \code{G_cs}. By default, no weighting is applied. (numeric)}
} 

\details{
\eqn{V_{\text{Add}}} is the (weighted) average of the diagonal elements of \code{G_cs}, \eqn{V_{\text{A}}} is the (weighted) average of all the elements of \code{G_cs} and \eqn{V_{\text{A}\times\text{E}}} is the difference between \eqn{V_{\text{Add}}} and  \eqn{V_{\text{A}}}. Finally, the efficient number of dimensions \eqn{n_{\text{eff}}} is the ratio of the sum of the eigen values of \code{G_cs} over its maximum eigen  value. Note that \eqn{n_{\text{eff}}} is returned for information, but is expected to be biased in practice due to an over-estimation of the maximum eigen value.
}
\value{This function yields \eqn{V_{\text{Add}}},  \eqn{V_{\text{A}}},  \eqn{V_{\text{A}\times\text{E}}} and \eqn{n_{\text{eff}}} as a one-row data.frame (data.frame, all numeric).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{rn_vgen}}, \code{\link{rn_gen_decomp}}
}
\examples{
G <- diag(10)

rn_cs_gen(G_cs = G)
}
