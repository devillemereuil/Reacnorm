\docType{data}
\name{dragon_continuous}
\alias{dragon_continuous}

\title{A simulated dataset on the reaction norm of aggressiveness and thermal performance of dragons measured directly in a wild population, with heterogenous temperatures}
\format{
A data frame with 1000 rows and 5 variables:
\describe{
\item{Ind}{The individual on which the measures were performed at the different temperatures.}
\item{Family}{The family to which each individual belongs to, i.e.\ which mother they are the offspring off. Offspring from the same mothers are considered as half-sibs.}
\item{Env}{The value of the experimental temperature at which the experiment was performed. For convenience, the value was standardised so that the average is 0.}
\item{Aggressiveness}{Aggressiveness of the dragon, measured using a very precise, continuous index, when presented with an armored knight to combat.}
\item{Performance}{Performance (sprint speed to rush toward a dummy princess) of each individual dragon evaluated at each of the experimental temperatures.}
}
}
\usage{
dragon_continuous
}
\description{
A simulated dataset on the reaction norm of aggressiveness and thermal performance of dragons measured directly in a wild population, with heterogenous temperatures. Offspring from the same mothers are considered as half-sibs.
}
\keyword{datasets}
