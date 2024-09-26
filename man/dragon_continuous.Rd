\docType{data}
\name{dragon_continuous}
\alias{dragon_continuous}

\title{A simulated dataset on the reaction norms of dragons in a continuous environment}
\format{
A data frame with 1000 rows and 5 variables:
\describe{
\item{Ind}{The individual on which the measures were performed at the different temperatures.}
\item{Family}{The family to which each individual belongs to, i.e.\ which mother they are the offspring off. Offspring from the same mothers are considered as half-sibs.}
\item{Temp}{The value of the recorded temperature at which the experiment was performed in the field. For convenience, the value was standardised so that the average is 0.}
\item{Aggressiveness}{Aggressiveness of the dragon, measured using a very precise, continuous index, when presented with an armored knight to combat in the field.}
\item{Performance}{Performance (sprint speed to rush toward a dummy princess) of each individual dragon evaluated directly on the field.}
}
}
\usage{
dragon_continuous
}
\description{
A simulated dataset on the reaction norm of aggressiveness and thermal performance of dragons measured directly in a wild population, with heterogenous temperatures. Offspring from the same mothers are considered as half-sibs.
}
\keyword{datasets}
