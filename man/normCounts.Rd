% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normCounts.R
\name{normCounts}
\alias{normCounts}
\title{Normalise a counts matrix to the median library size}
\usage{
normCounts(x, log = FALSE, prior.count = 0.5, lib.size = NULL)
}
\arguments{
\item{x}{a \code{DGEList} object or matrix of counts.}

\item{log}{logical, indicates whether the output should be on the log2 scale
or counts scale. Default is FALSE.}

\item{prior.count}{The prior count to add if the data is log2 normalised.
Default is a small count of 0.5.}

\item{lib.size}{a vector of library sizes to be used during the normalisation
step. Default is NULL and will be computed from the counts matrix.}
}
\value{
a matrix of normalised counts
}
\description{
This function takes a \code{DGEList} object or matrix of counts and
normalises the counts to the median library size. This puts the normalised
counts on a similar scale to the original counts.
}
\details{
If the input is a DGEList object, the normalisation factors in
\code{norm.factors} are taken into account in the normalisation. The prior
counts are added proportionally to the library size
}
\examples{
# Simulate some data from a negative binomial distribution with mean equal
# to 100 and dispersion set to 1. Simulate 1000 genes and 6 samples.
y <- matrix(rnbinom(6000, mu = 100, size = 1), ncol = 6)

# Normalise the counts
norm.y <- normCounts(y)

# Return log2 normalised counts
lnorm.y <- normCounts(y, log=TRUE)

# Return log2 normalised counts with prior.count = 2
lnorm.y2 <- normCounts(y, log=TRUE, prior.count=2)

par(mfrow=c(1,2))
boxplot(norm.y, main="Normalised counts")
boxplot(lnorm.y, main="Log2-normalised counts")

}
\author{
Belinda Phipson
}
