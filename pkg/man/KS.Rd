\encoding{utf8}
\name{KS}
\alias{pKS}
\alias{qKS}

\title{
The Kolmogorov Smirnov Distribution
}
\description{
Computes the distribution function and the quantiles of the two-sided Kolmogorov Smirnov (KS) test.
}

\usage{
pKS(q)
qKS(p)
}

\arguments{
	\item{q}{numeric vector of quantiles.}
	\item{p}{numeric value with one probability (not a vector of probabilities).}
}

\details{
The distribution equals that of the supremum of the absolute values of a standard Brownian bridge and therefore often occurs in change point testing. The function is only a wrapper for the C code used in the function \code{\link{ks.test}}.
}

\value{
Function \code{pKS} returns a numeric vector with values of the distribution function at the given quantiles.

Function \code{qKS} returns a numeric value of the desired quantile.
}

\author{
R Core Team and Alexander \enc{Dürre}{Duerre}
}


\seealso{
\code{\link{pKS}} for the distribution function.
}

\examples{
curve(pKS(x), from = 0, to = 2)

qKS(0.95)
}
