\encoding{utf8}
\name{reweightedQn}
\alias{reweightedQn}

\title{
Reweighted Variant of the Robust Scale Estimator Qn
}
\description{
Computes a reweighted variant of the robust scale estimator Qn.
}

\usage{
reweightedQn(x, boundq = 0.999)
}

\arguments{
	\item{x}{numeric vector of observations.}
	\item{boundq}{numeric value specifying the quantile of the chi square distribution for the reweighting step, see Details.}
}

\details{
In a first step the ordinary Qn estimator is computed. In a second step the empirical standard deviation is computed using all observation which are lower than the \code{boundq} quantile of a chi square distribution with one degree of freedom. This additional reweighting step improves the efficiency as compared to the ordinary Qn estimator, which was introduced by Ma and Genton (2001).
}

\value{
Numeric value with the estimated scale.
}

\references{
Ma, Y. and Genton M. (2001): Highly Robust Estimation of Dispersion Matrices, \emph{Journal of Multivariate Analysis}, vol. 78, 11-36, \doi{10.1006/jmva.2000.1942}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
\code{\link{Qn}} for the non-weighted Qn estimator.
}
