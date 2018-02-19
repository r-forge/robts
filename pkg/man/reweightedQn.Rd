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
In a first step the ordinary Qn estimator is computed. In a second step the empirical standard deviation is computed using all observations whose squared residuals based median and Qn are smaller than the \code{boundq} quantile of a chi square distribution with one degree of freedom. This additional reweighting step improves the efficiency as compared to the Qn estimator, which was introduced by Rousseeuw and Croux (1993).
}

\value{
Numeric value with the estimated scale.
}

\references{
Rousseeuw, P.J., Croux, C. (1993): Alternatives to the median absolute deviation, \emph{Journal of the American Statistical Association}, vol. 88, 1273--1283.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
\code{\link{Qn}} for the non-reweighted Qn estimator.
}
