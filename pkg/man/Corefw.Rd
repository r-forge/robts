\encoding{utf8}
\name{Corefw}
\alias{Corefw}

\title{
Fully Efficient Robust Correlation Estimation
}
\description{
Calculates a robust and asymptotically fully efficient correlation matrix, see Gervini (2003) for details.
}
\usage{
Corefw(data, boundq = 0.975, startestimator = "wMCD")
}

\arguments{
	\item{data}{matrix with the observations in rows.}
	\item{boundq}{quantile bound for comparison of theoretical and empirical distribution function, see details.}
	\item{startestimator}{start estimator to be used: \code{"rMCD"}, \code{"wMCD"} and \code{"S"}-estimator are possible, see details.}
}

\value{
Numeric correlation matrix.
}

\details{
This implementation calculates the (asymptotically) fully efficient scatter estimator proposed by Gervini (2003). Based on a an initial scatter estimator which can be determined by \code{startestimator}, residuals and their cumulative distribution function are calculated. This empirical distribution is compared with the theoretical one (a chi square distribution). Only if the empirical distribution function lays under the theoretic one, from the \code{boundq} quantile onwards, observations are identified as outliers, the number depends on the distance between the distribution functions. Based on all observations which are not marked as outliers, the usual correlation is returned.

This procedure has asymptotically the same efficiency as the usual empirical correlation and retains the breakdown point of the initial scale estimator.
}

\references{
Gervini, D. (2003): A robust and efficient adaptive reweighted estimator of multivariate location and scatter, \emph{Journal of multivariate analysis}, vol 84, 116--144.
}
\author{
Alexander \enc{Dürre}{Duerre}
}

\examples{
set.seed(1066)
x <- rnorm(100)
set.seed(1485)
X <- cbind(x, x + rnorm(100), x^2)
Corefw(X)
}
