\encoding{utf8}
\name{Corefw}
\alias{Corefw}

\title{
Robust Estimation of the Correlation Matrix
}
\description{
Calculates the correlation matrix robustly, possibly with a fully efficient reweighting step, see \enc{Dürre}{Duerre} et al. (2014) for details.
}
\usage{
Corefw(data, boundq = 0.975, startestimator = "wMCD")
}

\arguments{
	\item{data}{matrix with the observations in rows.}
	\item{boundq}{quantile-bound for comparism of theoretical and empirical distribution function.}
	\item{startestimator}{start estimator to be used: Minimum covariance determinant (\code{"rMCD"} for raw version, default \code{"wMCD"} for a reweighting step) or \code{"S"}-estimator.}
}

\value{
Numeric correlation matrix.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.
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