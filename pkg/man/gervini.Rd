\encoding{utf8}
\name{gervini}
\alias{gervini}

\title{
Robust estimation of the correlation matrix
}
\description{
Calculates the correlation matrix robustly, possibly with a fully efficient reweighting step, see \enc{Dürre}{Duerre} et al. (2014) for details.
}
\usage{
gervini(data, boundq = 0.975, startestimator = c("wMCD", "rMCD", "S"))
}

\arguments{
	\item{data}{Matrix with the observations in rows.}
	\item{boundq}{Quantile-bound for comparism of theoretical and empirical distribution function.}
	\item{startestimator}{Start estimator to be used: Minimum covariance determinant (\code{rMCD} for raw version, \code{wMCD} for a reweighting step) or \code{S}-estimator.}
}

\value{
Correlation matrix.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.
}
\author{
Alexander \enc{Dürre}{Duerre}
}

\examples{
set.seed(1066)
x <- rnorm(100)
set.seed(1485)
X <- cbind(x, x + rnorm(100), x^2)
gervini(X)
}