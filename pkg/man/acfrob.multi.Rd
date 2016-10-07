\encoding{utf8}
\name{acfrob.multi}
\alias{acfrob.multi}

\title{
Robust Autocorrelation Estimation Based on Correlation Matrices 
}

\description{
Robustly estimates the autocorrelation function of a time series based on robust estimators of the covariance matrix. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfmulti"}.
}

\usage{
acfrob.multi(x, lag.max, multi.method=c("weightedMCD", "rawMCD", "Stahel-Donoho", "S",
         "reweight", "Tyler", "M", "sscor"), ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{multi.method}{character string indicating the multivariate covariance estimator to be used. Must be one of \code{'weightedMCD'} (the default), \code{'rawMCD'}, \code{'Stahel-Donoho'}, \code{'S'}, \code{'M'}, \code{'reweight'} and \code{'Tyler'}, see Details.}
	\item{...}{additional arguments passed to the respective multivariate covariance estimator, see Details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{
This function estimates the autocorrelation function based on multivariate correlation estimators. This method computes the acf en bloc and is therefore very time consuming if \code{lag.max} is large (e.g. greater 15). The approach works very well in case of outliers occuring in blocks, see \enc{Dürre}{Duerre} et al. (2015) for more details. There are many robust correlation estimators available.

	If \code{multi.method = "rawMCD"}, one applies the MCD proposed by Rousseeuw (1985). This estimator is known to be very robust but also inefficient under normality. Therefore one often uses reweighted versions. See the help page of \code{\link{covMcd}} for more details.

	If \code{multi.method = "weightedMCD"}, one applies the weighted MCD, which is basically the empirical covariance of observations which are not identified as outliers. The identification is based on Mahalanobis distances calculated from the raw MCD, see \code{\link{covMcd}} for more details. For acf estimation the weighted MCD is preferable to the unweighted one, see \enc{Dürre}{Duerre} et al. (2014).

	If \code{multi.method = "reweight"}, one applies the efficient weighted MCD which was proposed by Gervini (2003). See the help page of \code{\link{Corefw}} for more details.

	If \code{multi.method = "Stahel-Donoho"}, one applies the Stahel-Donoho estimator. See the help page of \code{\link{CovSde}} for more details. Note that this estimator gets really time consuming if \code{lag.max} is large!

	If \code{multi.method = "S"}, one applies an S estimator of scatter. See the help page of \code{\link{CovSest}} for more details.

	If \code{multi.method = "M"}, one applies an M estimator of scatter. See the help page of \code{\link{mvhuberM}} for more details.

	If \code{multi.method = "Tyler"}, one applies the Tyler shape matrix, which is a special M estimator of scatter. See the help page of \code{\link{tyler.shape}} for more details.
}  

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

Gervini, D. (2003): A robust and efficient adaptive reweighted estimator of multivariate location and scatter, \emph{Journal of Multivariate Analysis}, vol. 84, 116--144.

Rousseeuw, P. (1985): Multivariate estimation with high breakdown point, \emph{Mathematical Statistics and Applications Vol. B}, 283--297.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.regression}}, \code{\link{acfrob.trim}}.

Robust multivariate scatter estimators: \code{\link{covMcd}}, \code{\link{Corefw}}, \code{\link{CovSde}}, \code{\link{CovSest}}, \code{\link{mvhuberM}}, \code{\link{tyler.shape}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "multi")
}
