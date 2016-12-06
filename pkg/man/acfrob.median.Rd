\encoding{utf8}
\name{acfrob.median}
\alias{acfrob.median}

\title{
Robust Autocorrelation Estimation Based on Median Correlation
}

\description{
Robustly estimates the autocorrelation function of a time series based on the median correlation. See Chakhchoukh (2010) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfmedian"}.
}

\usage{
acfrob.median(x, lag.max, biascorr = TRUE)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{biascorr}{logical indicating whether a bias correction of the result should be done. Note that this bias correction is only valid for Gaussian processes.}
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. The ARE is currently not available for this estimation approach and is therefore \code{NA}.}
}

\details{
This function estimates the autocorrelation function based on the median correlation, which is proposed by Chakhchoukh (2010). The acf is estimated for each lag individually. Let \eqn{X} denote the original time series and \eqn{Y} the lagged one, then the autocorrelation is estimated by
\deqn{median(XY)/median(X^2).}
The resulting value does not measure the usual autocorrelation. Therefore a nonlinear transformation based on a Monte Carlo simulation is applied, which is only valid for Gaussian processes.
}  

\references{
Chakhchoukh, Y. (2010): A new robust estimation method for ARMA models, \emph{IEEE Transactions on Signal Processing}, vol. 58, 3512--3522.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.bireg}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "median")
}
