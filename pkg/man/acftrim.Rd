\encoding{utf8}
\name{acftrim}
\alias{acftrim}

\title{
Robust Autocorrelation Estimation Based on Trimming
}

\description{
Robustly estimates the autocorrelation function of a time series based on the trimmed time series. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acftrim"}.
}

\usage{
acftrim(x, lag.max, trim=0.1, biascorr = TRUE)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{trim}{numeric indicating the amount of trimmed values, see Details.}
	\item{biascorr}{logical indicating whether a bias correction of the result should be done. Note that this bias correction is only valid for Gaussian processes.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{
The function basically first trims the smallest and largest observations of the time series and then calculates the usual sample acf. This procedure was proposed by Chan and Wei (1992).

\code{trim} determines the quantiles (\eqn{\alpha} and \eqn{1-\alpha}) at which the data is trimmed. Therefore \eqn{2\alpha} of the data is trimmed, see \enc{Dürre}{Duerre} et al. (2015) for details.

Since the estimator does not measure the usual autocorrelation, a nonlinear transformation based on simulations is applied, which guarantees consistency under gaussian processes.
}

\references{
Chan, W., and Wei, W. (1992): A comparison of some estimators of time series autocorrelations, \emph{Computational statistics & data analysis}, vol 14, .149--163.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrobfil}}, \code{\link{acfrank}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "acftrim")
}
