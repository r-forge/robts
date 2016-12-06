\encoding{utf8}
\name{acfrob.trim}
\alias{acfrob.trim}

\title{
Robust Autocorrelation Estimation Based on Trimming
}

\description{
Robustly estimates the autocorrelation function of a time series based on the trimmed time series. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acftrim"}.
}

\usage{
acfrob.trim(x, lag.max, trim = 0.1, biascorr = TRUE)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{trim}{numeric indicating the amount of trimmed values, see Details.}
	\item{biascorr}{logical indicating whether a bias correction of the result should be done. Note that this bias correction is only valid for Gaussian processes.}
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. The ARE is currently not available for this estimation approach and is therefore \code{NA}.}
}

\details{
The function basically first trims the smallest and largest observations of the time series and then calculates the usual sample acf. This procedure was proposed by Chan and Wei (1992).

\code{trim} determines the quantiles (\eqn{\alpha} and \eqn{1-\alpha}) at which the data is trimmed. Therefore \eqn{2\alpha} of the data is trimmed, see \enc{Dürre}{Duerre} et al. (2015) for details.

Since the estimator does not measure the usual autocorrelation, a nonlinear transformation based on simulations is applied, which guarantees consistency under Gaussian processes.
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

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.bireg}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "trim")
}
