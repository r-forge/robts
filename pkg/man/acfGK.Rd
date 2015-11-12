\encoding{utf8}
\name{acfGK}
\alias{acfGK}

\title{
Robust Autocorrelation Estimation Based on GK Approach
}
\description{
Robustly estimates the autocorrelation function of a time series based on GK estimators. See Ma and Genton (2000) for details.
}

\usage{
acfGK(timeseries,maxlag,method="Qn",...)
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
	\item{method}{character string indicating which scale estimator should be used. Must be one of 'Qn','Tau' or 'MAD', see details.}
	\item{...}{further arguments passed to the respective internal function. For instance tuning parameter for scale Tau, see details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

This function estimates the autocorrelation function based on the GK approach, which was proposed by Ma and Genton (2000). The acf is estimated for each lag individually. Let X denote the original timeseries and Y the lagged one, then the autocorrelation is estimated by
\deqn{\frac{\hat{\sigma}(X+Y)^2-\hat{\sigma}(X-Y)^2}{\hat{\sigma}(X+Y)^2+\hat{\sigma}(X-Y)^2}.}
If one uses the standard deviation as scale estimator \eqn{\hat{\sigma}}, one roughly gets the usual empirical acf. For robust estimation different scale estimators are available.
  
	
	If \code{method="Qn"} the Qn is used, which was proposed by Rousseeuw and Croux (1993). The estimator is known to be very robust and also rather efficient under normality. See the help of \code{\link{Qn}} for details.
	
	If \code{method="Tau"} the scale Tau estimator is used, which was proposed by Maronna and Zamar (2002). The estimator is also known to be robust and rather efficient under normality. The estimator includes two tuning-parameters, which are set as in the default of \code{\link{scaleTau2}} but can be changed using the ... argument. See the help of \code{\link{scaleTau2}} for more details.
	
	If \code{method="MAD"} one uses the usual median absolute deviation, which is known to be very robust and fast to compute. For more details see the help of \code{\link{mad}}.
}

\references{
Ma, Y. and Genton M. (2000): Highly robust estimation of the autocovariance function, \emph{Journal of time series analysis}, vol. 21, 663--684.

Maronna, R. and Zamar, R. (2002): Robust estimates of location and dispersion for high-dimensional datasets, \emph{Technometrics}, vol. 44, 307--317.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}.

Robust scale estimators: \code{\link{mad}}, \code{\link{Qn}}, \code{\link{scaleTau2}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfGK(tss,7)
}
