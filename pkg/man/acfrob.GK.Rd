\encoding{utf8}
\name{acfrob.GK}
\alias{acfrob.GK}

\title{
Robust Autocorrelation Estimation Based on the Gnanadesikan Kettenring Approach
}

\description{
Robustly estimates the autocorrelation function (acf) of a time series based on Gnanadesikan Kettenring (GK) estimators. See Ma and Genton (2000) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfGK"}.
}

\usage{
acfrob.GK(x, lag.max, scalefn = Qn, ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{scalefn}{function which calculates a scale estimator. Its first argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
	\item{...}{further arguments passed to the function \code{scalefn}.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{
This function estimates the autocorrelation function based on the GK approach, which was proposed by Ma and Genton (2000). The acf is estimated for each lag individually. Let \eqn{X} denote the original time series and \eqn{Y} the lagged one, then the autocorrelation is estimated by
\deqn{\frac{\hat{\sigma}(X+Y)^2-\hat{\sigma}(X-Y)^2}{\hat{\sigma}(X+Y)^2+\hat{\sigma}(X-Y)^2}.}
If one uses the empirical standard deviation as the scale estimator \eqn{\hat{\sigma}}, one roughly gets the usual empirical acf. For robust estimation different scale estimators are available.
  
	
	If \code{scalefn = Qn}, the Qn is used, which was proposed by Rousseeuw and Croux (1993). The estimator is known to be very robust and also rather efficient under normality. See the help page of \code{\link{Qn}} in the package \code{robustbase} for details.
	
	If \code{scalefn = scaleTau2}, the \eqn{\tau} estimator of scale is used, which was proposed by Maronna and Zamar (2002). The estimator is also known to be robust and rather efficient under normality. The estimator includes two tuning parameters, which can be changed using the \code{...} argument. See the help page of \code{\link{scaleTau2}} in the package \code{robustbase} for more details.
	
	If \code{scalefn = MAD}, one uses the usual median absolute deviation, which is known to be very robust and fast to compute, but not very efficient under normality. For more details see the help page of \code{\link{mad}} in the package \code{stats}.

	If \code{scalefn = reweightedQn}, one uses an efficiently reweighted Qn. ADD MORE DETAILS HERE OR LINK TO THE RESPECTIVE HELP PAGE.
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

Alternative acf subroutines: \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.regression}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "GK")
}
