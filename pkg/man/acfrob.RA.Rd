\encoding{utf8}
\name{acfrob.RA}
\alias{acfrob.RA}

\title{
Robust Autocorrelation Estimation Based on Residual Autocorrelation
}

\description{
Robustly estimates the autocorrelation function of a time series based on a robustly transformed time series. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfRA"}.
}

\usage{
acfrob.RA(x, lag.max, psi = c("huber", "bisquare"), k, locfn = median, scalefn = mad,
      biascorr = TRUE, ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{psi}{character string indicating the psi function to be used. Must be either \code{'huber'} for Huber's function or \code{'bisquare'} for Tukey's function, see Details.}
  \item{k}{numeric tuning parameter determining the psi function. If missing, the function chooses \code{k = 1.37} (for Huber's function) or \code{k = 4.68} (for Tukey's function).}
	\item{locfn}{function which calculates a location estimator. Its argument must be a vector of the data.}
	\item{scalefn}{function which calculates a scale estimator. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
	\item{biascorr}{logical indicating whether a bias correction of the result should be done. Note that this bias correction is only valid for Gaussian processes.}
	\item{...}{further arguments passed to the function \code{scalefn}.}
}

\details{

The function estimates the residual autocovariance, which is the usual acf of the robustly transformed time series. Using estimators for location and scale, which can be chosen by using the arguments \code{locfn} and \code{scalefn}, the time series is transformed by applying a psi function (currently only Huber and Tukey are possible). The tuning parameter for Huber's function is \code{k = 1.37} and for Tukey's function \code{k = 4.68} but both can be changed using the \code{...} argument. For the meaning of the tuning parameters, see e.g. \enc{Dürre}{Duerre} et al. (2015).

There is a simulation based consistency correction implemented for Gaussian time series and the preset tuning parameters.
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. The ARE is currently not available for this estimation approach and is therefore \code{NA}.}
}

\references{

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222, \doi{10.1002/wics.1351}.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.bireg}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "RA")
}

\keyword{Autocorrelation function}
