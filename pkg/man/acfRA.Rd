\encoding{utf8}
\name{acfRA}
\alias{acfRA}

\title{
Robust autocorrelation estimation based on residual autocorrelation estimators 
}
\description{
Robustly estimates the autocorrelation function of a time series based on a robustly transformed timeseries. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfRA(timeseries,maxlag,Psi="Huber",meanvalue=median,scattervalue=mad,...)
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
	\item{Psi}{character string indicating the used Psi function. Must be either '"Huber"' or '"Tukey"', see details.}
	\item{meanvalue}{function which estimates the location of the timeseries.}
	\item{scattervalue}{function which estimates the scale of the timeseries.}
	\item{...}{tuning parameters for Huber or Tukey Psi function, see details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

The function estimates the residual autocovariance, which is the usual acf of the robustly transformed timeseries. Using an estimator for location and scale which can be set using the arguments \code{meanvalue} and \code{scattervalue} the timeseries is transformed by applying a psi function, only Huber and Tukey are possible. The tuning parameter for the Huber function is \code{k=1.37} and for Tukey \code{k=4.68} but both can be changed using the ... argument. For the meaning of the parameters, see \enc{Dürre}{Duerre} et al. (2014).

There is a simulation based consistency correction implemented for Gaussian timeseries and the preset tuning parameters.
}  

\references{

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfRA(tss,7)
}
