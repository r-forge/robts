\encoding{utf8}
\name{acfrobfil}
\alias{acfrobfil}

\title{
Robust autocorrelation estimation based on robust filtering
}
\description{
Robustly estimates the autocorrelation function of a time series based on the robust filter algorithm. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfrobfil(timeseries,p,maxlag,robfiltype="filtered",aic=TRUE,psi.l=2,psi.0=3) 
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{p}{numeric determining the (maximal) order of the AR fit.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
	\item{psi.l}{numeric determining the used psi function, see details.}
	\item{psi.0}{numeric determining the used psi function, see details.}
	\item{robfiltype}{characterstring indicating whether the acf is calculated based on the empirical acf using robustly filtered values (\code{robfiltype="filtered"}) or based on an AR fit (\code{robfiltype="ar"}). }
	\item{aic}{logical indicating wether the AR order should be estimated by robust AIC criteria or should be equal to p}
	
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

The function estimates the acf using the robust filter algorithm, which is a kind of robust Kalman filter, as described in chapter 8 of Maronna et al. (2006). The method fits AR models of increasing order by a robust version of the Durbin Levinson algorithm. The AR coefficients are estimated by minimizing a robust scale (\code{\link{scaleTau2}}) of one step ahead residuals of robustly filtered predictors. The filtering process depends on the choosen Psi function and respects the estimated dependence structue of the timeseries. The order of the AR process can be determined by the user by the argument \code{p} if one chooses \code{aic=FALSE} or estimated by a robust aic criteria if \code{aic=TRUE}. In this case \code{p} represents the maximal possible AR order.

There are two possibilities to derive an acf estimation.

If \code{robfiltype="filtered"} one calculates the usual empirical acf from the robustly filtered timeseries/predictors.

If \code{robfiltype="ar"} one calculates the acf by the estimated AR coefficients using the function \code{\link{ARMAacf}}.

Note that the algorithm is very timeconsuming even for small timeseries, especially if one chooses large \code{p}.

More details cn be found in Maronna et al. (2006), \enc{Dürre}{Duerre} et al. (2014) and the help of the subroutine \code{\link{ARfilter}}.
}

\references{

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Maronna, R., Martin, D. and Yohai, V. (2006): \emph{Robust statistics}, Wiley, Chichester.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acftrim}}

The subroutine \code{\link{ARfilter}}

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrobfil(tss,p=5,maxlag=7)
}
