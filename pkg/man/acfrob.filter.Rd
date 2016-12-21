\encoding{utf8}
\name{acfrob.filter}
\alias{acfrob.filter}

\title{
Robust Autocorrelation Estimation Based on Robust Filtering
}

\description{
Robustly estimates the autocorrelation function of a time series based on the robust filter algorithm. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfrobfil"}.
}

\usage{
acfrob.filter(x, lag.max, order.max = lag.max, robfil.method = c("filtered", "ar"),
              aic = TRUE, aicpenalty=function(p) {2*p}, psi.l = 2, psi.0 = 3,
              partial = FALSE)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{order.max}{integer value determining the (maximal) order of the AR fit.}
	\item{robfil.method}{character string indicating whether the acf is calculated based on the empirical acf using robustly filtered values (\code{robfil.method = "filtered"}) or based on an AR fit (\code{robfil.method = "ar"}). }
	\item{aic}{logical indicating whether the AR order should be estimated by robust AIC criteria considering orders up to \code{order.max} or should be equal to \code{order.max}.}
	\item{aicpenalty}{function of the model order, indicating the penalty for a larger model.}
	\item{psi.l}{numeric value determining the used psi function, see Details.}
	\item{psi.0}{numeric value determining the used psi function, see Details.}
	\item{partial}{logical indicating whether the function should return the partial acf instead of the acf.}
}

\details{
The function estimates the acf using the robust filter algorithm, which is a kind of robust Kalman filter, as described in chapter 8 of Maronna et al. (2006). The method fits AR models of increasing order by a robust version of the Durbin Levinson algorithm. The AR coefficients are estimated by minimizing a robust scale (\code{\link{scaleTau2}}) of one step ahead residuals of robustly filtered predictors. The filtering process depends on the psi function and respects the estimated dependence structure of the time series. The order of the AR process can be determined by the argument \code{order.max} if one chooses \code{aic = FALSE} or estimated by a robust AIC criterion if \code{aic = TRUE}. In this case \code{order.max} represents the maximal possible AR order.

There are two possibilities to derive an acf estimation:

If \code{robfil.method = "filtered"}, one calculates the usual empirical acf from the robustly filtered time series/predictors.

If \code{robfil.method = "ar"}, one calculates the acf by the estimated AR coefficients using the function \code{\link{ARMAacf}}.

A psi function which is two times continuous differentiable is used here, see \enc{Dürre}{Duerre} et al. (2015) for a formal definition. The form depends on the tuning parameters \code{psi.l} and \code{psi.0}, see the help page of \code{\link{ARfilter}} for more details.

More details can be found in Maronna et al. (2006), \enc{Dürre}{Duerre} et al. ( 2015) and the help page of the subroutine \code{\link{ARfilter}}.
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated (partial) autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. The ARE is currently not available for this estimation approach and is therefore \code{NA}.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222, \doi{10.1002/wics.1351}.

Maronna, R., Martin, D. and Yohai, V. (2006): \emph{Robust statistics}, Wiley, Chichester, \doi{10.1002/0470010940}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.bireg}}, \code{\link{acfrob.trim}}.

The subroutine \code{\link{ARfilter}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "filter", order.max = 5)
}

\keyword{Autocorrelation function}
