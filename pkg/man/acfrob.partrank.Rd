\encoding{utf8}
\name{acfrob.partrank}
\alias{acfrob.partrank}

\title{
Robust Autocorrelation Estimation Based on Partial Autocorrelations 
}

\description{
Robustly estimates the autocorrelation function of a time series based on a robust estimation of the partial autocorrelation function. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfpartrank"}.
}

\usage{
acfrob.partrank(x, lag.max, cor.method = c("spearman", "kendall", "quadrant", "gaussian",
            "masarotto"), biascorr = TRUE, partial = FALSE)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{cor.method}{character string indicating the correlation estimator to be used. Must be one of \code{'spearman'} (the default), \code{'kendall'}, \code{'quadrant'}, \code{'gaussian'}, \code{'masarotto'}, see Details.}
	\item{biascorr}{logical indicating whether a bias correction of the result should be done. Note that this bias correction is only valid for Gaussian processes. There is no bias correction available for \code{partial.method = 'masarotto'}.}
	\item{partial}{logical indicating whether the function should return the partial acf instead of the acf.}
}

\details{
This function estimates the autocorrelation function based on partial autocorrelations. This procedure was first proposed by Masarotto (1987) and has the advantage that it produces positive definite and therefore valid acf estimations. To estimate the pacf, the estimation procedure described in \enc{Möttönen}{Moettoennen} et al. (1999) is used. It depends on successive correlation estimations of forward and backward residuals.

There are several possible rank estimators available. The argument \code{cor.method} specifies which of the following approaches is to be used, where \code{cor.method = "spearman"} is the default.

\describe{
\item{\code{"spearman"}}{One uses the well known Spearman's rho with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Spearman's rho.}

\item{\code{"kendall"}}{One uses the well known Kendall's tau with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Kendall's tau. Note that Kendall's tau can get rather time consuming for large time series (e.g. greater than 10000).}

\item{\code{"quadrant"}}{One uses the quadrant correlation with a consistency correction under normality. See Croux and Dehon (2010) for this correction and more information about the estimator. Note that this estimator has a rather low efficiency under normality and should be therefore only applied for time series with a high amount of outliers.}

\item{\code{"gaussian"}}{One uses the Gaussian rank correlation which is the usual correlation applied to transformed values. Denote \eqn{R(i)} the rank of \eqn{X(i)} and \eqn{\Phi}{Phi} the distribution function of the normal distribution, then the transformation is defined as
\deqn{\Phi^{-1}(R(i)/(n+1)).}
This estimator has the same asymptotic efficiency as the usual empirical correlation under normality, but is robust to some degree.}

\item{\code{"masarotto"}}{One follows the original proposal of Masarotto (1987), which is a kind of M estimator. Starting from initial estimations for the variance by the \code{\link{Qn}} and for the correlation by the GK estimator, (using the Qn) new variances and correlations are iteratively computed until convergence, using weights which are defined by the following psi function
\deqn{\psi(x) = 3/(1+x).}{psi(x) = 3/(1+x)}
Note that this is not really a rank or sign based estimator. See Masarotto (1987) for more details.}
}

The computation of the acf from the pacf is carried out by the algorithm described in \enc{Möttönen}{Moettoennen} et al. (1999).
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated (partial) autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. The ARE is currently not available for this estimation approach and is therefore \code{NA}.}
}

\references{
Croux, C. and Dehon, C. (2010): Influence functions of the Spearman and Kendall correlation measures, \emph{Statistical Methods & Applications}, vol. 19, 497--515, \doi{10.1007/s10260-010-0142-z}.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222, \doi{10.1002/wics.1351}.

Masarotto, D. (2003): Robust identification of autoregressive moving average models, \emph{Applied statistics}, vol. 36, 214--220, \doi{10.2307/2347553}.

\enc{Möttönen}{Moettoennen}, J., Koivunen, V. and Oja, H. (1999): Robust autocovariance estimation based on sign and rank correlation coefficients, In: \emph{Proceedings of the IEEE Signal Processing Workshop on Higher-Order Statistics}, 187--190, \doi{10.1109/HOST.1999.778722}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.bireg}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "partrank")
}

\keyword{Autocorrelation function}
