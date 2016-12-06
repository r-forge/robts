\encoding{utf8}
\name{acfrob.rank}
\alias{acfrob.rank}

\title{
Robust Autocorrelation Estimation Based on Signs and Ranks 
}

\description{
Robustly estimates the autocorrelation function of a time series based nonparametric estimators using ranks and signs. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfrank"}.
}

\usage{
acfrob.rank(x, lag.max, cor.method=c("gaussian", "spearman", "kendall", "quadrant",
        "masarotto"), biascorr = TRUE)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{cor.method}{character string indicating the used correlation estimator. Must be one of \code{'gaussian'} (the default), \code{'spearman'}, \code{'kendall'}, \code{'quadrant'} and \code{'masarotto'}, see Details.}
	\item{biascorr}{logical indicating whether a bias correction of the result should be done. Note that this bias correction is only valid for Gaussian processes. There is no bias correction available for \code{partial.method = 'masarotto'}.}
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution.}
}

\details{

The function estimates the acf using rank correlation estimators, which estimate the autocorrelations individually for every lag. The approach works also well for count time series, see Fried et al. (2014).

There are several possible rank estimators available. The argument \code{cor.method} specifies which of the following approaches is to be used, where \code{cor.method = "gaussian"} is the default.

\describe{
\item{\code{"gaussian"}}{One uses the Gaussian rank correlation which is the usual empirical correlation applied to transformed values. Denote \eqn{R(i)} the rank of \eqn{X(i)} and \eqn{\Phi} the distribution function of the normal distribution, then the transformation is defined as
\deqn{\Phi^{-1}(R(i)/(n+1)).}
This estimator has the same asymptotic efficiency as the usual empirical correlation under normality, but is robust to some degree.}

\item{\code{"spearman"}}{One uses the well known Spearman's rho with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Spearman's rho.}

\item{\code{"kendall"}}{One uses the well known Kendall's tau with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Kendall's tau. Note that Kendall's tau can get rather time consuming for large time series (e.g. greater than 10000).}

\item{\code{"quadrant"}}{One uses the quadrant correlation with a consistency correction under normality. See Croux and Dehon (2010) for this correction and more information about the estimator. Note that this estimator has a rather low efficiency under normality and should be therefore only applied for time series with a high amount of outliers.}

\item{\code{"masarott0"}}{One uses the M estimator proposed by Masarotto (1987). Starting from initial estimations for the variance by the \code{\link{Qn}} and for the correlation by the GK estimator (using the Qn), new variances and correlations are iteratively computed until convergence using weights which are defined by the following psi function:
\deqn{\psi(x) = 3/(1+x).}
Note that this is not really a rank or sign based estimator. See Masarotto (1987) for more details.}
}
}

\references{
Croux, C. and Dehon, C. (2010): Influence functions of the Spearman and Kendall correlation measures, \emph{Statistical Methods & Applications}, vol. 19, 497--515.

Fried, R., Liboschik, T., Elsaid, H., Kitromilidou, S. and Fokianos, K. (2014): On outliers and interventions in count time series following GLMs, \emph{Austrian Journal of Statistics}, vol. 43, 181-193.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

Masarotto, D. (2003): Robust identification of autoregressive moving average models, \emph{Applied statistics}, vol. 36, 214--220.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.bireg}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "rank")
}
