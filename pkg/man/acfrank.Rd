\encoding{utf8}
\name{acfrank}
\alias{acfrank}

\title{
Robust autocorrelation estimation based on signs and ranks 
}
\description{
Robustly estimates the autocorrelation function of a time series based nonparametric estimators using ranks and signs. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfrank(timeseries,maxlag,method="gaussian")
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
	\item{method}{character string indicating the used correlation estimator. Must be one of '"masarotto"', '"gaussian"','"spearman"','"kendall"','"quadrant"', see details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

The function estimates the acf using rank correlation estimators, which estimate the autocorrelations individually for every lag. The approach works also well for count time series, see Fried et al. (2014). There are several possible rank estimators available


	If \code{method="gaussian"} one uses the Gaussian rank correlation which is the usual correlation applied to transformed values. Denote \eqn{R(i)} the rank of \eqn{X(i)} and \eqn{\Phi} the distribution function of the normal distribution, then the transformation is defined as

\deqn{\Phi^{-1}(R(i)/(n+1)).}

This estimator has the same asymptotic efficiency as the usual empirical correlation under normality, but is robust to some degree.

	If \code{method="spearman"} one uses the well known Spearmans rho with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Spearmans rho.

	If \code{method="kendall"} one uses the well known Kendalls tau with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Kendalls tau. Note that Kendalls tau can get rather time consuming for large timeseries (e.g. greater than 10000).

	If \code{method="quadrant"} one uses the quadrant correlation with a consistency correction under normality. See Croux and Dehon (2010) for this correction and more information about the estimator. Note that this estimator has a rather low efficiency under normality and should be therfore only applied for timeseries with a high amount of outliers.

If \code{method="masarotto"} one uses and M estimator proposed by Masarotto (1987). Starting from initial estimations for the variance by the \code{\link{Qn}} and for the correlation by the GK (using the Qn) estimator, new variances and correlations are iteratively computed until convergence using weights which are defined by the following psi function

\deqn{\psi(x)=3/(1+x).}
Note that this is not really a rank or sign based estimator. See Masarotto (1987) for more details.
}

\references{

Croux, C. and Dehon, C. (2010): Influence functions of the Spearman and Kendall correlation measures, \emph{Statistical Methods & Applications}, vol 19, 497--515.

Fried, R., Liboschik, T., Elsaid, H., Kitromilidou, S. and Fokianos, K. (2014): On outliers and interventions in count time series following GLMs, \emph{Austrian Journal of Statistics}, vol 43, .181--193.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Masarotto, D. (2003): Robust identification of autoregressive moving average models, \emph{Applied statistics}, vol. 36, 214--220.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrank(tss,7)
}
