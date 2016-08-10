\encoding{utf8}
\name{acfpartrank}
\alias{acfpartrank}

\title{
Robust Autocorrelation Estimation Based on Partial Autocorrelations 
}
\description{
Robustly estimates the autocorrelation function of a time series based on a robust estimation of the partial autocorrelation function. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfpartrank(timeseries,maxlag,partial.method="spearman")
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{maxlag}{numeric value of maximum lag at which to calculate the acf.}
	\item{method}{character string indicating the used correlation estimator. Must be one of 'masarotto', 'gaussian','spearman','kendall','quadrant', see details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

This function estimates the autocorrelation function based on partial autocorrelations. This procedure was first proposed by Masarotto (1987) and has the advantage that it produces positive definite and therefore valid acf estimations. To estimate the pacf, the estimation procedure described in \enc{Möttönen}{Moettoennen} et al. (1999) is used. It depends on successive correlation estimations of forward and backward residuals. There are many possibilities to estimate these correlations.

	If \code{partial.method="masarotto"} one follows the original proposal of Masarotto (1987), which is a kind of M estimator. Starting from initial estimations for the variance by the \code{\link{Qn}} and for the correlation by the GK (using the Qn) estimator, new variances and correlations are iteratively computed until convergence, using weights which are defined by the following psi function

\deqn{\psi(x)=3/(1+x).}
See Masarotto (1987) for more details.


	If \code{method="gaussian"} one uses the Gaussian rank correlation which is the usual correlation applied to transformed values. Denote \eqn{R(i)} the rank of \eqn{X(i)} and \eqn{\Phi} the distribution function of the normal distribution, then the transformation is defined as

\deqn{\Phi^{-1}(R(i)/(n+1)).}

This estimator has the same asymptotic efficiency as the usual empirical correlation under normality, but is robust to some degree.

	If \code{method="spearman"} one uses the well known Spearman's rho with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Spearmans rho.

	If \code{method="kendall"} one uses the well known Kendall's tau with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Kendalls tau. Note that Kendall's tau can get rather time consuming for large timeseries (e.g. greater than 10000).

	If \code{method="quadrant"} one uses the quadrant correlation with a consistency correction under normality. See Croux and Dehon (2010) for this correction and more information about the estimator. Note that this estimator has a rather low efficiency under normality and should be therefore only applied for timeseries with a high amount of outliers.

The computation of the acf from the pacf is carried out by the algorithm described in \enc{Möttönen}{Moettoennen} et al. (1999)
}  

\references{
Croux, C. and Dehon, C. (2010): Influence functions of the Spearman and Kendall correlation measures, \emph{Statistical Methods & Applications}, vol 19, 497--515.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

Masarotto, D. (2003): Robust identification of autoregressive moving average models, \emph{Applied statistics}, vol. 36, 214--220.

\enc{Möttönen}{Moettoennen}, J., Koivunen, V. and Oja, H. (1999): Robust autocovariance estimation based on sign and rank correlation coefficients, In: \emph{Proceedings of the IEEE Signal Processing Workshop on Higher-Order Statistics}, 187--190.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfpartrank(tss,7)
}
