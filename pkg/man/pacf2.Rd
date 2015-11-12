\encoding{utf8}
\name{pacf2}
\alias{pacf2}

\title{
Robust Estimation of the pacf Using Signs and Ranks
}
\description{
Robustly estimates the partial autocorrelation function of a time series based on signs anr ranks. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
pacf2(x, lag.max = NULL, method = c("spearman", "gaussian",
    "kendall", "quadrant", "masarotto"),na.action = na.fail)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{numeric value of maximum lag at which to calculate the acf.}
	\item{method}{character string indicating the used correlation estimator. Must be one of 'masarotto', 'gaussian','spearman','kendall','quadrant', see details.}
	\item{na.action}{function to be called to handle missing values.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{

This function estimates the partial autocorrelations. To estimate the pacf, the estimation prodedure described in \enc{Möttönen}{Moettoennen} et al. (1999) is used. It depends on successive correlation estimations of forward and backward residuals. There are many possibilities to estimate these correlations.

	If \code{method="masarotto"} one follows the original proposal of Masarotto (1987), which is a kind of M estimator. Starting from initial estimations for the variance by the \code{\link{Qn}} and for the correlation by the GK estimator using the Qn, new variances and correlations are iteratively computed until convergence, using weights which are defined by the following psi function

\deqn{\psi(x)=3/(1+x).}
See Masarotto (1987) for more details.


	If \code{method="gaussian"} one uses the Gaussian rank correlation which is the usual correlation applied to transformed values. Denote \eqn{R(i)} the rank of \eqn{X(i)} and \eqn{\Phi} the distribution function of the normal distribution, then the transformation is defined as

\deqn{\Phi^{-1}(R(i)/(n+1)).}

This estimator has the same asymptotic efficiency as the usual empirical correlation under normality, but is robust to some degree.

	If \code{method="spearman"} one uses the well known Spearmans rho with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Spearmans rho.

	If \code{method="kendall"} one uses the well known Kendalls tau with a consistency correction under normality. See Croux and Dehon (2010) for this correction and \code{\link{cor}} for more information about Kendalls tau. Note that Kendalls tau can get rather time consuming for large timeseries (e.g. greater than 10000).

	If \code{method="quadrant"} one uses the quadrant correlation with a consistency correction under normality. See Croux and Dehon (2010) for this correction and more information about the estimator. Note that this estimator has a rather low efficiency under normality and should be therefore only applied for timeseries with a high amount of outliers.

}  

\references{
Croux, C. and Dehon, C. (2010): Influence functions of the Spearman and Kendall correlation measures, \emph{Statistical Methods & Applications}, vol 19, 497--515.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Masarotto, D. (2003): Robust identification of autoregressive moving average models, \emph{Applied statistics}, vol. 36, 214--220.

\enc{Möttönen}{Moettoennen}, J., Koivunen, V. and Oja, H. (1999): Robust autocovariance estimation based on sign and rank correlation coefficients, In: \emph{Proceedings of the IEEE Signal Processing Workshop on Higher-Order Statistics}, 187--190.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
pacf2(tss,7)
}
