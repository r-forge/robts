\encoding{utf8}
\name{rob.change}
\alias{rob.change}

\title{
Robust Change Point Detection
}
\description{
Robust and non robust tests for changes in location and scale based on U-statistics and U-quantiles under short range dependence. See Dehling et al. (2013) and Dehling et al. (2015) for details.
}

\usage{
rob.change(x, property = c("location", "scale"), procedure = c("HL", "Wilcoxon", "Cusum"),
    conf.level = 0.95, alternative = c("two.sided", "less", "greater"),
    varmethod = c("window", "acf", "acfextra"), overlapping = FALSE, shiftcorrect = TRUE, 
    borderN = 10, plot = FALSE, ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{property}{character string whether one tests against a change in location or scale}
	\item{procedure}{character string defining the applied change point test. Possible are the non robust cusum test, the Wilcoxon test and the Hodges-Lehmann test. See Details for more information about the tests.}
	\item{conf.level}{numeric indicating the confidence level of the test.}
	\item{alternative}{character string which specifies the alternative hypothesis. Must be one of 'two.sided', 'less' or 'greater'.}
	\item{varmethod}{character string defining the estimator of the long run variance. Must be one of 'window", 'acf' or 'acfextra', see Details.}
	\item{overlapping}{logical indicating whether block sums should be distinct or overlapping. Only relevant if \code{varmethod == "window"}, see Details.}
	\item{shiftcorrect}{logical, if TRUE the estimation of the long run variance considers a change in location by first estimating the change point and then subtracting the estimated jump height from the second part of the time series. This should considerably improve the power of the test.}
	\item{borderN}{numeric, if \code{shiftcorrect == TRUE} the change point for the estimation of the log run variance is not searched below the first and last \code{borderN} observations.}
	\item{plot}{logical, if TRUE a plot with the trajectory of the test statistic including confidence bands is drawn.}	
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
Object of class \code{"htest"} containing the following values:
	\item{statistic}{the value of the test statistic. Under the null the test statistic follows, dependent of the chosen alternative, (asymptotically) a Kolmogorov Smirnov distribution.}
	\item{p.value}{the p-value of the test.}
	\item{estimate}{the estimated time of change.}
	\item{null.value}{the jump height of the at most one change point model, which is under the null hypothesis always 0.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a characters string indicating the chosen test.}
	\item{data.name}{a character giving the names of the data.}
	\item{trajectory}{the cumulative sum on which the tests are based on. Could be used for additional plots.}
}

\details{
The cusum type test procedures are originally designed to test for a constant location but are also able to detect changes in scale using the following transformation
\deqn{\tilde{x}_i = log(|x_i-median(x)|)}
of the original time series \eqn{x}. As usual cusum type tests are very powerful if there is one structural break but might suffer in case of multiple changes, especially if they are in different directions. There are three different tests implemented.

If 'procedure' is 'Cusum' the usual Cusum test for a change in location is applied. Because of its sensitivity relative to particular large or small values one should be careful in case of outliers. Furthermore, the test loses a lot of power if the time series is heavy tailed. In both cases one should apply one of the following tests.

If 'procedure' is 'Wilcoxon' one basically applies a two sample Mann-Whitney-Wilcoxon test for every possible change-point and choose a kind of maximal test statistic. The exact definition and details can be found in Dehling et al. (2013). By using ranks the tests gets less sensitive to outliers and does not suffer under heavy tailed distributions. Even under Gaussian processes the power of the test is not much worse than that of the cusum test. However, the test has little power against changes far from the center of the time series.

If 'procedure' is 'HL' one basically applies a two sample Hodges-Lehmann test for every possible change-point and choose a kind of maximal test statistic. The exact definition and details can be found in Dehling et al. (2015). As the Mann-Whitney-Wilcoxon test the Hodges-Lehmann test is less sensitive to outliers and does not suffer under heavy tailed distributions. The power of the test is in fact comparable to the Wilcoxon based one if the change-point is in the center of the time series but outperforms the former if the change is nearer to the margins. The test procedure is computationally demanding and can take some time for large time series (with more than 2000 observations).

There are several possibilities to estimate the long run variance. Three of them with further variations are implemented.

If 'varmethod' is 'window' one uses a subsampling estimator, which is based on centered absolute moments of block sums. One can either use distinct blocks, which is proposed by Carlstein (1986) and is applied if one sets \code{"overlapping = FALSE"} or overlapping blocks, which is proposed by Peligrad and Shao (1995) and is applied in case of \code{"overlapping = TRUE"}. For more details how to set the block length and other tuning possibilities, see the help page of \code{\link{asymvar.window}}.

if 'varmethod' is 'acf' one uses a kernel estimator, which weights estimated autocovariances. For more details how to set the bandwidth and other tuning possibilities, see the help page of \code{\link{asymvar.acf}}.

if 'varmethod' is 'acfextra' one uses an extrapolation of the acf based on an AR fit. For more details how to set the bandwidth and other tuning possibilities, see the help page of \code{\link{asymvar.acfextra}}.

For the Hodges-Lehmann estimator one needs an additional estimation of the density at 0. See the help page of \code{\link{densdiff}} how to adjust tuning constants for that.
}

\references{

Carlstein, E. (1986): The use of subseries values for estimating the variance of a general statistic from a stationary sequence, \emph{The Annals of Statistics}, vol. 14, 1171--1179.

Dehling, H., Fried, R., Wendler, M. (2015): A robust method for shift detection in time series, preprint. 
\href{http://arxiv.org/abs/1506.03345}{arXiv 1506.03345}

Dehling, H., Fried, R., Garcia, I., Wendler, M. (2013): Change-point detection under dependence based on two-sample U-statistics, preprint. 
\href{http://arxiv.org/abs/1304.2479}{arXiv 1304.2479}

Peligrad, M., Shao, Q. (1995): Estimation of the Variance of Partial Sums for rho-Mixing Random Variables, \emph{Journal of Multivariate Analysis}, vol. 152, 140--157. 
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}


\seealso{
Further non robust change-point tests are provided by the packages \code{changepoint} and \code{strucchange}.

rob.change calls the subroutines \code{\link{strucchange.cusum}}, \code{\link{strucchange.wilcox}} or \code{\link{strucchange.HL}}.

The long run variance is calculated by \code{\link{asymvar.window}}, \code{\link{asymvar.acf}} or \code{\link{asymvar.acfextra}}.

The density at 0 for the Hodges-Lehmann estimator is calculated by \code{\link{densdiff}}
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
rob.change(tss,plot=TRUE)
}
