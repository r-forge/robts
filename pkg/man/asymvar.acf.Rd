\encoding{utf8}
\name{asymvar.acf}
\alias{asymvar.acf}

\title{
Calculation of the Long Run Variance Based on a Kernel Estimator
}
\description{
Computes the long run variance, which is required for change point testing, by a kernel estimation.
}

\usage{
asymvar.acf(x, obs = c("untransformed", "ranks"), cc = 1.4, K = 3,
            type = c("bartlett", "trapezoid"))
}

\arguments{
	\item{x}{numeric vector or univariate time series object.}
	\item{obs}{character string indicating whether the long run variance of a cusum statistic (\code{"untransformed"}) or of a Wilcoxon statistic (\code{"ranks"}) is calculated.}
	\item{cc}{numeric value which determines the chosen bandwidth, see Details.}
	\item{K}{numeric value which determines the chosen bandwidth, see Details.}
	\item{type}{character string determining the used kernel, see Details.}
}

\value{
List containing the following named elements:
\item{lrv}{estimated long run variance}
\item{bandwith}{used bandwidth}
}

\details{
Cusum-type change point tests require an estimation of the long run variance. One possibility is to use a kernel estimation which downweights empirical autocorrelations of larger lags which one cannot reliable estimate. There are different tuning options regarding the bandwidth \eqn{bw} and kernel. For the later only two are implemented: 

If \code{kernel} is \code{"bartlett"}, one uses the classical Bartlett kernel which is defined as
\deqn{k(x)= 1 - x/bn~~\mbox{for}~|x|\leq bw~~\mbox{and}~~ 0~~\mbox{for}~|x|>bw}

If \code{kernel} is \code{"trapezoid"}, the following flat-top kernel is used
\deqn{k(x)= 1 ~~\mbox{for}~ |x|\leq bw,~~ 1-|x-bw|/bw~~\mbox{for}~2bw\leq |x|>bw~~and 0~~\mbox{for}~x>2bw}

For the bandwidth a data-adaptive thumb rule depending on two tuning parameters is used. The basic idea is to determine after which lag the autocorrelation is negligible. Therefor one looks if the empirical acf is smaller than
\deqn{cc \sqrt{\log_{10}(N)/N}.}
If the next \eqn{K} autocorrelation also fulfill this criteria, this lag is considered as bandwidth.
}

\seealso{
The long run variance can be also estimated by \code{\link{asymvar.window}} and \code{\link{asymvar.acfextra}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.acf(tss)
}
