\encoding{utf8}
\name{asymvar.acf}
\alias{asymvar.acf}

\title{
Calculation of the long run variance based on a kernel estimator
}
\description{
Computes the long run variance, which is required for change point testing, by a kernel estimation.
}

\usage{
asymvar.acf(x,obs=c("untransformed","ranks"),cc=1.4,K=3,type=c("bartlett","trapezoid"),...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{obs}{character string indicating whether the lon run variance of a cusum statistic or a Wilcoxon statistic is calculated.}
	\item{cc}{numeric which determines the choosen bandwith, see details.}
	\item{K}{numeric which determines the choosen bandwith, see details.}
	\item{type}{character string determining the used kernel. Must be either '"bartlett"' or '"trapezoid"', see details.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
List containing the folowing values:
\item{lrv}{estimated long run variance}
\item{bandwith}{used bandwith}
}

\details{
Cusum-type changepoint tests require an estimation of the long run variance. One possibility is to use a kernel estimation wich downweights empirical autocorrelations of larger lags which one cannot reliable estimate. There are different tuning options regarding the bandwith \eqn{bn} and kernel. For the later only two are implemented. 

If 'kernel' is '"bartlett"' one uses the classical Bartlett kernel which is defined as
\deqn{k(x)= 1 - x/ln~~\mbox{for} |x|\leq ln~~\mbox{and} 0~~\mbox{for} |x|>ln}

If 'kernel' is '"trapezoid"' the following flat-top kernel is used
\deqn{k(x)= 1 ~~\mbox{for} |x|\leq ln~~\mbox{and}~~ 1-|x-ln|/ln~~\mbox{for} 2ln\leq |x|>ln~~and 0~~\mbox{for} x>2ln}

For the bandwith a data-adaptive thumbrule depending on two tuning-parameter is used. The basic idea is to determine after which lag the autocorrelation is negligible. There for one looks if the emprirical acf is smaller than
\deqn{cc \sqrt{\log(N)_{10}/N}}
if the next \eqn{K} autocorrelation also fuilfill this criteria this lag is considered as bandwith.
}


\seealso{

The long run variance can be also estimated by \code{\link{asymvar.window}} and \code{\link{asymvar.acfextra}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.acf(tss)
}
