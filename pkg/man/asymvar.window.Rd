\encoding{utf8}
\name{asymvar.window}
\alias{asymvar.window}

\title{
Calculation of the long run variance based on subsampling
}
\description{
Computes the long run variance, which is required for change point testing, based on supsampling. See Dehling et al. (2013) for details.
}

\usage{
asymvar.window(x=x,overlapping=FALSE,obs=c("untransformed","ranks"),dd=c("independent","carlstein-cor","carlstein-Qn"),momentp=1,...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{overlapping}{logical indicating whether blocksums should be distinct or overlapping, see details.}
	\item{obs}{character string indicating whether the lon run variance of a cusum statistic or a Wilcoxon statistic is calculated.}
	\item{dd}{character string indicating how the blocklength should be determined. Must be one of '"independent"', '"carlstein-cor"' or '"carlstein-Qn"', see details.}
	\item{momentp}{numeric indicating what absolute centered mean is used to estimate the long run variance.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
List containing the folowing values:
\item{lrv}{estimated long run variance}
\item{blocklength}{used blocklength of the subsamples}
}

\details{
Cusum-type changepoint tests require an estimation of the long run variance. One possibility is to use subsampling methods, where one splits the sample into blocks and estimates the variance by a centered absolute empirical moment of these blocksums. There are different tuning options: setting the block length \eqn{l}, choosing the centered absolute moment \eqn{p} and characterising the block building. For the later, there are two options. 

If 'overlapping' is 'TRUE' one uses a \eqn{n-l} successive overlapping blocks. This method was proposed by Peligrad and Shao (1995). Otherwise one uses \eqn{\lfloor n/l \rfloor} distict blocks which was proposed by Carlstein (1986).

There are three different methods to determine a working blocklength  \eqn{l}:

if 'dd' is '"carlstein-cor"' one uses Carlstein-rule for an optimal blocklength. Herefore it is assumed that the timeseries follows an AR(1) modell. See Carlstein (1986) for more Information.

if 'dd' is '"carlstein-Qn"' one basically applies Carlstein-rule with a robust estimation of the AR(1) coefficient. More in detail it is estimated by the autocorrelation of the first lag by \code{\link{acfGK}} .

if 'dd' is '"independent"' one uses a thumb rule only based on the timeseries length \eqn{n}:
\deqn{round((3N)^(1/3)+1).}

}

\references{

Carlstein, E. (1986): The use of subseries values for estimating the variance of a general statistic from a stationary sequence, \emph{The Annals of Statistics}, vol. 14, 1171--1179.

Dehling, H., Fried, R., Shapirov, O., Wornowitzki, M. (2013): Estimation of the variance of partial sums of dependent processes, \emph{Statistics and Probability Letters}, vol. 83, 141--147.

Peligrad, M., Shao, Q. (1995): Estimation of the Variance of Partial Sums for rho-Mixing Random Variables, \emph{Journal of Multivariate Analyis}, vol. 152, 140--157. 
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}


\seealso{

The long run variance can be also estimated by \code{\link{asymvar.acf}} and \code{\link{asymvar.acfextra}}.

the robust estimation of the AR(1) coefficient is performed by \code{\link{acfGK}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.window(tss)
}
