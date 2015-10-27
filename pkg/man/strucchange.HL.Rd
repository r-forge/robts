\encoding{utf8}
\name{strucchange.HL}
\alias{strucchange.HL}

\title{
Hodges-Lehmann cumulative sum
}
\description{
Computes the test trajectory for the Hodges-Lehmann change point test. See Dehling et al. (2015) for details.
}

\usage{
strucchange.HL(x,varmethod=c("window","acf","acfextra"),overlapping=TRUE,shiftcorrect=TRUE,...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{varmethod}{character string defining the estimator for the long run variance. Must be one of '"window"', '"acf"' or '"acfextra"', see details.}
	\item{overlapping}{logical indicating whether blocksums should be distinct or overlapping. Only relevant if \code{varmethod == window}, see details.}
	\item{shiftcorrect}{logical, if TRUE the estimation of the long run variance considers a change in location by first estimating the changepoint and then subtracting the estimated jump height from the second part of the time series}
	\item{borderN}{numeric, if \code{shiftcorrect == TRUE} the changepoint for the estimation of the log run variance is not searched below the first and last borderN observations.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
Numeric of length timeseries-length - 1 which contains the whole path of the test statistic.
}

\details{
The Hodges-Lehmann change point test is a robust and nonparametric test against a change in location under short range dependence. It is preferable to the usual Cusum test if the timeseries contains outliers or has heavy tails. See Dehling et al. (2015) for details. This funtion computes the testtrajectory which is basically a series of two sample tests, splitting the timeseries at every possible break point. The actual test is excuted by the funtion \code{\link{rob.change}}. 

The long run variance of the test statistic consists of two parts, the long run variance of the ranks of the timeseries \eqn{\sigma} and the density at 0 of the difference of two independent random variables which have the same marginal distribution as an observation of the timeseries \eqn{u(0)}. There are several possibilities to estimate both. For the long run variance of the ranks there are three options. 

If 'varmethod' is '"window"' one uses a subsampling estimator, which is based on centred absolute moments of block sums. One can either use distinct blocks, which was proposed by Carlstein (1986), and is applied if one sets \code{"overlapping=FALSE"} or overlapping blocks, which was proposed by Peligrad and Shao (1995), and works in case of \code{"overlapping=TRUE"}. For more details how to set the blocklength and other tuning possibilities, see the help page of \code{\link{asymvar.window}}.

if 'varmethod' is '"acf"' one uses a kernel estimator, which weights estimated autocovariances. For more details how to set the bandwith and other tuning possibilities, see the help page of \code{\link{asymvar.acf}}.

if 'varmethod' is '"acfextra"' one uses an extrapolation of the acf based on an ar-fit. For more details how to set the bandwith and other tuning possibilities, see the help page of \code{\link{asymvar.acfextra}}.

To estimate \eqn{u(0)} a classical kernel estimator is used, in fact the function \code{\link{density}} of the stats package. There are several possibilities of differences to use. See the help page of \code{\link{densdiff}} for more details. 
}

\references{

Carlstein, E. (1986): The use of subseries values for estimating the variance of a general statistic from a stationary sequence, \emph{The Annals of Statistics}, vol. 14, 1171--1179.

Dehling, H., Fried, R., Wendler, M. (2015): A robust method for shift detection in time series, preprint. 
\href{http://arxiv.org/abs/1506.03345}{arvix 1506.03345}

Peligrad, M., Shao, Q. (1995): Estimation of the Variance of Partial Sums for rho-Mixing Random Variables, \emph{Journal of Multivariate Analyis}, vol. 152, 140--157. 
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}


\seealso{
Further non robust change point tests are provided by the package \code{strucchange}

The wrapper function \code{\link{rob.change}}

The long run variance is calculated by \code{\link{asymvar.window}}, \code{\link{asymvar.acf}} or \code{\link{asymvar.acfextra}}

The density at 0 for the Hodges-Lehmann estimator is calculated by \code{\link{densdiff}}
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
trajectory <- strucchange.HL(tss)
plot(trajectory)
}
