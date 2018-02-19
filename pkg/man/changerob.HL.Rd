\encoding{utf8}
\name{changerob.HL}
\alias{changerob.HL}

\title{
Hodges-Lehmann Cumulative Sum Statistic
}
\description{
Computes the test trajectory for the Hodges-Lehmann change point test. See Dehling et al. (2015) for details.
}

\usage{
changerob.HL(x, var.method = c("window", "acf", "acfextra"),
             overlapping = TRUE, shiftcorrect = TRUE, borderN = 10, ...)
}

\arguments{
	\item{x}{numeric vector or univariate time series object.}
	\item{var.method}{character string defining the estimator of the long run variance, see Details.}
	\item{overlapping}{logical value indicating whether block sums should be distinct or overlapping. Only relevant if \code{var.method = "window"}, see Details.}
	\item{shiftcorrect}{logical value. If \code{TRUE}, the estimation of the long run variance considers a change in location by first estimating the change point and then subtracting the estimated jump height from the second part of the time series. This should considerably improve the power of the test.}
	\item{borderN}{integer value. If \code{shiftcorrect = TRUE}, the change point for the estimation of the log run variance is not searched below the first and last \code{borderN} observations.}
	\item{...}{further arguments passed to the respective internal function for calculation of the asymptotical long run variance.}
}

\details{
The Hodges-Lehmann change point test is a robust and nonparametric test against a change in location under short range dependence. It is preferable to the usual cusum test if the time series contains outliers or has heavy tails. See Dehling et al. (2015) for details. This function computes the test trajectory, which is basically a series of two sample tests, splitting the time series at every possible change point. The actual test is executed by the function \code{\link{changerob}}. 

The long run variance of the test statistic consists of two parts, the long run variance of the ranks of the time series \eqn{\sigma}{sigma} and the density at 0 of the difference of two independent random variables which have the same marginal distribution as an observation of the time series \eqn{u(0)}. There are several possibilities to estimate both. For the long run variance of the ranks \eqn{\sigma}{sigma} there are three options. More details on the different options for calculating the long run variance are given on the help page of the function \code{\link{changerob}}.

For estimating \eqn{u(0)}, a classical kernel estimator is used, in fact the function \code{\link{density}} of the stats package. There are several possibilities of differences to use. See the help page of \code{\link{densdiff}} for more details. 

Note that the procedure is computational demanding and can take some time for long time series (with more than 2000 observations).
}

\value{
Numeric vector which contains the whole path of the test statistic.
}

\references{

Carlstein, E. (1986): The use of subseries values for estimating the variance of a general statistic from a stationary sequence, \emph{The Annals of Statistics}, vol. 14, 1171--1179, \doi{10.1214/aos/1176350057}.

Dehling, H., Fried, R., Wendler, M. (2015): A robust method for shift detection in time series, preprint. 
\href{http://arxiv.org/abs/1506.03345}{arXiv 1506.03345}

Peligrad, M., Shao, Q. (1995): Estimation of the Variance of Partial Sums for rho-Mixing Random Variables, \emph{Journal of Multivariate Analysis}, vol. 152, 140--157, \doi{10.1006/jmva.1995.1008}. 
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}


\seealso{
The wrapper function \code{\link{changerob}}.

The long run variance is calculated by \code{\link{asymvar.window}}, \code{\link{asymvar.acf}} or \code{\link{asymvar.acfextra}}.

The density at 0 for the Hodges-Lehmann estimator is calculated by \code{\link{densdiff}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 0.5), n = 100)
trajectory <- changerob.HL(tss)
plot(trajectory)
}
