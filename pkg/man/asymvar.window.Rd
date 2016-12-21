\encoding{utf8}
\name{asymvar.window}
\alias{asymvar.window}

\title{
Calculation of the Long Run Variance Based on Subsampling
}
\description{
Computes the long run variance, which is required for change point testing, based on subsampling. See Dehling et al. (2013) for details.
}

\usage{
asymvar.window(x, overlapping = FALSE, obs = c("untransformed", "ranks"),
               dd = c("independent", "carlstein-cor", "carlstein-Qn"),
               momentp = 1)
}

\arguments{
	\item{x}{numeric vector or univariate time series object.}
	\item{overlapping}{logical indicating whether block sums should be distinct or overlapping, see Details.}
	\item{obs}{character string indicating whether the long run variance of a cusum statistic (\code{"untransformed"}) or of a Wilcoxon statistic (\code{"ranks"}) is calculated.}
	\item{dd}{character string indicating how the block length should be determined, see Details.}
	\item{momentp}{numeric value indicating what absolute centered mean is used to estimate the long run variance.}
}

\details{
Cusum-type change point tests require an estimation of the long run variance. One possibility is to use subsampling methods, where one splits the sample into blocks and estimates the variance by a centered absolute empirical moment of these block sums. There are different tuning options: setting the block length \eqn{l}, choosing the centered absolute moment \eqn{p} and characterising the block building. For the later, there are two options. 

If \code{overlapping} is \code{TRUE}, one uses a \eqn{n-l} successive overlapping blocks. This method was proposed by Peligrad and Shao (1995). Otherwise, one uses \eqn{\lfloor n/l \rfloor} distinct blocks which was proposed by Carlstein (1986).

There are three different methods to determine a working block length \eqn{l}, which are chosen by the argument \code{dd}:
\describe{
\item{\code{"carlstein-cor"}}{Use Carlstein's rule for an optimal block length. Therefore it is assumed that the time series follows an AR(1) model. See Carlstein (1986) for more Information.}

\item{\code{"carlstein-Qn"}}{Basically apply Carlstein's rule with a robust estimation of the AR(1) coefficient. To be more specific, it is estimated by the autocorrelation of the first lag using \code{\link{acfrob.GK}}.}

\item{\code{"carlstein-Qn"}}{independent}{Use a thumb rule only based on the time series length \eqn{n}: \deqn{round((3N)^{1/3}+1).}}
}
}

\value{
List containing the following named elements:
\item{lrv}{estimated long run variance}
\item{blocklength}{used block length of the subsamples}
}

\references{
Carlstein, E. (1986): The use of subseries values for estimating the variance of a general statistic from a stationary sequence, \emph{The Annals of Statistics}, vol. 14, 1171--1179, \doi{10.1214/aos/1176350057}.

Dehling, H., Fried, R., Shapirov, O., Wornowizki, M. (2013): Estimation of the variance of partial sums of dependent processes, \emph{Statistics and Probability Letters}, vol. 83, 141--147, \doi{10.1016/j.spl.2012.08.012}.

Peligrad, M., Shao, Q. (1995): Estimation of the Variance of Partial Sums for rho-Mixing Random Variables, \emph{Journal of Multivariate Analysis}, vol. 152, 140--157, \doi{10.1006/jmva.1995.1008}. 
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}


\seealso{
The long run variance can be also estimated by \code{\link{asymvar.acf}} and \code{\link{asymvar.acfextra}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
asymvar.window(tss)
}
