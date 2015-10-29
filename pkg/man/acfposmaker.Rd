\encoding{utf8}
\name{acfposmaker}
\alias{acfposmaker}

\title{
Transforms acf estimations into positive semidefinite ones 
}
\description{
Transforms acf estimations into positive semidefinite ones by iterative projections. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfposmaker(acfv,...)
}

\arguments{
	\item{acfv}{original acf estimation.}
	\item{...}{tuning parameters for internal functions, see details.}
}

\value{
Numeric vector of transformed autocorrelations.
}

\details{

The function transforms a vector into semidefinite one by a projection algorithm, proposed by Al-Homidan (2006). First based on the estimated acf the corresponding correlation matrix (of the lagged timeseries) is build. Then the best positive semidefinite matrix with toeplitz structure is computed by iteratively projecting onto positive semidefinite and toeplitz matrices. The algorithm stops either a maximal number of iterations is reached \code{maxit=100} or the changes of the Frobeniusnorm is smaller then a given threshold \code{tol
=10^(-8)}, both values can be changed by the ... argument. See \enc{Dürre}{Duerre} for more details of the algorithm.
}  

\references{

Al-Homidan, S.: (2006): Semidefinite and second-order cone optimization approach for the toeplitz matrix approximation problem, \emph{Journal of Numerical Mathematics}, vol. 14, 1--15.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfRA(tss,7)
}
