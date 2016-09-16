\encoding{utf8}
\name{acfposmaker}
\alias{acfposmaker}

\title{
Transforms acf Estimations into Positive Semidefinite ones 
}
\description{
Transforms acf estimations into positive semidefinite ones by iterative projections. See \enc{Dürre}{Duerre} et al. (2015) for details.
}

\usage{
acfposmaker(acfv, ...)
}

\arguments{
	\item{acfv}{original acf estimation.}
	\item{...}{tuning parameters for internal functions, see Details.}
}

\value{
Numeric vector of transformed autocorrelations.
}

\details{

The function transforms a vector into a semidefinite one by the projection algorithm described in Al-Homidan (2006). First, based on the estimated acf, the corresponding correlation matrix (of the lagged time series) is build. Then the best positive semidefinite matrix with Toeplitz structure is computed by iteratively projecting onto positive semidefinite and toeplitz matrices. The algorithm stops either a maximal number of iterations is reached \code{maxit = 100} or the changes of the Frobenius norm is smaller then a given threshold \code{tol
 = 10^(-8)}, both values can be changed by the \citet{...} argument. See \enc{Dürre}{Duerre} (2015) for more details of the algorithm.
}  

\references{

Al-Homidan, S.: (2006): Semidefinite and second-order cone optimization approach for the toeplitz matrix approximation problem, \emph{Journal of Numerical Mathematics}, vol. 14, 1--15.

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfRA(tss,7)
}
