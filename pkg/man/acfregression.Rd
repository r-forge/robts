\encoding{utf8}
\name{acfregression}
\alias{acfregression}

\title{
Robust Autocorrelation Estimation Based on Residual Autocorrelation
}

\description{
Robustly estimates the autocorrelation function of a time series based on... TO BE INCLUDED. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfregression"}.
}

\usage{
acfregression(x, lag.max, regression.method = c("MM", "LTS", "L1"))
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{regression.method}{character string indicating the regression method, see Details.}
}

\value{
Numeric vector of estimated autocorrelations.
}

\details{
TO BE INCLUDED.
}  

\references{

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
acfrob(tss, approach = "acfregression")
}
