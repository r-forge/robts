\encoding{utf8}
\name{acfrob.bireg}
\alias{acfrob.bireg}

\title{
Robust Autocorrelation Estimation Based on ...
}

\description{
Robustly estimates the autocorrelation function of a time series based on... TO BE INCLUDED. See \enc{Dürre}{Duerre} et al. (2015) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "acfregression"}.
}

\usage{
acfrob.bireg(x, lag.max, regression.method = c("MM", "LTS", "L1"))
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

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, approach = "regression")
}
