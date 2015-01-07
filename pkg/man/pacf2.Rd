\encoding{utf8}
\name{pacf2}
\alias{pacf2}

\title{Partial Autocorrelations}

\description{
Calculates the partial autocorrelations of a time series using rank estimators, see \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
pacf2(x, lag.max = NULL, method = c("spearman", "gaussian", "kendall", "quadrant", "masarotto"),
     plot = TRUE, na.action = na.fail)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{numeric value of maximum lag at which to calculate the acf. If \code{NULL} (default), \code{10 * log(length(x), 10)}, is chosen.}
	\item{method}{character string giving the correlation estimator to be used: default \code{"spearman"}, \code{"gaussian"}, \code{"kendall"}, \code{"quadrant"} or \code{"masarotto"}, see \enc{Dürre}{Duerre} et al. (2014) for details.}
	\item{plot}{logical. If \code{TRUE} (the default) the acf is plotted.}
	\item{na.action}{function to be called to handle missing values.}
}

\value{
Object of class \code{"acf"}, see \code{\link{acf}} for details.
}
\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.
}
\author{
Alexander \enc{Dürre}{Duerre}
}


\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
pacf2(tss, lag.max = 20, plot = FALSE)
}
