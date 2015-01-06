\encoding{utf8}
\name{acfrob}
\alias{acfrob}

\title{
Robust Autocorrelation Function Estimation
}
\description{
Robustly estimates the autocorrelation function of a time series, a wrapper function calling an internal autocorrelation function. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfrob(x, lag.max = NULL,
		fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA", "acfrank", "acfrobfil", "acftrim"),
		plot = TRUE, na.action = na.fail, p = NULL, robfiltype = c("emp", "pacf", "pacfMott"), ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{numeric value of maximum lag at which to calculate the acf.}
	\item{fun}{character string naming the function used for computation, see details.}
	\item{plot}{logical. If \code{TRUE} (the default) the acf is plotted.}
	\item{na.action}{function to be called to handle missing values.}
	\item{p}{numeric value of order of AR fit used if \code{fun == acfrobfil}.}
	\item{robfiltype}{character string naming the type of robust filtering if \code{fun == acfrobfil}, see details.}
	\item{...}{further arguments passed to the respective internal autocorrelation function.}
}

\value{
Object of class \code{"acf"}, see \code{\link{acf}} for details.
}

\details{
	Possible autocorrelation functions and respective further arguments:
	
	\code{acfGK}:
	
	\code{acfmedian}:
	
	\code{acfmulti}:
	
	\code{acfpartrank}:
	
	\code{acfRA}:
	
	\code{acfrank}:
	
	\code{acfrobfil}:
	
	\code{acftrim}:
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.
}
\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik, and Jonathan Rathjens
}


\seealso{
Classical, nonrobust estimation is provided by the function \code{\link{acf}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, plot = FALSE)
acfrob(tss, fun = "acfrobfil", p = 3)
}
