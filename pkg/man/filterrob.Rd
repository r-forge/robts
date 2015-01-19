\encoding{utf8}
\name{filterrob}
\alias{filterrob}

\title{
Robustly Filtered Time Series
}
\description{
Filters a given time series robustly. The filtered series can be applied to improve estimations. Optionally, the coefficients of a fitted AR model can be used. See \enc{Dürre}{Duerre} et al. (2014) and Maronna et al. (2006), chapter 8.6 for details.
}
\usage{
filterrob(x, filter = NULL, method = c("direct", "ar"), p = 0, psifunc = smoothpsi, na.action = na.fail)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{filter}{numeric vector of the coefficients of the AR model fitted to \code{x} if \code{method == "ar"}.}
	\item{method}{character string whether the filtering is performed directly or using model coefficients calculated previously.}
	\item{p}{order of AR model if \code{method == "direct"}.}
	\item{psifunc}{Psi function, currently no alternative.}
	\item{na.action}{function to be called to handle missing values.}
}

\value{
	Named list of 2:
	\item{ts}{time series: the robustly filtered version of \code{x}.}
	\item{residuals}{residuals of the fitted AR model if \code{method == "ar"}.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
filterrob(tss, p = 6)
phis <- arrob(tss, aic = FALSE, order.max = 3)$ar
filterrob(tss, filter = phis, method = "ar")
}
