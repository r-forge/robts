\encoding{utf8}
\name{filterrob}
\alias{filterrob}

\title{
Robustly Filtered Time Series
}
\description{
Filters a given time series robustly. The filtered series can be applied to improve estimations. Optionally, the coefficients of a fitted AR model can be used.  See \enc{Dürre}{Duerre} et al. (2014) and Maronna et al. (2006), chapter 8.6 for details.
}
\usage{
filterrob(x, filter = NULL, method = c("direct", "ar"), p = 0, na.action = na.fail, psi.l=2,psi.0=3)
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{filter}{numeric vector of the coefficients of the AR model fitted to \code{x} if \code{method == "ar"}.}
	\item{method}{character string whether the filtering is performed directly or using model coefficients calculated previously.}
	\item{p}{order of AR model if \code{method == "direct"}.}
	\item{na.action}{function to be called to handle missing values.}
	\item{psi.l}{numeric determining the used psi function, see details.}
	\item{psi.0}{numeric determining the used psi function, see details.}
}

\details{
This function filters a timeseries relative to outliers. It is assumed that the timeseries follows an AR-process with gaussian innovations. If one knows the model one can determine the AR coefficients by the argument \code{filter} and set \code{method="ar"}. If one does not know the coefficients but at least the model order \code{p}, one has to set \code{method="direct"}. In this case the parameter are estimated by the robust filter algorithm, see the help of \code{\link{ARfilter}} for details. If one does not even know the order one should use the function \code{arrob} which determines the model order by an AIC criteria.
Based on the model, residuals are calculated and smoothed by a psi-function, which is the identity for absolute values smaller than \code{psi.l} and 0 for absolute values larger than \code{psi.0}. See \enc{Dürre}{Duerre} et al. (2014) for a formal definition of the psi function. 

Note that this function can be also abused to handle missing values. If one replaces missing values with abnormal high ones, the algorithm will replace them by the best linear prediction. 
}

\value{
	Named list of 2:
	\item{ts}{time series: the robustly filtered version of \code{x}.}
	\item{residuals}{residuals of the fitted AR model if \code{method == "ar"}.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): \emph{Robust Statistics: Theory and Methods}, Wiley, chapter 8.
}

\seealso{
The function \code{\link{ARfilter}} estimates an AR modell based on robustly filtered values.
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
