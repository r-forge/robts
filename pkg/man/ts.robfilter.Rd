\encoding{utf8}
\name{ts.robfilter}
\alias{ts.robfilter}

\title{
Robust filtering of time series
}
\description{
Filters a time series robustly according to an AR model. The filtered series can be applied to improve estimations. See \enc{Dürre}{Duerre} et al. (2014) and Maronna et al. (2006), chapter 8.6 for details.
}
\usage{
ts.robfilter(ts, p, psifunc = super)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{p}{Order of the AR model used.}
	\item{psifunc}{Psi function, currently no alternative.}
}

\value{
Vector of the same length as \code{ts} containing the robustly filtered time series.
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Alexander \enc{Dürre}{Duerre}
}

\seealso{
	\code{\link{acfyohai}} for application.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
ts.robfilter(ts = tss, p = 6)
}
