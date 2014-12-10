\encoding{utf8}
\name{ARopt.filter}
\alias{ARopt.filter}

\title{
Optimal order of AR model calculated via robust filtering
}
\description{
Calculates the optimal order an AR(p) model via robust filtering. See \enc{Dürre}{Duerre} et al. (2014) and Maronna et al. (2006), chapter 8.6 for details.
}

\usage{
ARopt.filter(tss, pmax = NULL, psifunc = smoothpsi)
}

\arguments{
	\item{tss}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{pmax}{Maximum order of AR model to consider. If \code{NULL}, the default minimum of \eqn{(T - 1) / 4} and \eqn{10 * lg(T)} is used.}
	\item{psifunc}{Psi function, currently no alternative.}
}

\value{
The optimal p.
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R., and Liboschik, T (2014): Robust estimation of (partial) autocorrelation, working paper.

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}

\author{
Alexander \enc{Dürre}{Duerre} and Jonathan Rathjens
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
ARopt.filter(tss, pmax = 7)
}
