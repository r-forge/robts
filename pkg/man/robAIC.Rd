\encoding{utf8}
\name{robAIC}
\alias{robAIC}

\title{
Robust AIC of an AR(p) model
}
\description{
Calculates a robust version of Akaikes information criterion for the AR(p) model fit of a given time series. See Maronna et al. (2006), chapter 8.6.6 for details.
}
\usage{
robAIC(ts, p, psifunc = super)
}

\arguments{
	\item{ts}{Vector of observations sorted as time series, without \code{NA}s.}
	\item{p}{Order of the AR model.}
	\item{psifunc}{Psi function, currently no alternative.}
}

\value{
Value of the AIC.
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}
\author{
Alexander \enc{Dürre}{Duerre}
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
robAIC(ts = tss, p = 6)
}
