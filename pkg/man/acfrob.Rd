\encoding{utf8}
\name{acfrob}
\alias{acfrob}

\title{
Robust (Partial) Autocorrelation Function Estimation
}
\description{
Robustly estimates the (partial) autocorrelation function of a time series, a wrapper function calling an internal (partial) autocorrelation function. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
acfrob(x, lag.max = NULL,
    type = c("correlation", "partial"),
    fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank",
        "acfRA", "acfrank", "acfrobfil", "acftrim"),
    plot = TRUE, na.action = na.fail, p = NULL,
    robfiltype = c("emp", "pacf", "pacfMott"),
    partialtype = c("ranks", "durbin-levinson", "filter"), ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{numeric value of maximum lag at which to calculate the acf.  If \code{NULL} (default), \code{10 * log(length(x), 10)}, is chosen.}
	\item{type}{character string whether the partial autocorrelation is calculated.}
	\item{fun}{character string naming the function used for computation, see details.}
	\item{plot}{logical.  If \code{TRUE} (the default) the acf is plotted.}
	\item{na.action}{function to be called to handle missing values.}
	\item{p}{numeric value of order of AR fit used if \code{fun == acfrobfil}.}
	\item{robfiltype}{character string naming the type of robust filtering if \code{fun == "acfrobfil"}, see details.}
	\item{partialtype}{character string naming the type of pacf calculation if \code{type == "partial"}, see details.}
	\item{...}{further arguments passed to the respective internal autocorrelation function.}
}

\value{
Object of class \code{"acf"}, see \code{\link{acf}} for details.
}

\details{
	Possible autocorrelation functions and respective further arguments:
	
	\code{acfGK}: Gnanadesikan-Kettenring (GK) method.  Character string argument \code{method} giving the type of variance estimator to be used.  Possible values are default \code{"Qn"}, \code{"MAD"} or \code{"Tau"}.  Additional arguments passed to the function \code{\link[robustbase]{scaleTau2}} in package \code{robustbase}.
	
	\code{acfmedian}: based on mediancorrelation.  No further arguments.
	
	\code{acfmulti}: from multivariate correlation.  Character string argument \code{method} giving the correlation estimator to be used: default \code{"weightedMCD"}, \code{"rawMCD"}, \code{"Stahel-Donoho"}, \code{"S"}, \code{"M"}, \code{"reweight"}, or \code{"Tyler"}.  Additional arguments passed to the correlation estimator.
	
	\code{acfpartrank}: from partial autokorrelations and ranks.  Character string argument \code{method} giving the correlation estimator to be used: default \code{"spearman"}, \code{"gaussian"}, \code{"kendall"}, \code{"quadrant"} or \code{"masarotto"}.
	
	\code{acfRA}: RA approach.  Character string \code{Psi} giving Psi function to be used: default \code{"Huber"} or \code{"Tukey"}. Mean estimator function \code{meanvalue} (default \code{\link{median}}) and scale estimator function \code{scattervalue} (default \code{\link{mad}}).  Additional robustness parameters for \code{Psi} functions.
	
	\code{acfrank}: rank estimators.  Character string argument \code{method} giving the correlation estimator to be used: \code{"spearman"}, default \code{"gaussian"}, \code{"kendall"}, \code{"quadrant"} or \code{"masarotto"}.
	
	\code{acfrobfil}: using robustly filtered values.  Obligatory argument: Numeric value \code{p} of order of AR fit to be used.  Psi function \code{psifunc} to be used, defaults to \code{smoothpsi}.  Returned acf depends on \code{acfrob}'s argument \code{robfiltype}: default \code{"emp"} for the empirical acf of the robustly filtered timeseries, \code{"pacf"} for the estimated partial autocorrelations, \code{"pacfMott"} for the estimated partial autocorrelations after Mottonen et al. (1999).
	
	\code{acftrim}: trimmimg.  Numeric value \code{alpha} of part of largest and smallest observations to be trimmed, defaults to \code{0.1}.
	
	Types of pacf calculation, see \enc{Dürre}{Duerre} et al. (2014) for details:
	
	\code{ranks}: using rank estimators.
	
	\code{durbin-levinson}: The acf is calculated via one of the functions above.  Afterwards, the Durbin-Levinson algorithm is used to obtain the pacf.
	
	\code{filter}: using robustly filtered values.
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Mottonen, J., Koivunen, V. and Oja, H. (1999): Robust autocovariance estimation based on sign and rank correlation coefficients, In: \emph{Proceedings of the IEEE Signal Processing Workshop on Higher-Order Statistics}, 187--190.
}
\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
Classical, nonrobust estimation is provided by the function \code{\link{acf}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, plot = FALSE)
acfrob(tss, type = "partial", plot = FALSE)
}
