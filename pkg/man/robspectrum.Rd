\encoding{utf8}
\name{robspectrum}
\alias{robspectrum}

\title{
Robust Spectrum
}
\description{
Calculates the spectrum (or the spectral density) of a time series robustly via a periodogram or via the autocorrelation function, a wrapper calling an internal function.
}
\usage{
robspectrum(x, method = c("pgram", "ar"), psifunc = smoothpsi,
	truncation = log(length(x), 10) * 10, bandwidth = FALSE,
	acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank",
	"acfRA", "acfrank", "acftrim", "nonrob"), 
    arrob.method = c("yule-walker", "durbin-levinson", 
    "robustregression", "filter", "gm", "nonrob"), 
    var1 = FALSE, sdestim = c("Qn", "scaleTau2", "mad", "sd"),
    plot = TRUE, kernel = c("parzen", "daniell", "bartlett",
    "rectangular"), ...)
}

\details{
If \code{method == "pgram"}, the AR(p) model coefficients and the respective optimal p ist computed by \code{\link{arrob}}.  The smoothed periodogram is calculated from the filtered residuals of this model.  See Maronna et al. (2006), chapter 8.14.

If \code{method == "ar"}, the sprectrum is calculated directly from the acf.
}

\arguments{
\item{x}{numeric vector containing the time series.}
\item{method}{character string naming the calculation method used: \code{"pgram"} for the periodogram-based function, \code{"ar"} for calculation from the acf.}
\item{psifunc}{Psi function for filtering of residuals if \code{method == "pgram"}, currently no alternative.}
\item{truncation}{Truncation of the acf for spectrum estimation via fourier transform of acf. Otherwise this number is converted into the equivalent bandwidth.}
\item{bandwidth}{If FALSE the bandwidth is determined by the truncation. Otherwise the stated number is used and converted into the truncation if necessary. A bandwidth of 0 prevent smoothing or at least tries to minimize smoothing if the method is "ar".}
\item{acf.fun}{character string naming the function to calculate the autocorrelation function.  See \code{"\link{acfrob}"} for details.  Do not use \code{"acfrobfil"}.}
\item{arrob.method}{character string naming the calculation method passed to \code{\link{arrob}} if \code{method == "pgram"}, see there for details.}
\item{var1}{logical whether the variance estimator is fixed to 1, resulting in a spectral density estimation if \code{method == "ar"}.}
\item{sdestim}{character string naming the standard deviation estimator if \code{method == "ar"} and \code{var1 == FALSE}.}
\item{plot}{logical whether the result is plotted and retured invisible.}
\item{kernel}{The used smoothing kernel. Possible are Bartless, Rectangular, Parzen and Daniell for method "ar" and "pgram".}
\item{...}{further arguments passed to the acf calculation by \code{\link{acfrob}} if \code{method == "ar"}.}
}

\value{
An object of class \code{"spec"}, see \code{\link{spectrum}} for details.
}
\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): \emph{Robust Statistics: Theory and Methods}, Wiley, chapter 8.
}
\author{
Alexander \enc{Dürre}{Duerre} and Jonathan Rathjens
}

\seealso{
Classical, nonrobust estimation is provided by the function \code{\link{spectrum}}.
}
\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 50)
robspectrum(tss)
robspectrum(tss, method = "ar", plot = FALSE)
}
