\encoding{utf8}
\name{robspectrum}
\alias{robspectrum}

\title{
Robust Spectrum
}
\description{
Calculates the spectrum (or the spectral density) of a time series robustly via pre-whitening or via the autocorrelation function, a wrapper calling an internal function.
}
\usage{
robspectrum(x, method = c("pgram", "acf"), psifunc = smoothpsi,
	truncation = log(length(x), 10) * 10, bandwidth = FALSE,
	acf.fun = c("acfGK", "acfmedian", "acfmulti", "acfpartrank",
	"acfRA", "acfrank", "acftrim", "nonrob"), 
    arrob.method = c("yule-walker", "durbin-levinson", 
    "robustregression", "filter", "gm", "nonrob"), 
    var1 = FALSE, sdestim = c("Qn", "scaleTau2", "mad", "sd"),
    plot = TRUE, kernel = c("bartlett","parzen", "daniell",
    "rectangular"), ...)
}

\details{
There are two general methods to compute a robust spectrum.

If \code{method == "pgram"} the pre-whitening method as described in chapter 8.14 of Maronna et al. (2006) is applied. The method contains two spectrum estimations. First one fits an AR model robustly. See the help of \code{\link{arrob}} for an overview of possible methods which can be choosen by the argument \code{arrob.method}. Afterwards residuals are calculated and transformed by a given psi function \code{psifunc}, which should project exceptional large values to 0. The usual spectrum of this residuals is then combined with the theoretic spectrum of the AR fit, see Maronna et al. (2006) for details.

If \code{method == "ar"}, the sprectrum is calculated directly from acf using the well-known relationship between acf and spectrum. Here are also all methods of \code{\link{acfrob}} possible and can be choosen by \code{acf.fun}, see the help page of \code{\link{acfrob}} for more details and additional arguments. Because of computational reasons, it is usually not possible to estimate the whole acf up to lag \eqn{n-1}. Truncating the acf will always yield a smoothed estimate, which is in most applications also desired.

Since the autocorrelation does not contain informations about the scale of the timeseries, one can only derive the spectral density, which integrates to 1. If one is interested in the classical spectrum, one has to set \code{var1=FALSE}. In this case the estimated spectral density is transformed to a spectrum using a scale estimator which can be choosen using the argument \code{sdestim}.

The degree of smoothing can be either determined by setting a bandwith \code{bandwith} or a maximal lag of the acf \code{truncation}.
}

\arguments{
\item{x}{numeric vector containing the time series.}
\item{method}{character string naming the calculation method used: \code{"pgram"} for the pre-whitening method, \code{"acf"} for calculation from the acf, see details}
\item{psifunc}{psi function for filtering of residuals if \code{method == "pgram"}, see details.}
\item{truncation}{numeric indicating the truncation of the acf for spectrum estimation via fourier transform of acf. Otherwise this number is converted into the equivalent bandwidth.}
\item{bandwidth}{logical, if FALSE the bandwidth is determined by the truncation. Otherwise the stated number is used and converted into the truncation if necessary. A bandwidth of 0 prevent smoothing or at least tries to minimize smoothing if the method is "acf".}
\item{acf.fun}{character string naming the function to calculate the autocorrelation function.  See \code{"\link{acfrob}"} for details.  Do not use \code{"acfrobfil"}.}
\item{arrob.method}{character string naming the calculation method passed to \code{\link{arrob}} if \code{method == "pgram"}, see there for details.}
\item{var1}{logical whether the variance estimator is fixed to 1, resulting in a spectral density estimation if \code{method == "acf"} instead of a real spectrum.}
\item{sdestim}{character string naming the standard deviation estimator if \code{method == "acf"} and \code{var1 == FALSE}. Must be one of 'Qn', 'scaleTau2', 'mad' or 'sd'.}
\item{plot}{logical whether the result is plotted and returned invisible.}
\item{kernel}{character string naming the used smoothing kernel. Possible are Bartlett, Rectangular, Parzen and Daniell for method "acf" and "pgram".}
\item{...}{further arguments passed to the acf calculation by \code{\link{acfrob}} if \code{method == "acf"} or the AR calculation if \code{method == "pgram"}.}
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
robspectrum(tss, method = "acf", plot = FALSE)
}
