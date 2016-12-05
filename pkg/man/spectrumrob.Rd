\encoding{utf8}
\name{spectrumrob}
\alias{spectrumrob}
\alias{spectrumrob.pgram}
\alias{spectrumrob.acf}

\title{
Robust Spectrum
}                           
\description{
Calculates the spectrum (or the spectral density) of a time series robustly via pre-whitening or via the autocorrelation function, a wrapper calling an internal function.
}
\usage{
spectrumrob(x, method = c("pgram", "acf"), plot = TRUE, ...)

spectrumrob.pgram(x, psifn = function(x) M_psi(x, type="smooth"),
                  arrob.method = c("yw", "regression", "filter", "gm", "nonrobust"),
                  kernel = c("parzen", "bartlett", "rectangular", "daniell"),
                  truncation = round(log(length(x), 10) * 10), bandwidth = NULL,
                  na.action = na.fail, series = deparse(substitute(x)), ...)
                  
spectrumrob.acf(x, acf.approach = c("GK", "median", "multi", "partrank",
            "RA", "rank", "trim", "bireg", "nonrobust"),
            kernel = c("parzen", "bartlett", "rectangular", "daniell"),
            truncation = round(log(length(x), 10) * 10), bandwidth = NULL,
            var1 = FALSE, scalefn = Qn, na.action = na.fail,
            series = deparse(substitute(x)), ...)
}

\arguments{
\item{x}{numeric vector containing the time series.}
\item{method}{character string giving the calculation method used: \code{"pgram"} for the pre-whitening method and \code{"acf"} for calculation from a robustly estimated acf, see Details.}
\item{plot}{logical whether the result is plotted and returned invisibly (the default).}
\item{psifn}{psi function for filtering of residuals if \code{method = "pgram"}, see Details.}
\item{arrob.method}{character string giving the calculation method passed to \code{\link{arrob}} if \code{method = "pgram"}, see there for details. One can additionally choose a nonrobust approach by \code{"nonrobust"}.}
\item{acf.approach}{character string giving the approach for robustly estimating the autocorrelation function which is passed to \code{\link{acfrob}} if \code{method = "acf"}, see there for details. One can additionally choose a nonrobust approach by \code{"nonrobust"} but not the approach \code{"filter"}.}
\item{kernel}{character string naming the used smoothing kernel. Possible kernels are \code{"parzen"} (the default), \code{"bartlett"}, \code{"rectangular"} and \code{"daniell"}.}
\item{truncation}{integer value indicating the truncation of the acf for spectrum estimation via Fourier transform of acf. Otherwise this number is converted into the equivalent bandwidth.}
\item{bandwidth}{numeric value giving the bandwidth for smoothing as a proportion of the frequency range. The given bandwidth is converted into the corresponding value for the argument \code{truncation}, which overrides any value given by that argument. Only if this argument is \code{NULL}, the given truncation is converted into the corresponding bandwidth and used for smoothing. A bandwidth of zero prevents from smoothing or at least tries to minimize smoothing if the \code{method = "acf"}.}
\item{var1}{logical whether the variance estimator is fixed to 1, resulting in a spectral density estimation instead of a real spectrum. This only applies to the case \code{method = "acf"}.}
\item{scalefn}{function which calculates a scale estimator. Only relevant if \code{var1 = FALSE}. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
\item{na.action}{function to be called to handle missing values.}
\item{series}{the time series name.}
\item{...}{further arguments passed to the acf calculation by \code{\link{acfrob}} if \code{method == "acf"} or the AR calculation if \code{method == "pgram"}.}
}

\details{
There are two general methods to compute a robust spectrum.

If \code{method == "pgram"} the pre-whitening method as described in Chapter 8.14 of Maronna et al. (2006) is applied. The method combines two spectrum estimations. First one fits an AR model robustly. See the help page of \code{\link{arrob}} for an overview of possible methods which can be chosen by the argument \code{arrob.method}. Afterwards residuals of this fit are calculated and transformed by a given psi function \code{psifn}, which should project exceptional large values to 0. The usual spectrum of this residuals is then combined with the theoretical spectrum of the AR fit, see Maronna et al. (2006) for details.

If \code{method == "ar"}, the spectrum is calculated directly from acf using the well-known relationship between acf and spectrum. The acf can be estimated robustly by various methods, choosing a specific approach by the argument \code{acf.approach}. See the help page of \code{\link{acfrob}} for more details and additional arguments. Because of computational reasons, it is usually not possible to estimate the whole acf up to lag \eqn{n-1}. Truncating the acf will always yield a smoothed estimate. However, this is desired anyway in most applications.

Since the autocorrelation does not contain information about the scale of the time series, one can only derive the spectral density, which integrates to 1. If one is interested in the classical spectrum, one has to set \code{var1 = FALSE} (the default). In this case the estimated spectral density is transformed to a spectrum using a scale estimator which can be specified using the argument \code{scalefn}.

The degree of smoothing can be either determined by setting either a bandwidth \code{bandwith} or a maximal lag of the acf \code{truncation}. Note that a given bandwidth overrides any value given by the argument \code{truncation}.
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
spectrumrob(tss, method = "pgram", plot = FALSE)
spectrumrob(tss, method = "acf", plot = FALSE)
}
