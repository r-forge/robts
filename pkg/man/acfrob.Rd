\encoding{utf8}
\name{acfrob}
\alias{acfrob}

\title{
Robust Autocorrelation, Autocovariance or Partial Autocorrelation Function Estimation
}

\description{
Robustly estimates (and by default plots) the autocorrelation function (acf), the autocovariance function (acvf) or the partial autocorrelation function (pacf) of a time series. This is a wrapper function calling one of various subroutines. See \enc{Dürre}{Duerre} et al. (2015) for details.
}

\usage{
acfrob(x, lag.max = NULL, type = c("correlation", "covariance", "partial"),
  approach = c("acfGK", "acfmedian", "acfmulti", "acfpartrank", "acfRA",
  "acfrank", "acfrobfil", "acftrim", "acfregression"), plot = TRUE,
  na.action = na.fail, psd = TRUE, scalefn = Qn,
  partial.method = c("automatic", "durbin-levinson"), ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf. Default ist \code{10 * log(length(x), 10)}.}
	\item{type}{character string giving the type of acf to be computed. Allowed values are \code{"correlation"} (the default), \code{"covariance"} or \code{"partial"}.}
	\item{approach}{character string naming the approach used for computation, see Details.}
	\item{plot}{logical. If \code{TRUE} (the default) the acf is plotted.}
	\item{na.action}{function to be called to handle missing values.}
	\item{scalefn}{function which calculates a scale estimator. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
	\item{partial.method}{character string giving the method employed for calculation of the acf in case \code{type == "partial"}, see Details.}
	\item{psd}{logical. If \code{TRUE} the autocorrelation is enforced to be positive semidefinite, see Details.}
	\item{...}{further arguments passed to the internal function of the respective estimation approach.}
}

\details{
There are many approaches for robust autocorrelation estimation. They are described in the following.
	
	If \code{fun = "acfGK"}, one applies a Gnanadesikan-Kettenring (GK) estimator, which calculates the acf for every lag individually and is based on a scale estimator. This is method is quick and favourable if outliers are isolated. See the help page of \code{\link{acfGK}} for more details and additional arguments.
	
	If \code{fun = "acfmedian"}, one applies the so called median correlation, which calculates the acf for every lag individually. Basically this estimator substitutes the median for the mean in the definition of the empirical correlation. See the help page of \code{\link{acfmedian}} for details.
	
	If \code{fun = "acfmulti"}, one applies a robust estimator for the correlation matrix to compute the acf en bloc. This method can get quite time consuming if \code{lag.max} is large (e.g. larger than 15). It is favourable if outliers accour in blocks. See the help page of \code{\link{acfmulti}} for details and further arguments.
	
	If \code{fun = "acfpartrank"}, one first calculates the partial autocorrelations and then derive the acf. This procedure has the advantage that one automatically gets a positive definite acf. See the help page of \code{\link{acfpartrank}} for more details and further arguments.
	
	If \code{fun = "acfRA"}, one first robustly transform the time series using a psi function and calculates the usual acf afterwards. See the help page of \code{\link{acfRA}} for more details and further arguments.
	
	If \code{fun = "acfrank"}, one applies a rank correlation estimator like Spearmans rho, which calculates the acf for every lag individually. It is favourable if one considers integer valued processes. See the help page of \code{\link{acfrank}} for more details and further arguments.
	
	If \code{fun = "acfrobfil"}, one applies the robust filter algorithm. This is basically a robust version of the Kalman filter. The acf is then either computed by applying the usual acf to the filtered time series or by a robust AR fit. See the help page of \code{\link{acfrobfil}} for more details and further arguments.
	
	If \code{fun = "acftrim"}, one applies a trimmed correlation, which calculates the acf for every lag individually. This estimator basically substitutes the mean by a trimmed mean in the definition of the empirical correlation. See the help page of \code{\link{acftrim}} for more details and further arguments.

	If \code{fun = "acfregression"}, one applies a method proposed by Politis. WE NEED MORE DETAILS HERE.
	
For obtaining the robust autocovariance estimation (\code{type = "covariance"}), the acf is multiplied with a robust estimation of the marginal variance of the time series. This estimation is obtained by the squared scale estimator \code{scalefn}, which is by default \code{\link{Qn}} from the package \code{robustbase}.
	
The partial autocorrelation (\code{type = "partial"}) could be computed by employing several methods. If \code{partial.method = "durbin-levinson"} one first estimates the acf by the chosen estimation approach and then derives the pacf from this by applying the Durbin-Levinson algorithm. This method works for all estimation approaches.
For some approaches there exists a direct way to obtain the pacf. If \code{partialtype = "automatic"} (the default) such a direct way is used wherever it is available and otherwise the Durbin-Levinson algorithm is applied. For \code{approach = "acfpartrank"} and \code{approach = "acfrobfil"} the pacf is calculated anyway and will be used directly if \code{partialtype = "automatic"}.

With exception of the approaches \code{"acfpartrank"}, \code{"acfrobfil"}, \code{"acftrim"}, \code{"acfRA"}, a positive definite acf is not guaranteed, although this will often be the case. If one needs a positive definite and therefore valid acf, one can enforce this property by setting \code{psd = TRUE}. Usually this will also improve the acf estimation. See the help page of \code{\link{make_acf_psd}} for more details. 
}

\value{
Object of class \code{"acf"}, see \code{\link{acf}} for details.
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
Classical, nonrobust estimation is provided by the function \code{\link{acf}}.

The subroutines which actually estimate the acf: \code{\link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, plot = FALSE)
acfrob(tss, type = "partial", plot = FALSE)
acfrob(tss, type = "covariance", plot = FALSE)
}
