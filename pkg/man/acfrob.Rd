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
  approach = c("GK", "median", "multi", "partrank", "RA", "rank", "filter",
  "trim", "bireg"), ..., plot = TRUE, na.action = na.fail, psd = TRUE,
  scalefn = Qn, partial.method = c("automatic", "durbin-levinson"))
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf. Default is \code{10 * log(length(x), 10)}.}
	\item{type}{character string giving the type of acf to be computed. Possible values are \code{"correlation"} (the default), \code{"covariance"} or \code{"partial"}.}
	\item{approach}{character string naming the approach used for computation, see Details.}
	\item{...}{further arguments passed to the internal function of the respective estimation approach.}
	\item{plot}{logical. If \code{TRUE} (the default) the acf is plotted.}
	\item{na.action}{function to be called to handle missing values. Default is \code{\link{na.fail}}. Other reasonable possibilities are \code{\link{na.omit}} (omits all missing values at the beginning and the end of the time series and fails if there are more missing values), \code{\link{na.contiguous}} (uses longest contiguous stretch of non-missing observations) and \code{\link{na.extremify}} (experimental, sets missing values to very extreme values in the expectation that they have only little influence on the outcome of robust estimation methods).}
	\item{psd}{logical. If \code{TRUE} the autocorrelation is enforced to be positive semidefinite, see Details.}
	\item{scalefn}{function which calculates a scale estimator. Its argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
	\item{partial.method}{character string giving the method employed for calculation of the acf in case \code{type == "partial"}, see Details.}
}

\details{
There are many approaches for robust autocorrelation estimation. The argument \code{approach} specifies which of the following approaches is to be used, where \code{approach = "GK"} is the default.
\describe{	
\item{\code{"GK"}}{One applies a Gnanadesikan-Kettenring (GK) estimator, which calculates the acf for every lag individually and is based on a scale estimator. This is method is quick and favourable if outliers are isolated. See the help page of \code{\link{acfrob.GK}} for more details and additional arguments.}
	
\item{\code{"median"}}{One applies the so called median correlation, which calculates the acf for every lag individually. Basically this estimator substitutes the median for the mean in the definition of the empirical correlation. See the help page of \code{\link{acfrob.median}} for details.}
	
\item{\code{"multi"}}{One applies a robust estimator for the correlation matrix to compute the acf en bloc. This method can get quite time consuming if \code{lag.max} is large (e.g. larger than 15). It is favourable if outliers occur in blocks. See the help page of \code{\link{acfrob.multi}} for details and further arguments.}
	
\item{\code{"partrank"}}{One first calculates the partial autocorrelations and then derive the acf. This procedure has the advantage that one automatically gets a positive definite acf. See the help page of \code{\link{acfrob.partrank}} for more details and further arguments.}
	
\item{\code{"RA"}}{One first robustly transform the time series using a psi function and calculates the usual acf afterwards. See the help page of \code{\link{acfrob.RA}} for more details and further arguments.}
	
\item{\code{"rank"}}{One applies a rank correlation estimator like Spearman's rho, which calculates the acf for every lag individually. It is favourable if one considers integer valued processes. See the help page of \code{\link{acfrob.rank}} for more details and further arguments.}
	
\item{\code{"filter"}}{One applies the robust filter algorithm. This is basically a robust version of the Kalman filter. The acf is then either computed by applying the usual acf to the filtered time series or by a robust AR fit. See the help page of \code{\link{acfrob.filter}} for more details and further arguments.}
	
\item{\code{"trim"}}{One applies a trimmed correlation, which calculates the acf for every lag individually. This estimator basically substitutes the mean by a trimmed mean in the definition of the empirical correlation. See the help page of \code{\link{acfrob.trim}} for more details and further arguments.}

\item{\code{"bireg"}}{One applies a method proposed by Politis. HERE MORE DETAILS HAVE TO BE INCLUDED. See the help page of \code{\link{acfrob.bireg}} for more details and further arguments.}
}

For obtaining the robust autocovariance estimation (\code{type = "covariance"}), the acf is multiplied with a robust estimation of the marginal variance of the time series. This estimation is obtained by the squared scale estimator \code{scalefn}, which is by default \code{\link{Qn}} from the package \code{robustbase}.
	
The partial autocorrelation (\code{type = "partial"}) could be computed by employing several methods. If \code{partial.method = "durbin-levinson"} one first estimates the acf by the chosen estimation approach and then derives the pacf from this by applying the Durbin-Levinson algorithm. This method works for all estimation approaches.
For some approaches there exists a direct way to obtain the pacf. If \code{partial.method = "automatic"} (the default) such a direct way is used wherever it is available and otherwise the Durbin-Levinson algorithm is applied. For \code{approach = "partrank"} and \code{approach = "filter"} the pacf is calculated anyway and will be used directly if \code{partial.method = "automatic"}.

With exception of the approaches \code{"partrank"}, \code{"filter"}, \code{"trim"}, \code{"RA"}, a positive definite acf is not guaranteed, although this will often be the case. If one needs a positive definite and therefore valid acf, one can enforce this property by setting \code{psd = TRUE}. Usually this will also improve the acf estimation. See the help page of \code{\link{make_acf_psd}} for more details.
}

\value{
Object of classes \code{"acfrob"} and \code{"acf"}. This is a list which includes all elements of an object of class \code{"acf"} (see \code{\link{acf}} for details) plus the following additional elements:

\item{approach}{character string naming the estimation approach. Coincides with argument \code{approach}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. If the ARE is not available for the chosen estimation approach, it is \code{NA}.}
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2015): Robust estimation of (partial) autocorrelation, \emph{Wiley Interdisciplinary Reviews: Computational Statistics}, vol. 7, 205--222, \doi{10.1002/wics.1351}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
S3 method \code{\link[=plot.acfrob]{plot}}.

Classical, nonrobust estimation is provided by the function \code{\link{acf}}.

The subroutines which actually estimate the acf: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.median}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.bireg}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, plot = FALSE)
acfrob(tss, type = "partial", plot = FALSE)
acfrob(tss, type = "covariance", plot = FALSE)
}

\keyword{Autocorrelation function}
