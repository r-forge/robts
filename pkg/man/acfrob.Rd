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
    partialtype = c("ranks", "durbin-levinson", "filter"),
    posdef = TRUE, ...)
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{numeric value of maximum lag at which to calculate the acf.  If \code{NULL} (default), \code{10 * log(length(x), 10)}, is chosen.}
	\item{type}{character string whether the partial autocorrelation is calculated.}
	\item{fun}{character string naming the function used for computation, see details.}
	\item{plot}{logical.  If \code{TRUE} (the default) the acf is plotted.}
	\item{na.action}{function to be called to handle missing values.}
	\item{p}{numeric value of order of AR fit used if \code{fun == acfrobfil}, see details.}
	\item{robfiltype}{character string naming the type of robust filtering if \code{fun == "acfrobfil"}, see details.}
	\item{partialtype}{character string naming the type of pacf calculation if \code{type == "partial"}, see details.}
	\item{posdef}{logical whether autocorrelation is enforced to be positive semidefinite, see details.}
	\item{...}{further arguments passed to the respective internal autocorrelation function.}
}

\value{
Object of class \code{"acf"}, see \code{\link{acf}} for details.
}

\details{

This function robustly estimates the (partial) autocorrelation function. If one is interested in a robust autocovariance estimation, one can multiply the acf with a robust estimation of the marginal variance of the timeseries, e.g. using the \code{\link{Qn}} of the \code{robustbase} package. There are many approaches for robust autocorrelation estimation. They are described in the following.
	
	If \code{fun="acfGK"} one applies a Gnanadesikan-Kettenring (GK) estimator, which calculates the acf for every lag individually and is based on a scale estimator. This is method is quick and favourable if autliers are isolated.  See the help page of \code{link{acfGK}} for more details and additional arguments.
	
	If \code{fun="acfmedian"} one applies the so called mediancorrelation, which calculates the acf for every lag individually No further arguments. Basically this estimator substitutes the median for the mean in the definition of the empirical correlation. See the help of \code{\link{acfmedian}} for details.
	
	If \code{fun="acfmulti"} one applies a robust estimator for the correlation matrix to compute the acf en bloc. This method can get quite time consuming if lag.max is large (e.g. larger than 15). It is favourable if outliers accour in blocks. See \code{\link{acfmulti}} for details and further arguments.
	
	If \code{fun="acfpartrank"} one first calculates the partial autokorrelations and then derive the acf. This procedure has the advantage that one automatically gets a positive definite acf. See the help of \code{\link{acfpartrank}} for more details and further arguments.
	
	If \code{fun="acfRA"} one first robustly transform the timeseries using a psi-function and calculates the usual acf afterwards. See the help of \code{\link{acfRA}} for more details and further arguments.
	
	If \code{fun="acfrank"} one applies a rank correlation estimator like Spearmans-tau, which calculates the acf for every lag individually. It is favourable if one considers integer valued processes. See the help of \code{\link{acfrank}} for more details and further options.
	
	If \code{fun="acfrobfil"} one applies the robust filter algorithm. This is basically a robust version of the Kalman filter. The acf is then either computed by applying the usual acf to the filtered timeseries or by a robust AR-fit. See the help of \code{\link{acfrobfil}} for more details and further arguments.
	
	If \code{fun="acftrim"} one applies a trimmed correlation, which calculates the acf for every lag individually. This estimator basically substitutes a trimmed mean for the mean in the definition of the empirical correlation.  See the help of \code{\link{acftrim}} for more details and further arguments.
	
The partial autocorrelation could be either computed by first estimating the acf by one of the above methods (by setting the related \code{fun} argument and \code{partialtype="durbin-levinson"}) and deriving the pacf by the Durbin-Levinson algorithm or by one of the following methods. 
	
	If \code{partialtype="ranks"} one uses a rank based correlation estimator. See the help of \code{\link{pacf2}} for more information and further arguments.
	
	If \code{partialtype="filter"} one uses the robust filter algorithm. This is basically a robust version of the Kalman filter and provides a robust AR fit which is then used to derive the partial autocorrelation. See the help of \code{\link{ARfilter}} for more details.

	With exception of the methods \code{fun=c("acfpartrank","acfrobfil","acftrim","acfRA")} a positive definite acf is not guaranteed, although it will be in most of the cases. If one needs a positive definite and therefore valid acf one can enforce this property by a projection algorithm. Usually this will also improve the acf estimation. See the help of \code{\link{acfposmaker}} for more details. 
}

\references{
\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
Classical, nonrobust estimation is provided by the function \code{\link{acf}}.

The subroutines which actually estimates the acf: \code{link{acfGK}}, \code{\link{acfmedian}}, \code{\link{acfmulti}}, \code{\link{acfpartrank}}, \code{\link{acfRA}}, \code{\link{acfrank}}, \code{\link{acfrobfil}}, \code{\link{acftrim}}.

The subroutines which actually estimates the pacf: \code{\link{pacf2}}, \code{\link{ARfilter}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
acfrob(tss, plot = FALSE)
acfrob(tss, type = "partial", plot = FALSE)
}
