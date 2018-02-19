\encoding{utf8}
\name{acfrob.bireg}
\alias{acfrob.bireg}

\title{
Robust Autocorrelation Estimation Based on Bivariate Regression
}

\description{
Robustly estimates the autocorrelation function of a time series based on robust regression of covariates with the lagged observations. See Chang and Politis (2016) for details.

This function is intended for internal usage only. Users should rather use the wrapper function \code{\link{acfrob}} with argument \code{approach = "bireg"}.
}

\usage{
acfrob.bireg(x, lag.max, regression.method = c("MM", "LTS", "L1"))
}

\arguments{
	\item{x}{univariate numeric vector or time series object.}
	\item{lag.max}{integer value giving the maximum lag at which to calculate the acf.}
	\item{regression.method}{character string indicating the regression method, see Details.}
}


\details{
This function estimates the autocorrelation function based on robust regression estimators, which was proposed by Chang and Politis (2016). The acf is estimated for each lag individually. Let \eqn{Y} denote the original time series and \eqn{X} the lagged one, then the autocorrelation is estimated by regressing \eqn{Y} on \eqn{X}. Three robust regression estimators are proposed in Chang and Politis (2016) which are all implemented.

\describe{
\item{\code{MM}}{MM estimators were proposed in Yohai 1987 and are known to be robust and highly efficient under normality. See the help page of \code{\link{lmrob}} in the package \code{robustbase} for details.}

\item{\code{LTS}}{The least trimmed squares estimator is proposed by Rousseeuw 1984 and is highly robust but inefficient under the recommended trimming amount, see Maronna et al. (2016). See the help page of \code{\link{ltsReg}} in the package \code{robustbase} for details.}

\item{\code{L1}}{The median regression minimizes the sum of the absolute values of the residuals. It is not robust against outliers in the regressor space (see Maronna et al. 2016) which in the autocorrelation context is the timeseries itself. Therefore this method should be only applied if one is sure that there are no large outliers.}
}
}

\value{
A named list of the following elements:

\item{acfvalues}{Numeric vector of estimated autocorrelations at the lags 1,...,\code{lag.max}.}
\item{are}{numeric value giving the asymptotic relative efficiency (ARE) of the estimator as compared to the classical nonrobust estimator, under the assumption that the observations are uncorrelated and from a Gaussian distribution. The ARE is currently not available for this estimation approach and is therefore \code{NA}.}
} 

\references{

Chang, C., Politis, D. (2016): Robust autocorrelation estimation, \emph{Journal of Computational and Graphical Statistics}, vol. 25, 144--166, \doi{10.1080/10618600.2014.969431}.

Maronna, R., Martin, D., Yohai, V. (2006): Robust Statistics - Theory and Methods.

Rousseeuw, P. (1984): Least Median of Squares Regression, \emph{Journal of the American Statistical Association}, vol. 79, 871--880, \doi{10.1080/01621459.1984.10477105} .

Yohai, V. (1987): High Breakdown-Point and High Efficiency Robust Estimates for Regression, \emph{The Annals of Statistics}, vol. 15, 642--656, \url{http://www.jstor.org/stable/2241331}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{
The wrapper function \code{\link{acfrob}}.

Alternative acf subroutines: \code{\link{acfrob.GK}}, \code{\link{acfrob.filter}}, \code{\link{acfrob.multi}}, \code{\link{acfrob.partrank}}, \code{\link{acfrob.RA}}, \code{\link{acfrob.rank}}, \code{\link{acfrob.trim}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 0.5), n = 100)
acfrob(tss, approach = "bireg")
}

\keyword{Autocorrelation function}
