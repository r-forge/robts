\encoding{utf8}
\name{arrob.gm}
\alias{arrob.gm}

\title{
Robust Fitting of Autoregressive Models Based on GM Estimation
}
\description{
Robust fit of an autoregressive model to a time series by generalized M estimators as described in Maronna et al. (2006).
}
\usage{
arrob.gm(x, order.max = NULL, aic = TRUE, aicpenalty=function(p) 2*p,
    na.action = na.fail, series = deparse(substitute(x)), k1=1.37, k2=1,
    delta=1/2, maxit=10^3, epsilon=10^(-4))        
}

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{aic}{logical indicating whether the AR order should be estimated by robust AIC criteria considering orders up to \code{order.max} or should be equal to \code{order.max}.}
	\item{order.max}{integer value determining the (maximal) order of the AR fit.}
	\item{na.action}{function to be called to handle missing values.}
	\item{series}{the time series name.}
	\item{aicpenalty}{function of the model order defining the penalty term of the model selection criteria. The default results in a robust AIC.}
  \item{k1}{numeric value giving the tuning parameter for the Huber weights.}
  \item{k2}{numeric value giving the tuning parameter for the psi function applied to the regressors.}
  \item{delta}{numeric value for \eqn{\delta}, see Details.}
  \item{maxit}{integer value giving rhe maximal number of iterations for the iterative weighting procedure of the M estimators}
  \item{epsilon}{numeric value giving the desired accuracy of the iterated solution of the M estimator.}
}

\details{
This procedure fits an AR model to a time series. The AR coefficients are determined by robust regression, using the time series as dependent variable and the lagged time series as independent one. So if there are outliers in the time series, they are also in the regressors. That is why generalized M estimators are proposed which downweight also the explanatory variables. To minimize the complexity of the problem, the Durbin-Levinson algorithm is applied, fitting a sequence of AR processes of increasing order, which break the estimation down to \code{maxp} regressions with one prediction variable instead of one regression with \code{maxp} predictors.

The GM estimator is of Mallows type with a Huber function for the regression residuals and a bisquare function for the regressors. Both can be tuned by \code{k1} respectively \code{k2}. The weights \eqn{v_i} of the regressors are based on Malahanobis distances with the covariance matrix build of autocovariances derived by the so far estimated AR parameters, see Maronna et al. (2006) chapter 8.5 for details. 

The GM estimator is computed as an iteratively weighted sum. Based on an initial estimator, which is here the least trimmed squares estimator \code{\link{ltsReg}} of the \code{robustbase} package, residuals \eqn{r_i} and its related Huber \eqn{w_i} and bisquare weights \eqn{u_i} are calculated. The updated regression estimate is then
\deqn{\beta=\frac{\sum_{i=1}^n w_iv_i x_iy_i}{\sum_{i=1}^n w_iv_i x_i^2}}
and the scale estimate
\deqn{\sigma^2_{new}=\frac{\sigma^2_{old}}{n \delta}\sum_{i=1}^n u_i r_i^2,}
where \eqn{\delta} determines the breakdown point of the scale estimator. For \eqn{\delta<1/2} it is equivalent to the breakpoint (in fact the explosion breakdown point) and for \eqn{\delta>1/2} the (implosion) breakdown point equals \eqn{1-\delta}.

The iteration stops if either the maximal change of the residuals is smaller than \code{epsilon} or the maximal number of iterations \code{maxit} is reached. Subsequently the scale estimate is transformed to be consistent under normality.
}

\value{
An object of class \code{"ar"}, see \code{\link{ar}} for details.

Note that list element \code{aic} gives the value of the robust information criterion and not its difference with the lowest information criterion of all considered models, as it is returned by the function \code{\link{ar}}.
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
The wrapper function \code{\link{arrob}}.

Classical, nonrobust fitting is provided by the function \code{\link{ar}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
arrob(tss, method="gm")
}

\keyword{Autoregressive models}
