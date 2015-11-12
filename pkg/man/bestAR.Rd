\encoding{utf8}
\name{bestAR}
\alias{bestAR}

\title{
Robust Estimation of Autoregressive Models by GM Estimator
}
\description{
	Robust fit of an autogression model to a time series by generalized M estimators as described in Maronna et al. (2006).
}
\usage{
bestAR(timeseries,maxp,maxit=10^3,delta=1/2,epsilon=10^(-4),
    k1=1.37,k2=1,aicpenalty=function(p) {return(2*p)})
}

\arguments{
	\item{timeseries}{numeric vector of a univariate time series.}
	\item{maxp}{numeric value indicating the maximum AR model order considered.}
	\item{maxit}{numeric value indicating the maximal number of iterations to determine the GM estimator, see details.}
	\item{delta}{numeric value indicating the breakpoint of the GM estimator. Values larger than 0.5 results in an implosion breakpoint of 1-\code{delta}, see details.}
	\item{epsilon}{numeric value indicating the stopping rule of the iteration algorithm to determine the GM estimator, see details.}
	\item{k1}{numeric defining the tuning parameter of the huber weights, see details.}
	\item{k2}{numeric defining the tuning parameter of the bisquare weights of the regressors, see details.}
	\item{aicpenalty}{function of the model order defining the penalty term of the model selection criteria. The default results in a robust AIC.}
}

\details{

This procedure fits an AR model to a timeseries. The AR coefficients are determined by robust regression, using the timeseries as dependent variable and the lagged timeseries as independent one. So if there are outliers in the timeseries, they are also in the regressors. That is why generalized M estimators are proposed which downweight also the explanatory variables. To minimize the complexity of the problem, the Durbin-Levinson algorithm is applied, fitting a sequence of AR processes of increasing order, which break the estimation down to \code{maxp} regressions with one prediction variable instead of one regression with \code{maxp} predictors.

The GM estimator is of Mallows type with a Huber function for the regression residuals and a bisquare function for the regressors. Both can be tuned by \code{k1} respectively \code{k2}. The weights \eqn{v_i} of the regressors are based on malahanobis distances with the covariance matrix build of autocovariances derived by the so far estimated AR parameters, see Maronna et al. (2006) chapter 8.5 for details. 

The GM estimator is computed as an iteratively weighted sum. Based on an initial estimator, which is here the least trimmed squares estimator \code{\link{ltsReg}} of the \code{robustbase} package, residuals \eqn{r_i} and its related huber \eqn{w_i} and bisquare weights \eqn{u_i} are calculated. The updated regression estimate is then
\deqn{\beta=\frac{\sum_{i=1}^n w_iv_i x_iy_i}{\sum_{i=1}^n w_iv_i x_i^2}}
and the scale estimate
\deqn{\sigma^2_{new}=\frac{\sigma^2_{old}}{n \delta}\sum_{i=1}^n u_i r_i^2,}
where \eqn{\delta} determines the breakpoint of the scale estimator. For \eqn{\delta<1/2} it is equivalent to the breakpoint (in fact the explosion breakpoint) and for \eqn{\delta>1/2} the (implosion) breakpoint equals \eqn{1-\delta}.

The iteration stops if either the maximal change of the residuals is smaller than \code{epsilon} or the maximal number of iterations \code{maxit} is reached. Subsequently the scale estimate is transformed to be consistent under normality.

}

\value{
List containing the following values

\item{phima}{matrix representing the successive estimated AR coefficients of the model. In the first column are the coefficients of the AR model or order 1 and so on.}
\item{aic}{numeric of aic values of successive AR modells.}
\item{var.pred}{numeric containing the successive variance estimations of the residuals of the successive AR models. The first entry is for the AR(0) model and so on.}
\item{x.mean}{numeric of estimated location of the timeseries.}
\item{residuals}{matrix containing the residuals of successive AR models. The first column corresponds to the AR(1) model and so on.}
}

\references{

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8.

}
\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
The wrapper function \code{\link{arrob}}.
}
\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
bestAR(tss, maxp=5)
}
