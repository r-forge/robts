\encoding{utf8}
\name{armarob}
\alias{armarob}


\title{
Robust Fitting of ARMA Models
}
\description{
	Robust fit of an autoregressive moving average model to a time series. Its order is either optimized with respect to a robustified AIC or specified by the user.
}
\usage{
arrob(x, arorder, maorder, aic = FALSE, aicpenalty = function(p) {2*p}, na.action = na.fail, series = deparse(substitute(x)))

\arguments{
	\item{x}{numeric vector of a univariate time series.}
	\item{arorder}{integer value determining the (maximal) order of the AR part.}
	\item{maorder}{integer value determining the (maximal) order of the MA part.}
\item{aic}{logical indicating whether the ARMA order should be estimated by robust AIC criteria considering orders up to \code{arorder} respectively \code{maorder} (if \code{aic = TRUE}) or should be equal to \code{arorder} respectively \code{maorder} (if \code{aic = FALSE}). Note that fitting the model using the AIC criterion can take some time if the order if the AR or MA part is large.}
	\item{aicpenalty}{function of the model order defining the penalty term of the model selection criteria. The default results in a robust AIC.}
	\item{na.action}{function to be called to handle missing values. Default is \code{\link{na.fail}}. Other reasonable possibilities are \code{\link{na.omit}} (omits all missing values at the beginning and the end of the time series and fails if there are more missing values), \code{\link{na.contiguous}} (uses longest contiguous stretch of non-missing observations) and \code{\link{na.extremify}} (experimental, sets missing values to very extreme values in the expectation that they have only little influence on the outcome of robust estimation methods).}
	\item{series}{the time series name.}
}

\details{
The model fit follows the filtering approach described in Maronna et al. (2006). First an autoregressive process is fitted using the AIC criterion. From this the estimated variance of the innovations and the estimated AR coefficients are used to determine start values for the optimization of a robustified likelihood, which yields the final parameter estimation.

If \code{aic = TRUE}, every combination of AR parts of order \code{p} and MA parts of order \code{q} up to \code{arorder} respectively \code{maorder} are estimated and the one with the smallest model selection criterion returned. The criterion is the sum of the logarithm of a robust variance estimation of the residuals and a penalty term \eqn{pen(p)} divided by \eqn{n-p} depending on the AR order \eqn{p}, i.e.,
\deqn{log(\hat{\sigma}^2)+\frac{pen(p)}{n-p}.}{log(sigma_hat^2) + (pen(p))/(n-p)}
The penalty term can be set by the argument \code{aicpenalty}.

}

\value{
Object of class \code{"armarob"} containing:

    \item{order}{numeric vector, first entry represents the AR part, second the MA part.}
    \item{ar}{numeric vector of the estimated AR parameters.}
    \item{ma}{numeric vector of the estimated MA parameters.}
    \item{var.pred}{estimated variance of the innovations.}
    \item{aic}{either the robustified aic (if \code{aic = FALSE}) or a table aic values, containing models with AR parts up to order \code{arorder} and MA parts up to order \code{maorder}.}
}

\references{

Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8, \doi{10.1002/0470010940}.

}

\author{
Alexander \enc{Dürre}{Duerre}
}

\seealso{
Classical, nonrobust fitting is provided by the function \code{\link{ar}}.

S3 methods \code{\link{residuals.arrob}}, \code{\link{fitted.arrob}}, \code{\link{filtered.arrob}} and \code{\link{predict.arrob}} for objects of class \code{"arrob"}. 
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 0.5), n = 100)
armarob(tss, arorder=2, maorder=2, aic=TRUE)

}

