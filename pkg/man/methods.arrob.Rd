\encoding{utf8}
\name{methods.arrob}
\alias{residuals.arrob}
\alias{fitted.arrob}
\alias{filtered}
\alias{filtered.arrob}
\alias{predict.arrob}

\title{
S3 Methods for Robustly Fitted Autoregressive Models
}
\description{
	Methods for objects of class \code{"arrob"} as it is returned by a call to the function \code{\link{arrob}}. There are methods to extract (robustly filtered) residuals, fitted values or the robustly filtered observations themselves. There is also a method for (robust) prediction. 
}
\usage{
\method{residuals}{arrob}(object, x = object$x,
                method = c("recursive", "statespace", "nonrobust"),
                na.action = na.fail, ...)

\method{fitted}{arrob}(object, x = object$x,
             method = c("recursive", "statespace", "nonrobust"),
             na.action = na.fail, ...)

filtered(object, ...)
\method{filtered}{arrob}(object, x = object$x,
               method = c("recursive", "statespace", "nonrobust"), na.action, ...)

\method{predict}{arrob}(object, newdata, n.ahead = 1, se.fit = TRUE,
              method = c("recursive", "statespace", "nonrobust"), ...)
}

\arguments{
	\item{object}{an object of class \code{"arrob"}.}
	\item{x}{numeric vector of a univariate time series.}
	\item{method}{character string naming the robust filtering method to be used, see Details.}
	\item{newdata}{numeric vector of a univariate time series of which future values sould be predicted.}
	\item{n.ahead}{integer value giving the number of steps ahead at which to predict.}
	\item{se.fit}{logical value. Should the estimated standard errors of the prediction be returned?}
	\item{na.action}{function to be called to handle missing values.}
	\item{...}{additional arguments passed to the function \code{\link{filterrob}}.}	
}

\details{
If \code{method = "recursive"} or \code{method = "statespace"}, the observations in argument \code{x} (respectively \code{newdata}) are robustly filtered. See the function \code{\link{filterrob}} for details.
}

\value{
Numeric vector of residuals, fitted values or filtered observations.

For \code{predict.arrob}, a time series of predictions, or if \code{se.fit = TRUE}, a list with components \code{pred}, the predictions, and \code{se}, the estimated standard errors.
}

\references{
Maronna, R. A., Martin, R. D., and Yohai, V. J. (2006): Robust Statistics: Theory and Methods, Wiley, chapter 8, \doi{10.1002/0470010940}.
}

\author{
Alexander \enc{Dürre}{Duerre} and Tobias Liboschik
}

\seealso{
Robust fitting of autoregressive models with \code{\link{arrob}}.

Plotting method \code{\link{plot.arrob}} for fitted autoregressive models.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
fit <- arrob(tss, method = "yw")
residuals(fit)
fitted(fit)
filtered(fit)
predict(fit)
}

\keyword{Autoregressive models}
