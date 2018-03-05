\encoding{utf8}
\name{summary.arrob}
\alias{summary.arrob}
\alias{print.summary.arrob}

\title{
Summary of a Fitted Autoregressive Model
}
\description{
Summary method for class \code{"arrob"}.}
\usage{

## S3 method for class \code{arrob}
\method{summary}{arrob}(object, correlation=FALSE, symbolic.cor=FALSE, ...)

## S3 method for class \code{summary.arrob}
\method{print}{summary.arrob}(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...)
}

\arguments{
	\item{object}{an object of class \code{"arrob"}, usually, a result of a call to \code{arrob}.}
	\item{x}{an object of class \code{"summary.arrob"}, usually, a result of a call to \code{summary.arrob}.}
	\item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated AR parameters is returned and printed.}
	\item{digits}{the number of significant digits to use when printing.}
	\item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form (see \code{symnum}) rather than as numbers.}
	\item{signif.stars}{logical. If \code{TRUE}, ‘significance stars’ are printed for each coefficient.}
    \item{...}{further arguments are currently ignored. Only for compatibility with generic function.}
}

\details{
Standard errors and correlation matrix of the coefficients can only be calculated if one sets \code{asyvar=TRUE} in the function \code{arrob}. In this case they are determined by a block bootstrap. Theoretical results to the consistency of this method are missing yet, so these values should be taken with care just as the printed p-values which are based on the asymptotic normality of the estimators which is also not formally proofed yet. 
}

\value{
Object of classes \code{"summary.arrob"}. This is a list consisting of:

\item{residuals}{numeric vector of estimated residuals.}

\item{method}{character string indicating the used eestimation method.}

\item{order}{numeric value giving the order of the fitted autoregressive model.}

\item{estsig}{numerical value giving the estimated scsle of the innovations.}

\item{mean}{numerical value giving the estimated location of the time series.}

\item{coefvar}{matrix of the estimated covariance of the AR parameters.}

\item{n}{numerical value giving the number of observations.}

\item{corp}{logical indicating whether correlation between the estimaed paramters should be determined.}

\item{cor}{matrix of estimated correlation betweeb the AR parameters.}

\item{symbolic.cor}{logical indicating whether \code{cor} should be printed as symbols.}


}

\references{
Lahiri, S.N. (1999): Theoretical Comparisons of Block Bootstrap Methods, \emph{The annals of Statistics}, vol. 27, 386--404. 
}

\author{
Alexander \enc{Dürre}{Duerre}
}

\seealso{
Robust estimation of autoregressive processes \code{\link{arrob}}.
}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 0.5), n = 100)
summary(arrob(tss, method = "yw"))
}

\keyword{Autoregressive models}
