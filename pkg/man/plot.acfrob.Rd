\encoding{utf8}
\name{plot.acfrob}
\alias{plot.acfrob}

\title{
Plot Robust Autocorrelation, Autocovariance or Partial Autocorrelation Functions
}

\description{
Plot method for objects of class \code{"acfrob"}.
}

\usage{
\method{plot}{acfrob}(x, ci = 0.95, ci.col = "blue", ci.type = "white",
           ylim = NULL, ...)
}

\arguments{
	\item{x}{an object with classes \code{"acfrob"} and \code{"acf"}.}
	\item{ci}{coverage probability for confidence interval. Plotting of the confidence interval is suppressed if \code{ci} is zero or negative.}
	\item{ci.col}{colour to plot the confidence interval lines.}
	\item{ci.type}{character string giving the type of confidence limits. Currently the only possible type is \code{"white"}, where a white noise input is assumed.}
  \item{ylim}{numeric vector of length two giving the vertical limits for the plot.}
	\item{...}{further arguments passed to \code{\link{plot.acf}}. See that help page for details.}
}

\details{
Most of the plotting is done by the function \code{\link{plot.acf}}. Only the conficence limits of a robustly fitted (partial) autocorrelation function need some extra care.

The confidence interval plotted in \code{plot.acfrob} is based on an \emph{uncorrelated} series and should be treated with appropriate caution. 
Plotting confidence limits is currently only available for some of the approaches of the estimation function \code{\link{acfrob}}, for which the element \code{ci.factor} of their output is not \code{NA}. Otherwise plotting of the confidence limits are omitted.  
For auto\emph{covariance} functions there are currently no confidence limits available.
}


\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\seealso{
\code{\link{acfrob}} which calls \code{plot.acfrob} by default. Classical, nonrobust estimation is provided by the function \code{\link{acf}}.
}

\keyword{Autocorrelation function}
