\encoding{utf8}
\name{corGK}
\alias{corGK}

\title{
Robust Correlation Estimation Based on the Gnanadesikan Kettenring Approach
}
\description{
Robustly estimates the correlation of two samples based on an Gnanadesikan Kettenring (GK) estimator.

This function is intended for internal usage only and is called by the function \code{\link{acfrob.GK}}.
}

\usage{
corGK(x, y, scalefn = Qn, ...)
}

\arguments{
	\item{x, y}{numeric vectors of the same length with the two samples.}
	\item{scalefn}{function which calculates a scale estimator. Its first argument must be a vector of the data. One could for example use one of the functions \code{\link{Qn}} (the default), \code{\link{Sn}}, \code{\link{scaleTau2}} (all from the package \code{robustbase}), \code{\link{mad}} or (not robust) \code{\link{sd}}.}
	\item{...}{further arguments passed to the function \code{scalefn}.}
}

\details{
The definition of this functions can be accessed by typing \code{corGK}.
}  

\value{
Numeric value with the estimated correlation.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}
