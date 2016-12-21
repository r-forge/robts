\encoding{utf8}
\name{M_estimation}
\alias{M_psi}
\alias{M_wgt}
\alias{smoothpsi}

\title{
Psi and Weight Functions for Robust M Estimation
}
\description{
These psi and corresponding weight functions are employed for robust M estimation and robust filtering.

These functions are intended for internal use and are called for example by the functions \code{\link{acfrob.RA}}, \code{\link{arrob.filter}}, \code{\link{filterrob}} and \code{\link{spectrumrob.pgram}}.
}

\usage{
M_psi(x, type = c("huber", "bisquare", "smooth"), k)

M_wgt(x, type = c("huber", "bisquare"), k)

smoothpsi(x, k = c(2, 3))
}

\arguments{
	\item{x}{numeric for which the respective function should be returned. This can either be a single value or a vector.}
	\item{type}{character string giving the type of psi respectively weight function of interest.}
	\item{k}{numeric giving the tuning parameter(s).}
}

\details{
The definition of these functions can be accessed by typing \code{M_psi}, \code{M_wgt} or \code{smoothpsi}.
}  

\value{
Numeric with the specified psi respectively weight function evaluated at \code{x}.
}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}

\examples{
curve(M_psi(x, type = "huber", k = 2), from = -4, to = 4,
  main = expression(paste("Huber's ", Psi, " function with tuning constant 2")))

curve(M_psi(x, type = "bisquare", k = 3), from = -4, to = 4,
  main = expression(paste("Tukey's ", Psi, " function with tuning constant 3")))

curve(M_psi(x, type = "smooth", k = c(2, 3)), from = -4, to = 4,
  main = expression(paste("Smooth ", Psi, " function with tuning constants 2 and 3")))
}
