\encoding{utf8}
\name{na.extremify}
\alias{na.extremify}
\alias{na.extremify.ts}

\title{
Handing Missing Values in Objects for Application of Robust Procedures
}

\description{
This generic function is useful for dealing with \code{\link{NA}}s. It returns the object and replaces missing values by a vale which is very large with respect to the range of values in that object. This can be reasonable if a robust estimation procedure is applied which is influenced not at all or only very little by extreme observations.
}

\usage{
na.extremify(object, ...)

\method{na.extremify}{ts}(object, ...)
}

\arguments{
	\item{object}{an object. Currently there is only a methods for numeric vectors.}
	\item{...}{currently ignored.}
}

\details{
Missing values are replaced by the value \code{max(object) + 1e+03 * diff(range(object))} rounded up to the next power of 10.

This approach is only sensible if the proportion of missing values is lower than the breakdown point of the robust procedure which will be applied to the data afterwards. Note that all missing values are replaced by the same large value. The induced discreteness and asymmetry might cause problems for some robust estimation procedures.
}

\value{
The object without missing values. The class of \code{object} will be preserved. If there are missing values in \code{object}, the output has attribute \code{"na.action"} which provides information on the missing values and how these have been handled.
}

\author{
Tobias Liboschik
}

\seealso{
\code{\link{na.action}}; \code{\link{options}} with argument \code{na.action} for setting NA actions.

\code{\link{na.fail}}, \code{\link{na.omit}}, \code{\link{na.exclude}}, \code{\link{na.pass}} and \code{\link{na.contiguous}} for other methods for handling missing values.
}

\examples{
x <- ts(c(1, 2, 3, 4, NA, 6, 7))
na.extremify(x)
median(na.extremify(x))
}
