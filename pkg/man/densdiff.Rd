\encoding{utf8}
\name{densdiff}
\alias{densdiff}

\title{
Density Estimation at Zero of the Random Variable X - Y
}
\description{
Estimates the density of the difference of two samples at 0. This is necessary for the Hodges-Lehmann (change point) test. See Dehling et al. (2015) for details.
}

\usage{
densdiff(x, y, type2 = c("all", "within", "between"), adjust = 1, 
    kernelused = c("gaussian", "epanechnikov", "rectangular", 
    "triangular", "biweight","cosine", "optcosine"), ...)
}

\arguments{
	\item{x}{univariate numeric vector representing the first sample.}
	\item{y}{univariate numeric vector representing the second sample.}
	\item{type2}{character string indicating which differences should be used, to estimate the density. Must be one of 'all', 'within' and 'between', see Details.}
	\item{adjust}{numeric which determines the used bandwidth, see the help of \code{\link{density}} for Details.}
	\item{kernelused}{character string indicating the used kernel, see the help of \code{\link{density}} for Details.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\details{
The Hodges-Lehmann two sample test requires an estimation of the density of the difference of both random variables X and Y at 0. This is also true for the the Hodges-Lehmann change point test since it consists of a sequence of two samples tests. One approach to estimate this density is to apply a kernel estimator at differences of values of both samples. This is basically the task of this function.

There are different possibilities to build differences:

If \code{type2} is \code{"all"} one uses all possible differences of the merged sample.

If \code{type2} is \code{"within"} one uses all possible differences within the first sample and all possible differences within the second.

If \code{type2} is \code{"between"} one uses all possible differences between one observation of the first sample and one observation of the second.

The density estimation itself is computed by the function \code{\link{density}} of the \code{stats} package.
}

\value{
The estimated density at 0.
}

\references{
Dehling, H., Fried, R., Wendler, M. (2015): A robust method for shift detection in time series, preprint. \href{http://arxiv.org/abs/1506.03345}{arXiv 1506.03345}
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}

\seealso{
The change point detection functions \code{\link{changerob}} and \code{\link{changerob.HL}}.
}

\examples{
set.seed(1066)
x <- rnorm(100)
y <- rnorm(20)
densdiff(x, y)
}
