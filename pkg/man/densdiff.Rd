\encoding{utf8}
\name{densdiff}
\alias{densdiff}

\title{
Density Estimation at 0 of the Random Variable X-Y
}
\description{
Estimates the density of the difference of two samples at 0. This is necessary for the Hodges-Lehmann (changepoint) test. See Dehling et al. (2015) for details.
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
	\item{adjust}{numeric which determines the used bandwith, see the help of \code{\link{density}} for Details.}
	\item{kernelused}{character string indicating the used kernel, see the help of \code{\link{density}} for Details.}
	\item{...}{further arguments passed to the respective internal functions.}
}

\value{
The estimated density at 0.
}

\details{
The Hodges-Lehmann two sample test requires an estimation of the density of the difference of both random variables X and Y at 0. This is also true for the the Hodges-Lehmann changepoint test since it consists of a sequence of two samples tests. One approach to estimate this density is to apply a kernel estimator at differences of values of both samples. This is basically the task of this function.

There are different possibilities to build differences.

If 'type2' is 'all' one uses all possible differences of the merged sample.

If 'type2' is 'within' one uses all possible differences within the first sample and all possible differences within the second.

If 'type2' is 'between' one uses all possible differences between one observation of the first sample and one observation of the second.

The density estimation itself is computed by the function \code{\link{density}} of the stats package.
}

\references{

Dehling, H., Fried, R., Wendler, M. (2015): A robust method for shift detection in time series, preprint. 
\href{http://arxiv.org/abs/1506.03345}{arvix 1506.03345}
}

\author{
Roland Fried and Alexander \enc{Dürre}{Duerre}
}


\seealso{

The change point detection functions \code{\link{rob.change}} and \code{\link{strucchange.HL}}

}

\examples{
set.seed(1066)
x <- rnorm(100)
y <- rnorm(20)
densdiff(x,y)
}
