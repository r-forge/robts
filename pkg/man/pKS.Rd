\encoding{utf8}
\name{pKS}
\alias{pKS}

\title{
Distribution Function of the Kolmogorov Smirnov Distribution
}
\description{
Computes the asymptotic distribution of the two sided Kolmogorov Smirnov test
}

\usage{
pKS(x)
}

\arguments{
	\item{x}{vector of quantiles.}
}

\value{
Values of the distribution function of the given quantiles.
}

\details{
The distribution equals that of the supremum of the absolute values of a standard brownian bridge and often occurs in change point testing. The function is only a wrapper for the C code used in the function \code{\link{ks.test}} .
}  

\author{
R Core Team and Alexander \enc{Dürre}{Duerre}
}


\seealso{
\code{\link{ks.test}}
}

\examples{
curve(pKS,from=0,to=2)
}
