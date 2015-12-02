\encoding{utf8}
\name{qKS}
\alias{qKS}

\title{
Quantile Function of the Kolmogorov Smirnov Distribution
}
\description{
Computes the quantiles of the two sided Kolmogorov Smirnov test
}

\usage{
qKS(x)
}

\arguments{
	\item{x}{one! probability.}
}

\value{
value of the corresponding quantile
}

\details{
The distribution equals that of the supremum of the absolute values of a standard brownian bridge and therefore often occurs in change point testing. The function is only a wrapper for the C code used in the function \code{\link{ks.test}}.
}  

\author{
R Core Team and Alexander \enc{Dürre}{Duerre}
}


\seealso{
\code{\link{pKS}} for the distribution function.
}

\examples{
qKS(0.95)
}
