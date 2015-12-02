\encoding{utf8}
\name{ARfilter}
\alias{ARfilter}

\title{
Estimation of AR Models by Robust Filtering
}
\description{
Robustly estimates AR models using a robust filter algorithm. See \enc{Dürre}{Duerre} et al. (2014) for details.
}

\usage{
ARfilter(timeseries,p,aicpenalty=function(p) {return(2*p)},psi.l=2,psi.0=3)
}

\arguments{
	\item{timeseries}{univariate numeric vector or time series object.}
	\item{p}{numeric determining the order of the AR fit.}
	\item{psi.l}{numeric determining the used psi function, see details.}
	\item{psi.0}{numeric determining the used psi function, see details.}
	\item{aicpenalty}{function of p, indicating the penalty for a larger model, see details.}
	}

\value{
List containing the following values

\item{partial autocorrelations}{matrix representing the successive estimated partial correlations of the model. In the first column are the partial correlations of the AR model of order 1 and so on.}
\item{variance of innovations}{numeric containing the successive scale estimation of the residuals of the successive AR models. The first entry is for the AR(1) model and so on.}
\item{autocorrelation}{numeric of estimated acf up until lag p.}
\item{aic}{numeric of aic values of successive AR models.}
\item{robustly filtered timeseries}{matrix containing the robustly filtered time series of successive AR models. The first column corresponds to the AR(1) model and so on.}
\item{AR-coefficients}{matrix containing the estimated AR coefficients of successive AR models. The first column corresponds to the AR(1) model and the first row to the first AR coefficient and so on.}
}

\details{

The function fits AR models of increasing order by a robust version of the Durbin Levinson algorithm as described in chapter 8 of Maronna et al. (2006). The AR coefficients are estimated by minimizing a robust scale (\code{\link{scaleTau2}}) of one step ahead residuals of robustly filtered predictors. The filter process basically compares the predicted value under the estimated model with the observed one and transformes the corresponding residual with the psi function. The transformed value is then set to the sum of the predicted value and the transformed residual. For more details see Maronna et al. (2006) and \enc{Dürre}{Duerre} et al. (2014).

Following Maronna et al. (2006), the psi function should fulfill two properties. It should be the identity for small absolute values and 0 for large absolute values. Both thresholds can be determined by \code{psi.l} and \code{psi.0}. Here the psi function is choosen to be two times continuous differentiable, see \enc{Dürre}{Duerre} et al. (2014) for a formal definition.

}

\references{

\enc{Dürre}{Duerre}, A., Fried, R. and Liboschik, T. (2014): Robust estimation of (partial) autocorrelation, working paper, \url{http://hdl.handle.net/2003/33011}.

Maronna, R., Martin, D. and Yohai, V. (2006): \emph{Robust statistics}, Wiley, Chichester.

}

\author{
Alexander \enc{Dürre}{Duerre}, Tobias Liboschik and Jonathan Rathjens
}


\seealso{

The wrapper function \code{\link{arrob}}.

}

\examples{
set.seed(1066)
tss <- arima.sim(model = list(ar = 0.3, ma = 5), n = 100)
ARfilter(tss,p=5)
}
