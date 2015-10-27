### robust changepoint tests ###
## robust tests for a change in location or scale
# input:
# x: timeseries
# property: whether to test against a change in location or scale
# procedure: which test to use: usual Cusum test, Wilcoxon-test and Hodges-Lehmann test are possible
# conf.level: desired confidence level of the test
# alternative: power against which alternative
# 	possibilities:	two.sided: classical approach
#			greater: mean increases
#			less: mean decreases
# varmethod: procedure to estimate the long run variance
#	possibilities:	window: subsampling procedure
#			acf: kernel estimator
#			acfextra: extrapolates acf from ar-fit
# overlapping: should distinct or overlapping blocks be used to calculate the lonmg run variance, only used if varmethod=window
# shiftcorrect: If TRUE the estimation of the lon run variance considers a location change
# plot: if true trajectory of test is drawn

rob.change <- function(x,property=c("location","scale"),procedure=c("HL","Wilcoxon","Cusum"),conf.level=0.95,alternative=c("two.sided","less","greater"),varmethod=c("window","acf","acfextra"),overlapping=FALSE,shiftcorrect=TRUE,borderN=10,plot=FALSE,...) {
n <- length(x)
if (shiftcorrect) {
	if (borderN > n/4) warning(paste("Setting borderN to ",floor(n/4)," because it was to large."))
	borderN <- min(floor(n/4),borderN)
}

procedure <- match.arg(procedure)
property <- match.arg(property)
alternative <- match.arg(alternative)
varmethod <- match.arg(varmethod)

if (property=="scale") {
	x <- abs(x-median(x))
	x0 <- x==0
	xmin <- min(x[!x0])
	x[x0] <- xmin 
	}

if (procedure=="Cusum") {
	testtrajectory <- strucchange.cusum(x=x,varmethod=varmethod,overlapping=overlapping,shiftcorrect=shiftcorrect,borderN=borderN,...)
}
if (procedure=="Wilcoxon") {
	testtrajectory <- strucchange.wilcox(x=x,varmethod=varmethod,overlapping=overlapping,shiftcorrect=shiftcorrect,borderN=borderN,...)
}
if (procedure=="HL") {
	testtrajectory <- strucchange.HL(x=x,varmethod=varmethod,overlapping=overlapping,shiftcorrect=shiftcorrect,borderN=borderN,...)
}

data.name <- names(x)
if (alternative=="two.sided") {
	p.value <- 1-.Call(stats:::C_pKS2, p = max(abs(testtrajectory)), 10^(-6))
	statistic <- max(abs(testtrajectory))
	estimate <- which.max(abs(testtrajectory))
		if (plot==TRUE) {
			border <- quantileKS(1-conf.level)
			mini <- min(c(-border,testtrajectory))
			maxi <- max(c(border,testtrajectory))
			plot(testtrajectory,type="l",main=data.name,ylim=c(mini,maxi))
			abline(h=border,lty="dotted")
			abline(h=-border,lty="dotted")
			abline(v=estimate,lty="dashed")
			text(x=estimate-0.12*n,y=mini+0.12*(maxi-mini),paste("cp: ",estimate,"
p-value: ",round(p.value,3)))
		}
	}

if (alternative=="greater") {
	p.value <- exp(-max(c(0,testtrajectory))^2)
	statistic <- max(testtrajectory)
	estimate <- which.max(testtrajectory)
		if (plot==TRUE) {
			border <- sqrt(-log(1-conf.level))
			mini <- min(testtrajectory)
			maxi <- max(c(testtrajectory,border))
			plot(testtrajectory,type="l",main=data.name,ylim=c(mini,maxi))
			abline(h=border,lty="dotted")
			abline(v=estimate,lty="dashed")
			text(x=estimate-0.12*n,y=mini+0.12*(maxi-mini),paste("cp: ",estimate,"
p-value: ",round(p.value,3)))
		}
	}

if (alternative=="less") {
	p.value <- exp(-min(c(0,testtrajectory))^2)
	statistic <- min(testtrajectory)
	estimate <- which.min(testtrajectory)
		if (plot==TRUE) {
			border <- -sqrt(-log(1-conf.level))
			mini <- min(c(border,testtrajectory))
			maxi <- max(testtrajectory)
			plot(testtrajectory,type="l",main=data.name,ylim=c(mini,maxi))
			abline(h=border,lty="dotted")
			abline(v=estimate,lty="dashed")
			text(x=estimate-0.12*n,y=mini+0.12*(maxi-mini),paste("cp: ",estimate,"
p-value: ",round(p.value,3)))
		}
	}
names(estimate) <- "estimated change point"
names(statistic) <- "test statistic"
null.value <- 0
names(null.value) <- "possible level shift"
erg <- list(statistic=statistic,p.value=p.value,estimate=estimate,null.value=null.value,
	alternative=alternative,method=paste(procedure," for change in ",property),data.name=data.name,trajectory=testtrajectory)
class(erg) <- "htest"
return(erg)
}

