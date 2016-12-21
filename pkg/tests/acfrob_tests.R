library(robts)

set.seed(1837)


## Creating some test data:
depmodel <- list(ar=c(0.6, 0.3))
depdata <- arima.sim(model = depmodel, n = 100)
datalist <- list(
  iid = rnorm(100), #iid data
  dep = depdata, #dependent data 
  short = arima.sim(model = depmodel, n = 10), #short time series
  long = arima.sim(model = depmodel, n = 1000), #long time series
  missings = c(depdata[1:2], NA, NA, depdata[5:100]) #missing values
)

maxlag <- 10
for(scenario in names(datalist)){
  print(acfrob(datalist[[scenario]], plot = FALSE, na.action = na.contiguous))
}


## Check the subroutines of all approaches:
scenario <- "dep"

acfrob.GK(datalist[[scenario]], lag.max = maxlag)
acfrob.GK(datalist[[scenario]], lag.max = maxlag, scalefn = Qn)
acfrob.GK(datalist[[scenario]], lag.max = maxlag, scalefn = Sn)
acfrob.GK(datalist[[scenario]], lag.max = maxlag, scalefn = scaleTau2)
acfrob.GK(datalist[[scenario]], lag.max = maxlag, scalefn = mad)
acfrob.GK(datalist[[scenario]], lag.max = maxlag, scalefn = sd)

acfrob.median(datalist[[scenario]], lag.max = maxlag)
acfrob.median(datalist[[scenario]], lag.max = maxlag, biascorr = TRUE)
acfrob.median(datalist[[scenario]], lag.max = maxlag, biascorr = FALSE)

acfrob.multi(datalist[[scenario]], lag.max = maxlag)
for(multi.method in c("weightedMCD", "rawMCD", "Stahel-Donoho", "S", "reweight", "Tyler", "M", "sscor")){
  print(paste("multi.method =", multi.method))
  print(acfrob.multi(datalist[[scenario]], lag.max = maxlag, multi.method = multi.method))
}

acfrob.partrank(datalist[[scenario]], lag.max = maxlag)
for(cor.method in c("spearman", "kendall", "quadrant", "gaussian", "masarotto")){
print(paste("cor.method =", cor.method))
  print(acfrob.partrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = TRUE, partial = FALSE))
  print(acfrob.partrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = FALSE, partial = FALSE))
  print(acfrob.partrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = TRUE, partial = TRUE))
  print(acfrob.partrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = FALSE, partial = TRUE))
}

acfrob.RA(datalist[[scenario]], lag.max = maxlag)
acfrob.RA(datalist[[scenario]], lag.max = maxlag, psi = "huber", locfn = median, scalefn = mad, biascorr = TRUE)
acfrob.RA(datalist[[scenario]], lag.max = maxlag, psi = "huber", locfn = median, scalefn = scaleTau2, biascorr = TRUE)
acfrob.RA(datalist[[scenario]], lag.max = maxlag, psi = "huber", locfn = median, scalefn = scaleTau2, c1 = 3, c2 = 2, biascorr = TRUE)
acfrob.RA(datalist[[scenario]], lag.max = maxlag, psi = "huber", k=1.5, locfn = median, scalefn = mad, biascorr = TRUE)
acfrob.RA(datalist[[scenario]], lag.max = maxlag, psi = "huber", k=1.5, locfn = median, scalefn = mad, biascorr = FALSE)
acfrob.RA(datalist[[scenario]], lag.max = maxlag, psi = "bisquare")

acfrob.rank(datalist[[scenario]], lag.max = maxlag)
for(cor.method in c("spearman", "kendall", "quadrant", "gaussian", "masarotto")){
print(paste("cor.method =", cor.method))
  print(acfrob.rank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = TRUE))
  print(acfrob.rank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = FALSE))
}

acfrob.filter(datalist[[scenario]], lag.max = maxlag)
for(robfil.method in c("filtered", "ar")){
print(paste("robfil.method =", robfil.method))  
  print(acfrob.filter(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = TRUE, psi.l = 2, psi.0 = 3, partial = FALSE))
  print(acfrob.filter(datalist[[scenario]], lag.max = maxlag, order.max = 5, robfil.method = "filtered", aic = TRUE, psi.l = 2, psi.0 = 3, partial = FALSE))
  print(acfrob.filter(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = FALSE, psi.l = 2, psi.0 = 3, partial = FALSE))
  print(acfrob.filter(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = TRUE, psi.l = 1.5, psi.0 = 2, partial = FALSE))
  print(acfrob.filter(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = TRUE, psi.l = 2, psi.0 = 3, partial = TRUE))
}

acfrob.trim(datalist[[scenario]], lag.max = maxlag)
acfrob.trim(datalist[[scenario]], lag.max = maxlag, trim = 0.1, biascorr = TRUE)
acfrob.trim(datalist[[scenario]], lag.max = maxlag, trim = 0.1, biascorr = FALSE)
acfrob.trim(datalist[[scenario]], lag.max = maxlag, trim = 0, biascorr = TRUE)
acfrob.trim(datalist[[scenario]], lag.max = maxlag, trim = 0.25, biascorr = TRUE)

acfrob.bireg(datalist[[scenario]], lag.max = maxlag)


## Check cases with missing data:
scenario <- "missings"

acfrob(datalist[[scenario]], plot = FALSE, na.action = na.contiguous)
acfrob(datalist[[scenario]], plot = FALSE, na.action = na.extremify)


## Check type="covariance" and type="partial":
scenario <- "dep"

acfrob(datalist[[scenario]], plot = FALSE, type = "covariance")

acfrob(datalist[[scenario]], plot = FALSE, type = "partial")
