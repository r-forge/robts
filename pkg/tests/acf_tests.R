## Creating some test data:
depmodel <- list(ar=0.5, ma=0.3)
depdata <- arima.sim(model = depmodel, n = 100)
datalist <- list(
  iid = rnorm(100), #iid data
  dep = depdata, #dependent data 
  short = arima.sim(model = depmodel, n = 10), #short time series
  long = arima.sim(model = depmodel, n = 1000) #long time series
)

library(robts)

maxlag <- 10
for(scenario in names(datalist)){
  print(acfrob(datalist[[scenario]], plot = FALSE))
}


## Check the subroutines of all approaches:
scenario <- "dep"

acfGK(datalist[[scenario]], lag.max = maxlag)
acfGK(datalist[[scenario]], lag.max = maxlag, scalefn = Qn)
acfGK(datalist[[scenario]], lag.max = maxlag, scalefn = Sn)
acfGK(datalist[[scenario]], lag.max = maxlag, scalefn = scaleTau2)
acfGK(datalist[[scenario]], lag.max = maxlag, scalefn = mad)
acfGK(datalist[[scenario]], lag.max = maxlag, scalefn = sd)

acfmedian(datalist[[scenario]], lag.max = maxlag)
acfmedian(datalist[[scenario]], lag.max = maxlag, biascorr = TRUE)
acfmedian(datalist[[scenario]], lag.max = maxlag, biascorr = FALSE)

acfmulti(datalist[[scenario]], lag.max = maxlag)
for(multi.method in c("weightedMCD", "rawMCD", "Stahel-Donoho", "S", "reweight", "Tyler", "M", "sscor")){
  print(paste("multi.method =", multi.method))
  print(acfmulti(datalist[[scenario]], lag.max = maxlag, multi.method = multi.method))
}

acfpartrank(datalist[[scenario]], lag.max = maxlag)
for(cor.method in c("spearman", "kendall", "quadrant", "gaussian", "masarotto")){
print(paste("cor.method =", cor.method))
  print(acfpartrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = TRUE, partial = FALSE))
  print(acfpartrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = FALSE, partial = FALSE))
  print(acfpartrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = TRUE, partial = TRUE))
  print(acfpartrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = FALSE, partial = TRUE))
}

acfRA(datalist[[scenario]], lag.max = maxlag)
acfRA(datalist[[scenario]], lag.max = maxlag, psi = "huber", locfn = median, scalefn = mad, biascorr = TRUE)
acfRA(datalist[[scenario]], lag.max = maxlag, psi = "huber", locfn = median, scalefn = scaleTau2, biascorr = TRUE)
acfRA(datalist[[scenario]], lag.max = maxlag, psi = "huber", locfn = median, scalefn = scaleTau2, c1 = 3, c2 = 2, biascorr = TRUE)
acfRA(datalist[[scenario]], lag.max = maxlag, psi = "huber", k=1.5, locfn = median, scalefn = mad, biascorr = TRUE)
acfRA(datalist[[scenario]], lag.max = maxlag, psi = "huber", k=1.5, locfn = median, scalefn = mad, biascorr = FALSE)
acfRA(datalist[[scenario]], lag.max = maxlag, psi = "bisquare")

acfrank(datalist[[scenario]], lag.max = maxlag)
for(cor.method in c("spearman", "kendall", "quadrant", "gaussian", "masarotto")){
print(paste("cor.method =", cor.method))
  print(acfrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = TRUE))
  print(acfrank(datalist[[scenario]], lag.max = maxlag, cor.method = cor.method, biascorr = FALSE))
}

acfrobfil(datalist[[scenario]], lag.max = maxlag)
for(robfil.method in c("filtered", "ar")){
print(paste("robfil.method =", robfil.method))  
  print(acfrobfil(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = TRUE, psi.l = 2, psi.0 = 3, partial = FALSE))
  print(acfrobfil(datalist[[scenario]], lag.max = maxlag, p = 5, robfil.method = "filtered", aic = TRUE, psi.l = 2, psi.0 = 3, partial = FALSE))
  print(acfrobfil(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = FALSE, psi.l = 2, psi.0 = 3, partial = FALSE))
  print(acfrobfil(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = TRUE, psi.l = 1.5, psi.0 = 2, partial = FALSE))
  print(acfrobfil(datalist[[scenario]], lag.max = maxlag, robfil.method = "filtered", aic = TRUE, psi.l = 2, psi.0 = 3, partial = TRUE))
}

acftrim(datalist[[scenario]], lag.max = maxlag)
acftrim(datalist[[scenario]], lag.max = maxlag, trim = 0.1, biascorr = TRUE)
acftrim(datalist[[scenario]], lag.max = maxlag, trim = 0.1, biascorr = FALSE)
acftrim(datalist[[scenario]], lag.max = maxlag, trim = 0, biascorr = TRUE)
acftrim(datalist[[scenario]], lag.max = maxlag, trim = 0.25, biascorr = TRUE)

