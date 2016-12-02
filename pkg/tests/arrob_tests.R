set.seed(1654)

library(robts)
x <- arima.sim(n=200, model=list(ar=c(0.6,0.3)))

ARfilter(x, order.max=3, aicpenalty = function(p) 2*p, psi.l = 2, psi.0 = 3)


output <- arrob.gm(x, order.max = 10, aic = TRUE, aicpenalty=function(p) 2*p)
output
output[]

output <- arrob.regression(x, order.max = 1, aic = TRUE, aicpenalty=function(p) 2*p)
output
output[]

output <- arrob.filter(x, order.max = 3, aic = FALSE, aicpenalty=function(p) 2*p, psi.l = 2, psi.0 = 3)
output
output[]

output <- arrob.yw(x, order.max = 1, aic = TRUE, aicpenalty=function(p) 2*p, acf.approach="GK")
output
output[]

stats:::ar.burg(x, order.max=30, aic=T)[]
