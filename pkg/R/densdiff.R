## estimates density around of differences at 0
## input:
# x: first sample
# y: second sample
# type2: what differences should be considered?
# 	possible are: 	between (all possible pairs where minuend is from first sample and subtrahend from the second) 
#			all (all possible pairs except for the same values)
#			within (merge all possible pairs in first sample and all possible pairs in second sample)
# adjust: determine the bandwith used by density (see density)
## output:
# density estimation at 0
densdiff <- function(x, y, type2 = c("all", "within", "between"), adjust = 1, kernelused = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"), ...){
 kernelused <- match.arg(kernelused)
 type2 <- match.arg(type2)
 if (type2 == "between"){  # only pairs of x and y
   dif <- as.numeric(outer(x, y, "-"))}
 if (type2 == "all"){ # pairs of all
   z <- c(x, y)
   dif <- as.numeric(dist(z))}
 if (type2 == "within"){  # only pairs within x or within y
   dif1 <- as.numeric(dist(x))
   dif2 <- as.numeric(dist(y))
   dif <- c(dif1, dif2)}
 out <- density(dif, na.rm = TRUE, kernel = kernelused, adjust = adjust, from = 0, to = 0, n = 1)
  return(out$y)
}
