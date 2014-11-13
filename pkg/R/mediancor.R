####################
# auxiliary function: calculate two-dimensional correlation using median
# input
# x: first Variable (demeaned) to calculate correlation as vector
# y: second Variable (demeaned) to calculate correlation as vector
# output: something like a correlation (biased!)
####################

 
mediancor <- function(x,y) {
	return(median(x*y)/median(x^2))
}