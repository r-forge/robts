#####################
# auxiliary function: computes autocorrelations from partial correlations
# input: partial autocorrelations as vector
# output: autocorrelations
#####################

solveYuleII <- function(partial) {
p <- length(partial)

# solving the yule-walker equations
# building a linear equation system 

A <- diag(rep(-1,p))
for (i in 1:(p)) {
	for (j in 1:p) {
		if(i!=j){
			A[i,abs(i-j)] <- A[i,abs(i-j)]+partial[j]
			}
		}
	}

# solving the system

erg <- try(solve(A,-partial),silent=TRUE)
if (inherits(erg,"try-error")){
	warning("Solving the Yule-Walker equations failed.")
	return(NA)
	}

return(erg)
}
