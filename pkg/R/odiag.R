####################
# auxiliary function
# extract offdiagonal elements (from a former R-package demogR)
# input
# A: matrix
# at: number of off-diagonal
# output: offdiagomnal as vector
####################

odiag <- function (A, at = 0) 
{
    if (is.matrix(A)) {
        y <- A[row(A) == col(A) - at]
        return(y)
    }
    else {
        len <- length(A)
        B <- matrix(0, nrow = len + abs(at), ncol = len + abs(at))
        B[row(B) == col(B) - at] <- A
        return(B)
    }
}