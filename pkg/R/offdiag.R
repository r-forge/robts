####################
# auxiliary function
# extracts offdiagonal elements of symmetric matrix or builds matrix with a given offdiagonal
# input:
# 		x: matrix or vector
# 		at: number of off-diagonal
# output: offdiagonal as vector
####################

offdiag <- function (x, at = 0) {
  if (is.matrix(x)) {
    result <- x[row(x) == col(x) - at]
  } else {
    len <- length(x)
    result <- matrix(0, nrow = len + abs(at), ncol = len + abs(at))
    result[row(result) == col(result) - at] <- x
  }
  return(result)
}
