
## Projection onto Toeplitz matrices
## input: symmetric matrix
## output: symmetrix Toeplitz matrix
toep <- function(mat){
  #Returns projection of a symmetric matrix mat on the subspace of Toeplitz matrices
  d <- ncol(mat)
  mean_offdiag <- function(diagonal, mat) mean(mat[row(mat) == col(mat)-diagonal]) #mean of all objects on the off-diagonal 'diagonal' of the matrix 'mat' (diagonal=0 is the principal diagonal, negative values use diagonals below, positive values diagonals above the principal diagonal)
  offdiag_vector <- sapply(0:(d-1), mean_offdiag, mat=mat)
  result <- matrix(offdiag_vector[as.vector(abs(row(mat)-col(mat))+1)], ncol=d)
  return(result)
}

## Projection onto positive semidefinit matrices
## input: symmetric matrix
## output: positiv semidefinit matrix
psd <- function(mat){
  #Returns projection of a symmetric matrix mat on the subspace of positive semidefinite matrices
  specdecomp <- eigen(mat, symmetric=TRUE)
  specdecomp$values[specdecomp$values<0] <- 0
  result <- specdecomp$vectors %*% diag(specdecomp$values) %*% t(specdecomp$vectors)
  return(result)
}

is_psd <- function(mat) all(eigen(mat, symmetric=TRUE, only.values=TRUE)$values>=0) #Is a symmetric matrix mat positive semidefinite?

## Projection onto positive semidefnit Toeplitz matrices (Algorithm 2.2 of Al Homidan, 2006)
## input:
# mat: symmetric matrix
# maxit: maximal number of iterations of the projectionalgorithm
# tol: stopp iteration if difference of Frobenius-Norm is smaller tol
## output: 
# projection: positive semidefinite Toeplitz matrix
# original: input matrix
# frobenius: frobenius distance between in put and output matrix
# iteration: number of iterationsteps
psd_toep <- function(mat, maxit=100, tol=1e-8){
  #Returns projection of matrix mat on the subspace of positive semidefinite Toeplitz matrices applying Algorithm 2.2 in Al-Homidan (2006)
    #maxit: Integer. Maximal number of iterations.
    #tol: Numeric. Minimal change of the Frobenius norm between two iterations.
  F_j <- mat
  for(j in 0:maxit){
    add <- psd(toep(F_j)) - toep(F_j)
    if(sqrt(sum(add^2)) < tol) break #stop if Frobenius norm is below tol
    F_j <- F_j + add
  }
  final <- psd(toep(F_j))
  result <- list(projection=final, original=mat, frobenius=sqrt(sum((mat-final)^2)), iter=j, convergence=(j<maxit))
  return(result)
}

## Searches the nearest positive semidefinite autocorrelation function (in terms of the frobenius norm of the associated correlation matrix)
## input: autocorrelation function
## output: positive definite autocorrelation function
make_acf_psd <- function(acfvalues, ...) {
  mat <- toeplitz(c(1, acfvalues))
  pdmat <- psd_toep(mat, ...)$projection
  acfvalues_psd <- pdmat[1, -1]
  return(acfvalues_psd)
}
