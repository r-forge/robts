pKS <- function(q) c(.Call("pKS2", q, tol = 10^(-6)))

qKS <- function(p) {
  pKSm <- function(y) return(pKS(y)-p)
  quan <- uniroot(pKSm, lower = 0.2, upper = 3)
  res <- quan$root
  return(res)
}
