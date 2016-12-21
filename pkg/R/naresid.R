naresid.extremify <- function(omit, x, ...) {
  x[omit] <- NA
  return(x)
}

naresid.omit <- function(omit, x, ...) {
  attr(omit, "class") <- "exclude"
  x <- naresid(omit, x, ...)
  return(x)
}
  