napredict.extremify <- function(omit, x, ...) {
  return(x)
}

napredict.omit <- function(omit, x, ...) {
  attr(omit, "class") <- "exclude"
  x <- naresid(omit, x, ...)
  return(x)
}
 