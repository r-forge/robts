na.extremify <- function(object, ...) UseMethod("na.extremify")

na.extremify.ts <- function(object, ...) {
  replacement <- seq(along=object)[is.na(object)]
  if (any(replacement > 0L)) {
    roundUp <- function(x) 10^ceiling(log10(x))
    amplifier <- 1e+02
    extreme <- roundUp(max(object, na.rm=TRUE) + amplifier*diff(range(object, na.rm=TRUE)))
    object[replacement] <- extreme
    proportionNA <- length(replacement)/length(object)
    attr(replacement, "class") <- "extremify"
    attr(object, "na.action") <- replacement
  }
  return(object)
}
