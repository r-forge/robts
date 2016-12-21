handle_missings_ts <- function(x, na.action) {
  if(anyNA(x)) {
    class_original <- class(x)
    class(x) <- "ts"    
    x <- na.action(x)
  	na.attr <- attr(x, "na.action")
  	if(!is.null(na.attr)) {
      na.proportion <- length(na.attr)/length(x)
      if(attr(na.attr, "class") == "extremify") {
        if(na.proportion >= 0.25) stop("The data contain more than 25% missing values. Handling of missing values with\nthe method 'na.extremify' does not work for proportions higher than the\nbreakdown point of the estimation procedure.")
        warning(paste("The proportion of NAs which have been replaced by an extreme value is about ",  round(na.proportion*100), "%\nof the data.", sep=""))
      }
    class(x) <- class_original
  	}
  	if(anyNA(x)) stop("Please use a different method for handling missing values (specified by the\nargument 'na.action') which does not leave missing values in the data.")
	}
	return(x)
}
