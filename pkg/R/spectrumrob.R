spectrumrob <- function(x, method = c("pgram", "acf"), plot = TRUE, ...){
  method <- match.arg(method)
  res <- switch(method, pgram = spectrumrob.pgram(x, ...), acf = spectrumrob.acf(x, ...))
	if (plot) {
		plot(res)
		return(invisible(res))
	} else {
    return(res)
  }
}
