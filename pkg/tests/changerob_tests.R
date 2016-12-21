library(robts)

set.seed(1654)

## Creating some test data:
depmodel <- list(ar=c(0.6, 0.3))
x <- arima.sim(model = depmodel, n = 100)


## Check the different methods:
test <- "HL"
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = TRUE, plot = TRUE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = TRUE, shiftcorrect = TRUE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = FALSE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = TRUE)

test <- "Wilcoxon"
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = TRUE, plot = TRUE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = TRUE, shiftcorrect = TRUE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = FALSE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = TRUE)

test <- "CUSUM"
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = TRUE, plot = TRUE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = TRUE, shiftcorrect = TRUE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = FALSE)
changerob(x, property = "location", test = test, alternative = "two.sided", var.method = "window", overlapping = FALSE, shiftcorrect = TRUE)
