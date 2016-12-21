library(robts)

set.seed(1654)

## Creating some test data:
depmodel <- list(ar=c(0.6, 0.3))
x <- arima.sim(model = depmodel, n = 100)


## Check the different methods:

arrob(x, method = "yw")
arrob(x, aic = FALSE, order.max = 20, method = "yw")

arrob(x, method = "regression")
arrob(x, aic = FALSE, order.max = 20, method = "regression")

arrob(x, method = "filter")
arrob(x, aic = FALSE, order.max = 20, method = "filter")
ARfilter(x, order.max = 5)
filterrob(x, method = "statespace")
filterrob(x, method = "recursive")

arrob(x, method = "gm")
arrob(x, aic = FALSE, order.max = 20, method = "gm")


## Check cases with missing data:

x_missing <- c(x[1:2], NA, NA, x[5:100])

# Use longest complete strech of observations:

arrob(x_missing, method = "yw", na.action = na.contiguous)

arrob(x_missing, method = "regression", na.action = na.contiguous)

arrob(x_missing, method = "filter", na.action = na.contiguous)
filterrob(x_missing, method = "statespace", na.action = na.contiguous)
filterrob(x_missing, method = "recursive", na.action = na.contiguous)

arrob(x_missing, method = "gm", na.action = na.contiguous)

# Set missings to extreme values:

arrob(x_missing, method = "yw", na.action = na.extremify)

arrob(x_missing, method = "regression", na.action = na.extremify)

arrob(x_missing, method = "filter", na.action = na.extremify)
filterrob(x_missing, method = "statespace", na.action = na.extremify)
filterrob(x_missing, method = "recursive", na.action = na.extremify)

arrob(x_missing, method = "gm", na.action = na.extremify)
