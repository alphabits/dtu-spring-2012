# The number of customers arriving at a cafe per 10 minutes
# are given by

dat = c(4, 6, 3, 7, 2, 4)

# We assume that the data are an iid sample from the Poisson 
# distribution. 
# We want to plot the log likelihood function so we write a 
# function to calculate the log likelihood from a dataset

log.likelihood = function (lambda, x) {
    return(sum(dpois(x, lambda, log=TRUE)))
}

# To actually plot the log likelihood we write another function

plot.log.likelihood = function (x) {
    m = mean(x)
    s = m
    rng = seq(max(m-3*s, 0), m+3*s, length=200)
    vals = sapply(rng, log.likelihood, x=x)
    plot(rng, vals)
}

# We also need to examine the likelihood so we write a general 
# function

likelihood = function (lambda, x) {
    return(prod(dpois(x, lambda)))
}

# And then realize that we need the plotting function from
# above with only minor changes


plot.likelihood.function = function (x, likelihood.fns) {
    m = mean(x)
    s = m
    rng = seq(max(m-3*s, 0), m+3*s, length=200)
    vals = sapply(rng, likelihood.fns, x=x)
    plot(rng, vals)
}

# So now we can write the previous plot.log.likelihood as

plot.log.likelihood = function (x) {
    plot.likelihood.function(x, log.likelihood)
}

plot.likelihood = function (x) {
    plot.likelihood.function(x, likelihood)
}
