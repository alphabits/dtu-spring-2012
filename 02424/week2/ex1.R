dat = c(71, 74, 82, 76, 91, 82, 82, 75, 79, 82, 72, 90)

gaussian.likelihood = function (mean, sd, y) {
    return(prod(dnorm(y, mean=mean, sd=sd)))
}

plot.likelihood = function (y) {
    dat.mean = mean(y)
    dat.sd = sd(y)
    x = seq(dat.mean-2*dat.sd, dat.mean+2*dat.sd, length=200)
    likelihood = sapply(x, gaussian.likelihood, dat.sd, y)
    plot(x, likelihood/max(likelihood))
}
