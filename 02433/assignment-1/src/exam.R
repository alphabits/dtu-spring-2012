source('src/loaddata.R')
source('src/functions.R')
source('src/A1.R')
source('src/A2.R')


SAVEPLOTS = TRUE

off.diag.to.test = c(0.01, 0.05, 0.1)
lambdas.2 = list(c(2,8), c(3,7), c(4,6))
lambdas.3 = list(c(1,5,9), c(2,5,8), c(3,5,7))
lambdas.4 = list(c(1, 4, 6, 9), c(2, 4, 6, 8), c(1, 4, 8, 12))

time.2.ml = system.time(assign("fit.2.ml", 
        calculate.fit.result.ml(lambdas.2, off.diag.to.test, 2)))
time.2.em = system.time(assign("fit.2.em", 
        calculate.fit.result.em(lambdas.2, off.diag.to.test, 2, 
                                delta=rep(1,2)/2)))

time.3.ml = system.time(assign("fit.3.ml", 
        calculate.fit.result.ml(lambdas.3, off.diag.to.test, 3)))
time.3.em = system.time(assign("fit.3.em", 
        calculate.fit.result.em(lambdas.3, off.diag.to.test, 3, 
                                delta=rep(1,3)/3)))

time.4.ml = system.time(assign("fit.4.ml", 
        calculate.fit.result.ml(lambdas.4, off.diag.to.test, 4)))
time.4.em = system.time(assign("fit.4.em", 
        calculate.fit.result.em(lambdas.4, off.diag.to.test, 4, 
                                delta=rep(1,4)/4)))
fit.types = c("ml", "em")
for (num.state in 2:4) {
    for (fit in fit.types) {
        var.name = sprintf("time.%s.%s", num.state, fit)
        file.name = sprintf("results/timing-%s-%s.tex", num.state, fit)
        sink(file.name)
        print(get(var.name))
        sink()
    }
}

sim.pois.HMM = function (n, lambda, gamma) {
    state = 1
    x = rep(0, n)
    for (i in 1:n) {
        x[i] = rpois(1, lambda[state])
        state = which(rmultinom(1, 1, gamma[state,])==1)
    }
    return(x)
}

for (i in 1:4) {
    pdf(sprintf('plots/3-state-sim-%s.pdf', i), 12, 7)
    plot(sim.pois.HMM(242, fit.3$lambda, fit.3$gamma), type="l",
         main="Simulated data from 3-state HMM", xlab="Week", ylab="Soap sale",
         cex.axis=1.5, cex.main=1.5, cex.lab=1.5)
    dev.off()
}

ppois.HMM = function (x, delta, lambda) {
    res = 0
    for (y in 0:x) {
        for (state in 1:length(lambda)) {
            res = res + delta[state]*dpois(y, lambda[state])
        }
    }
    return(res)
}
