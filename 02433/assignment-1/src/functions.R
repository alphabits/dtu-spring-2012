source('src/utils.R')
source('src/export-helpers.R')
source('src/plot-helpers.R')


statdist = function (gamma) {
    m = dim(gamma)[1]
    one.row = rep(1, m)
    U = matrix(rep(one.row, m), nrow=m)
    Id = diag(one.row)
    inv = solve(Id - gamma + U)
    return(one.row %*% inv)
}

p.HMM = function (x, lambda, gamma, delta=FALSE) {
    if (delta == FALSE) {
        delta = statdist(gamma)
    }

    return(sum(dpois(x, lambda=lambda)*delta))
}

# From exercise 9 in week 2
pois.HMM.moments = function (m, lambda, gamma, lag.max=10) {
    delta = statdist(gamma)
    lambda = matrix(lambda)
    Lambda = diag(c(lambda))
    expectation = delta %*% lambda
    expectation.sq = expectation^2
    M = delta %*% Lambda %*% lambda
    variance = expectation + M - expectation.sq
    correlations = matrix(c(1, rep(0, lag.max)), nrow=1)
    tmp.gamma = gamma
    for (i in 2:(lag.max+1)) {
        cov = delta %*% Lambda %*% tmp.gamma %*% lambda - expectation.sq
        correlations[i] = cov/variance
        tmp.gamma = tmp.gamma %*% gamma
    }

    return(list(
        expectation=expectation,
        variance=variance,
        correlations=correlations
    ))
}

get.initial.gamma = function (m, off.diag=0.1) {
    on.diag = 1 - (m-1)*off.diag
    dat = c()
    for (column in 1:m) {
        if (column > 1) {
            dat = c(dat, rep(off.diag, column-1))
        }
        dat = c(dat, on.diag)
        if (column < m) {
            dat = c(dat, rep(off.diag, m-column))
        }
    }
    return(matrix(dat, nrow=m))
}

calculate.fit.result = function (fit.func, lambdas, off.diags, m, ...) {
    fit.res = list()
    for (lambda in lambdas) {
        for (off.diag in off.diags) {
            gamma.init = get.initial.gamma(m, off.diag)
            fit = fit.func(x, m, lambda, gamma.init, ...)
            fit$lambda.init = lambda
            fit$gamma.init = gamma.init
            fit$moments = pois.HMM.moments(m, fit$lambda, fit$gamma)
            fit.res[[(length(fit.res)+1)]] = fit
        }
    }
    return(fit.res)
}

calculate.fit.result.ml = function (lambdas, off.diags, m) {
    return(calculate.fit.result(pois.HMM.mle, lambdas, off.diags, m))
}

calculate.fit.result.em = function (lambdas, off.diags, m, ...) {
    return(calculate.fit.result(pois.HMM.EM, lambdas, off.diags, m, ...))
}
