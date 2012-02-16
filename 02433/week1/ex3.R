# Write a function to minimize minus log-likelihood
# of a Poisson mixture model with m components, 
# using the nonlinear minimizer nlm

x = c(13, 14, 8, 10, 16, 26, 32, 27, 18, 32, 36, 24, 22, 23, 22, 18, 25, 21, 21, 14, 8, 11, 14, 23, 18, 17, 19, 20, 22, 19, 13, 26, 13, 14, 22, 24, 21, 22, 26, 21, 23, 24, 27, 41, 31, 27, 35, 26, 28, 36, 39, 21, 17, 22, 17, 19, 15, 34, 10, 15, 22, 18, 15, 20, 15, 22, 19, 16, 30, 27, 29, 23, 20, 16, 21, 21, 25, 16, 18, 15, 18, 14, 10, 15, 8, 15, 6, 11, 8, 7, 18, 16, 13, 12, 13, 20, 15, 16, 12, 18, 15, 16, 13, 15, 16, 11, 11)

transform.params = function (lambda, delta) {
    eta = log(lambda)
    deltasum = sum(delta[2:length(delta)])
    tau = log(delta/(1 - deltasum))
    return(list(eta=eta, tau=tau[2:length(tau)]))
}

recover.params = function (eta, tau) {
    lambda = exp(eta)
    tausum = sum(exp(tau))
    delta = exp(tau)/(1+tausum)
    delta.one = 1 - sum(delta)
    return(list(delta=c(delta.one, delta), lambda=lambda))
}

get.eta.tau.from.param = function (p) {
    m = length(p)
    divide.at = ceiling(m/2)
    eta = p[1:divide.at]
    tau = p[(divide.at+1):m]
    return(list(eta=eta, tau=tau))
}

minus.log.likelihood = function (p, x) {
    eta.tau = get.eta.tau.from.param(p)
    params = recover.params(eta.tau[["eta"]], eta.tau[["tau"]])
    lambda = params[["lambda"]]
    delta = params[["delta"]]
    return(-sum(log(outer(x, lambda, dpois)%*%delta)))
}

fit.data = function (m, x) {
    #delta.init = rep(1/(m+3), m)
    #lambda.init = rep(3, m)
    #trans = transform.params(lambda.init, delta.init)
    #params.init = c(trans[["eta"]], trans[["tau"]])
    params.init = rep(10, 2*m-1)
    return(fit.data.init.params(params.init, x))
}

fit.data.init.params = function (init.params, x) {
    res = nlm(minus.log.likelihood, init.params, x=x)
    eta.tau = get.eta.tau.from.param(res[["estimate"]])
    params = recover.params(eta.tau[["eta"]], eta.tau[["tau"]])
    return(params)
}
