get.mllk.and.mle = function (m) {
    m = 250
    w = 6.3/m
    i = 1:m
    p.i = 2.1 + w*(i - 0.5)
    b.i = 2.1 + w*i
    theta.ind = 1
    r0.ind = 2
    K.ind = 3
    Q.ind = 4
    R.ind = 5

    calc.gamma = function (params) {
        theta = params[theta.ind]
        r0 = params[r0.ind]
        K = params[K.ind]
        Q = params[Q.ind]

        mu = p.i + r0*(1 - (exp(p.i)/K)^theta)
        trans.prob = matrix(0, m, m)

        for (i in 1:m) {
            row = (dnorm(b.i-w, mu[i], sqrt(Q)) + dnorm(b.i, mu[i], sqrt(Q)))*w/2
            trans.prob[i,] = row/sum(row)
        }

        return(trans.prob)
    }

    HMM.mllk = function(working.params, x, ...) {
        params = exp(working.params)
        R = params[R.ind]
        T = length(x)

        Ps = outer(x, p.i, dnorm, sqrt(R))
        Ps[is.na(Ps)] = 1
        gamma = calc.gamma(params)

        phi = Ps[1,]
        phi.sum = sum(phi)
        llk = log(phi.sum)
        phi = phi/phi.sum

        for (t in 2:T) {
            phi = phi %*% gamma * Ps[t,]
            phi.sum = sum(phi)
            llk = llk + log(phi.sum)
            phi = phi/phi.sum
        }

        return(-llk)
    }

    HMM.mle = function(x, params, ...) {
        working.params = log(params)
        res = nlm(HMM.mllk, working.params, x=x, hessian=TRUE)
        params = exp(res$estimate)
        mllk = res$minimum
        num.params = length(params)
        AIC = 2*(mllk+num.params)
        num.obs = sum(!is.na(x))
        BIC = 2*mllk+num.params*log(num.obs)
        list(params=params, code=res$code, mllk=mllk, AIC=AIC, BIC=BIC,
             hessian=res$hessian)
    }

    return(list(mllk=HMM.mllk, mle=HMM.mle))
}

inv.hessian.from.working.hessian = function (working.hessian, minimum) {
    inv.working.hessian = solve(working.hessian)
    M = diag(minimum)
    return(t(M) %*% inv.working.hessian %*% M)
}


# After some trial and error these params might be good starting params
params = c(0.5, 0.2, 1000, 0.01, 0.05)


