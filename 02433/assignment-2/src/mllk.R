inv.hessian.from.working.hessian = function (working.hessian, minimum) {
    inv.working.hessian = solve(working.hessian)
    M = diag(minimum)
    return(t(M) %*% inv.working.hessian %*% M)
}

get.HMM.model = function (x, m, state.bounds=c(2.1, 8.4)) {
    T = length(x)
    w = diff(state.bounds)/m
    i = 1:m
    p.i = state.bounds[1] + w*(i - 0.5)
    b.i = state.bounds[1] + w*i
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

    HMM.mllk = function(working.params) {
        params = exp(working.params)
        R = params[R.ind]

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

    HMM.mle = function(initial.params) {
        initial.working.params = log(initial.params)
        res = nlm(HMM.mllk, initial.working.params, x=x, hessian=TRUE)
        mle = exp(res$estimate)
        mllk = res$minimum
        num.params = length(params)
        AIC = 2*(mllk+num.params)
        num.obs = sum(!is.na(x))
        BIC = 2*mllk+num.params*log(num.obs)
        var.cov = inv.hessian.from.working.hessian(res$hessian, params)
        list(mle=mle, code=res$code, mllk=mllk, AIC=AIC, BIC=BIC,
             var.cov=var.cov)
    }

    return(list(mle=HMM.mle, mllk=HMM.mllk))
}




