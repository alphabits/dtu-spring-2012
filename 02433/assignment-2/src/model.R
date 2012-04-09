inv.hessian.from.working.hessian = function (working.hessian, M) {
    inv.working.hessian = solve(working.hessian)
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

    HMM.mllk = function(working.params, ...) {
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

    HMM.mle = function(initial.params, ...) {
        initial.working.params = log(initial.params)
        res = nlm(HMM.mllk, initial.working.params, x=x, hessian=TRUE, ...)
        mle = exp(res$estimate)
        mllk = res$minimum
        num.params = length(mle)
        AIC = 2*(mllk+num.params)
        num.obs = sum(!is.na(x))
        BIC = 2*mllk+num.params*log(num.obs)
        # M is the coordinate shift matrix as defined
        # in equation 3.2 in Zucchini 2009
        M = diag(mle)
        var.cov = inv.hessian.from.working.hessian(res$hessian, M)
        list(mle=mle, code=res$code, mllk=mllk, AIC=AIC, BIC=BIC, 
             var.cov=var.cov, nlm.res=res)
    }

    HMM.log.forward.backward = function (params) { 
        R = params[R.ind]
        log.alpha = log.beta = matrix(NA, m, T)

        Ps = outer(x, p.i, dnorm, sqrt(R))
        Ps[is.na(Ps)] = 1

        gamma = calc.gamma(params)

        # Calculate alpha
        phi = Ps[1,]
        phi.sum = sum(phi)
        lscale = log(phi.sum)
        phi = phi/phi.sum
        log.alpha[,1] = log(phi)+lscale
        for (t in 2:T) {
            phi = phi%*%gamma*Ps[t,]
            phi.sum = sum(phi)
            lscale = lscale+log(phi.sum)
            phi = phi/phi.sum
            log.alpha[,t] = log(phi)+lscale
        }

        # Calculate beta
        log.beta[,T] = rep(0, m)
        phi = rep(1/m, m)
        lscale = log(m)
        for (t in (T-1):1) {
            phi = gamma%*%(Ps[t+1,]*phi)
            log.beta[,t] = log(phi)+lscale
            phi.sum = sum(phi)
            phi = phi/phi.sum
            lscale = lscale+log(phi.sum)
        }

        return(list(log.alpha=log.alpha, log.beta=log.beta))
    }

    HMM.cond.state.probs = function(params) {
        fb = HMM.log.forward.backward(params)
        la = fb$log.alpha
        lb = fb$log.beta
        c = max(la[,T])
        llk = c+log(sum(exp(la[,T]-c)))
        stateprobs = matrix(NA, ncol=T, nrow=m)
        for (i in 1:T) {
            stateprobs[,i] = exp(la[,i]+lb[,i]-llk)
        }
        return(stateprobs)
    }

    HMM.local.decoding = function(params) {
        stateprobs = HMM.cond.state.probs(params)
        ild = rep(NA, T)
        for (t in 1:T) {
            ild[t] = which.max(stateprobs[,t])
        }
        return(ild)
    }

    return(list(mle=HMM.mle, mllk=HMM.mllk, gamma=calc.gamma, 
                local.decoding=HMM.local.decoding, p.i=p.i, b.i=b.i,
                data=x))
}

discrete.dnorm = function (x, mu, sd, resolution) {
    return(pnorm(x+0.5*resolution, mu, sd) - pnorm(x-0.5*resolution, mu, sd))
}
