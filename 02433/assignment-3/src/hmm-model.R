inv.hessian.from.working.hessian = function (working.hessian, M) {
    inv.working.hessian = solve(working.hessian)
    return(t(M) %*% inv.working.hessian %*% M)
}

statdist = function (gamma) {
    m = dim(gamma)[1]
    one.row = rep(1, m)
    U = matrix(rep(one.row, m), nrow=m)
    Id = diag(one.row)
    inv = solve(Id - gamma + U)
    return(one.row %*% inv)
}

get.M.from.natural.params = function (natural.params) {
    mu = natural.params$mu
    sd = natural.params$sd
    gamma = as.vector(natural.params$gamma)
    m = length(mu)
    gamma.submatrix = matrix(NA, m*(m-1), m^2)
    rows = m*(m+1)
    cols = m*(m+2)
    M = matrix(0, rows, cols)
    M[1:m,1:m] = diag(mu)
    M[(m+1):(2*m),(m+1):(2*m)] = diag(sd)

    for (row.idx in 1:(m*(m-1))) {
        for (col.idx in 1:m^2) {
            fake.row.idx = row.idx + ceiling(row.idx/m)
            if (fake.row.idx == col.idx) {
                value = gamma[col.idx]*(1-gamma[col.idx])
            } else if (abs(fake.row.idx-col.idx) == m) {
                value = -gamma[fake.row.idx]*gamma[col.idx]
            } else {
                value = 0
            }
            gamma.submatrix[row.idx, col.idx] = value
        }
    }

    M[(2*m+1):(m*(m+1)),(2*m+1):(m*(m+2))] = gamma.submatrix

    return(M)
}

my.log = function (x) {
    return(log(max(x, 1e-16)))
}

get.HMM.model = function (x, m) {
    T = length(x)

    HMM.gamma.natural.to.working = function (gamma) {
        working.matrix = log(gamma/diag(gamma))
        # Return off diagonal elements
        return(working.matrix[!diag(m)])
    }

    HMM.gamma.working.to.natural = function (working.gamma) {
        natural.gamma = diag(m)
        # Fill the off diagonal elements
        natural.gamma[!natural.gamma] = exp(working.gamma)
        natural.gamma = natural.gamma/rowSums(natural.gamma)
        return(natural.gamma)
    }

    HMM.natural.params.to.working = function (mu, sd, gamma) {
        return(c(log(mu), log(sd), HMM.gamma.natural.to.working(gamma)))
    }

    HMM.working.params.to.natural = function (parvector, return.list=TRUE) {
        mu = exp(parvector[1:m])
        sd = exp(parvector[(m+1):(2*m)])
        gamma = HMM.gamma.working.to.natural(parvector[(2*m+1):(m*(m+1))])
        if (return.list) {
            return(list(mu=mu, sd=sd, gamma=gamma))
        } else {
            return(c(mu, sd, as.vector(gamma)))
        }
    }

    HMM.get.probs.matrix = function (mu, sd, dat=NULL) {
        if (!is.null(dat)) {
            x = dat
            T = length(x)
        }

        probs = matrix(0, T, m)
        for (state in 1:m) {
            probs[,state] = discrete.dnorm(x, mu[state], sd[state], 1)
        }
        return(probs)
    }

    HMM.mllk = function(working.params, debug.info=FALSE, ...) {
        natural.params = HMM.working.params.to.natural(working.params)
        mu = natural.params$mu
        sd = natural.params$sd
        gamma = natural.params$gamma

        # For debugging
        print(mu)
        print(sd)
        print(gamma)

        Ps = HMM.get.probs.matrix(mu, sd)
        Ps[is.na(Ps)] = 1

        phi = Ps[1,]
        phi.sum = sum(phi)
        llk = my.log(phi.sum)
        phi = phi/max(phi.sum, 1e-16)

        for (t in 2:T) {
            phi = phi %*% gamma * Ps[t,]
            phi.sum = sum(phi)
            llk = llk + my.log(phi.sum)
            phi = phi/max(phi.sum, 1e-16)
        }

        return(-llk)
    }

    HMM.mle = function(mu0, sd0, gamma0, ...) {
        initial.working.params = HMM.natural.params.to.working(mu0, sd0, gamma0)
        res = nlm(HMM.mllk, initial.working.params, x=x, hessian=TRUE, ...)
        mle = HMM.working.params.to.natural(res$estimate, return.list=FALSE)
        mle.list = HMM.working.params.to.natural(res$estimate)
        mllk = res$minimum
        num.params = length(mle)
        AIC = 2*(mllk+num.params)
        num.obs = sum(!is.na(x))
        BIC = 2*mllk+num.params*log(num.obs)
        # M is the coordinate shift matrix as defined
        # in equation 3.2 in Zucchini 2009
        list(mle=mle, code=res$code, mllk=mllk, AIC=AIC, BIC=BIC, 
             nlm.res=res, mle.list=mle.list)
    }

    HMM.log.forward.backward = function (params, dat=NULL) { 
        if (!is.null(dat)) {
            x = dat
            T = length(x)
        }

        log.alpha = log.beta = matrix(NA, m, T)

        mu = params$mu
        sd = params$sd
        gamma = params$gamma

        Ps = HMM.get.probs.matrix(mu, sd, x)
        Ps[is.na(Ps)] = 1

        # Calculate alpha
        delta = statdist(gamma)
        phi = delta*Ps[1,]
        phi.sum = sum(phi)
        lscale = my.log(phi.sum)
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

    HMM.one.step.predictions = function (data.start.t, predict.x, params) {
        alpha.Ts = HMM.log.forward.backward(params, c(power.production[data.start.t:8000], predict.x))$log.alpha
        dist.resolution = 200
        xs = seq(0, 21000, length.out=dist.resolution)
        mu = params$mu
        sd = params$sd
        gamma = params$gamma
        num.preds = length(predict.x)
        pred.distributions = matrix(0, dist.resolution, num.preds)

        for (pred.num in 1:num.preds) {

            alpha.T = t(alpha.Ts[,(8000+pred.num-data.start.t)])
            L.T = sum(alpha.T)
            phi.T = alpha.T/L.T

            for (idx in 1:dist.resolution) {
                pred.distributions[idx, pred.num] = sum((phi.T%*%gamma)*discrete.dnorm(xs[idx], mu, sd, 1))
            }
        }

        return(pred.distributions)
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

    return(list(mle=HMM.mle, mllk=HMM.mllk, 
                local.decoding=HMM.local.decoding,
                data=x, 
                natural.params.to.working=HMM.natural.params.to.working,
                working.params.to.natural=HMM.working.params.to.natural,
                probs.matrix=HMM.get.probs.matrix,
                log.forward.backward=HMM.log.forward.backward,
                cond.state.probs=HMM.cond.state.probs,
                one.step.predictions=HMM.one.step.predictions))
}

discrete.dnorm = function (x, mu, sd, resolution) {
    return(pnorm(x+0.5*resolution, mu, sd) - pnorm(x-0.5*resolution, mu, sd))
}
