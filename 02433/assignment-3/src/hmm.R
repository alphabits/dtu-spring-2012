source('src/loaddata.R')
source('src/functions.R')
source('src/hmm-model.R')


hmm.model.2 = get.HMM.model(power.production, 2)
hmm.model.3 = get.HMM.model(power.production, 3)
hmm.model.4 = get.HMM.model(power.production, 4)

mu2 = c(4000, 9000)
sd2 = c(2000, 4000)
gamma2 = matrix(c(0.9, 0.1, 0.1, 0.9), nrow=2)

mu3 = c(16000, 10000, 4000)
sd3 = c(1800, 1800, 2000)
gamma3 = matrix(c(rep(c(0.9, rep(0.05, 3)), 2), 0.9), nrow=3)

mu4 = c(18000, 14000, 8000, 2000)
sd4 = c(1200, 1200, 1200, 1000)
gamma4 = matrix(c(rep(c(0.91, rep(0.03, 4)), 3), 0.91), nrow=4)

timing.2 = system.time(assign("mle.2", hmm.model.2$mle(mu2, sd2, gamma2)))
timing.3 = system.time(assign("mle.3", hmm.model.3$mle(mu3, sd3, gamma3)))
timing.4 = system.time(assign("mle.4", hmm.model.4$mle(mu4, sd4, gamma4)))

M.2 = get.M.from.natural.params(mle.2$mle.list)
var.cov.2 = inv.hessian.from.working.hessian(mle.2$nlm.res$hessian, M.2)
M.3 = get.M.from.natural.params(mle.3$mle.list)
var.cov.3 = inv.hessian.from.working.hessian(mle.3$nlm.res$hessian, M.3)
M.4 = get.M.from.natural.params(mle.4$mle.list)
var.cov.4 = inv.hessian.from.working.hessian(mle.4$nlm.res$hessian, M.4)

get.decoding = function (model, params) {
    decod = model$local.decoding(params)
    m = length(params$mu)
    for (i in 1:m) {
        decod[decod==i] = params$mu[i]
    }
    return(decod)
}

hmm.res.2 = list(
    time.info=timing.2,
    initial.values=list(mu=mu2, sd=sd2, gamma=gamma2),
    mle=mle.2,
    var.cov=var.cov.2,
    decoding=get.decoding(hmm.model.2, mle.2$mle.list)
)

hmm.res.3 = list(
    time.info=timing.3, 
    initial.values=list(mu=mu3, sd=sd3, gamma=gamma3),
    mle=mle.3,
    var.cov=var.cov.3,
    decoding=get.decoding(hmm.model.3, mle.3$mle.list)
)

hmm.res.4 = list(
    time.info=timing.4, 
    initial.values=list(mu=mu4, sd=sd4, gamma=gamma4),
    mle=mle.4,
    var.cov=var.cov.4,
    decoding=get.decoding(hmm.model.4, mle.4$mle.list)
)

save.mle.res(hmm.res.2, '2-state-normal')
save.mle.res(hmm.res.3, '3-state-normal')
save.mle.res(hmm.res.4, '4-state-normal')
