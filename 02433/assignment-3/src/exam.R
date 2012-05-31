source('src/loaddata.R')
source('src/hmm-model.R')

library('forecast')

arima.model = auto.arima(power.production)
coefs = arima.model$coef

idx = seq(1, 8000, by=2)

cex = 1.5

for (i in 1:4) {
    xs = arima.sim(list(ar=coefs[1], ma=c(coefs[2], coefs[3]), order=c(1,1,2)), 
                   8000, sd=sqrt(arima.model$sigma2))
    pdf(sprintf('plots/arima-sim-%s.pdf', i), 14, 6)
    plot(as.vector(xs)[idx], cex.lab=cex, cex.main=cex, cex.axis=cex, type="l",
         xlab="Time", ylab="y", main="Simulated dataset")
    dev.off()
}


mu4 = c(18000, 14000, 8000, 2000)
sd4 = c(1200, 1200, 1200, 1000)
gamma4 = matrix(c(rep(c(0.91, rep(0.03, 4)), 3), 0.91), nrow=4)
hmm.model.4 = get.HMM.model(power.production, 4)

timing.4 = system.time(assign("mle.4", hmm.model.4$mle(mu4, sd4, gamma4)))

params  = mle.4$mle.list
gamma = params$gamma
mu = params$mu
sd = params$sd
delta = statdist(gamma)

sim.hmm = function (n, mu, sd, gamma) {
    xs = rep(0, n)
    current.state = 2
    for (i in 1:n) {
        xs[i] = rnorm(1, mean=mu[current.state], sd=sd[current.state])
        current.state = which(rmultinom(1, 1, gamma[current.state,])==1)
    }
    return(xs)
}

for (i in 1:4) {
    xs = sim.hmm(8000, mu, sd, gamma)
    pdf(sprintf('plots/4-state-hmm-sim-%s.pdf', i), 14, 6)
    plot(as.vector(xs)[idx], cex.lab=cex, cex.main=cex, cex.axis=cex, type="l",
         xlab="Time", ylab="y", main="Simulated dataset")
    dev.off()
}

x.plot = seq(0, 21000, length.out=800)
probs = matrix(0, 4, 800)
for (i in 1:4) {
    probs[i,] = dnorm(x.plot, mean=mu[i], sd=sd[i])
}
# Phi is almost uniform for t=8000
pred.dist = matrix(rep(0.25, 4), nrow=1)%*%probs
pdf('plots/one-step-dist.pdf', 12, 6)
plot(x.plot, pred.dist, type="l", ylab="p(x)", xlab="Power production",
     cex.lab=cex, cex.axis=cex, cex.main=cex, main="One-step prediction distribution of HMM",
     lwd=2)
abline(v=9734.5)
dev.off()
