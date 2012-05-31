source('src/loaddata.R')

# Hard coded estimated parameters for dataset 1
theta = 0.4615
r0 = 0.1423
K = 822.9574
Q = 0.009
R = 0.0407

a0 = 2.1
a1 = 8.4
m = 250

w = (a1-a0)/m

b.i = a0 + w*(1:m)
p.i = b.i - 0.5*w

state = which.min(abs(p.i-dat1[1]))


for (sim.num in 1:4) {
    p = p.i[state]
    xs = rep(0, 199)
    for (i in 1:199) {
        xs[i] = rnorm(1, mean=p, sd=sqrt(R))
        p = rnorm(1, mean=p + r0*(1-(exp(p)/K)^theta), sd=sqrt(Q))
    }
    pdf(sprintf('plots/sim-%s.pdf', sim.num), 12, 5)
    plot(xs, type="l", lwd=2, main="Simulation of log population",
         xlab="Time", ylab="Log population size")
    dev.off()
}

par(mfrow=c(1,1))
plot(xs, type="l")
plot(dat1, type="l")
