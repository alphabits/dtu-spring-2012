theta = 0.5
r0 = 0.15
K = 823
Q = 0.00905
R = 0.0407

ps = c(3)
xs = c(3+rnorm(1, sd=R))

for (i in 2:250) {
    ps[i] = ps[i-1] + r0*(1 - (exp(ps[i-1])/K)^theta) + rnorm(1, sd=sqrt(Q))
    xs[i] = ps[i] + rnorm(1, sd=sqrt(R))
}
