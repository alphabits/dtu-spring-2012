# Include statdist
source('../utils.R')

x = c(2, 8, 6, 36, 1, 0, 0, 4, 7)
n = length(x)
gamma = matrix(c(0.9, 0.2, 0.1, 0.8), nrow=2)
delta = statdist(gamma)
lambda = c(1, 5)

get.B = function (x) {
    gamma %*% diag(lambda^x*exp(-lambda)/factorial(x))
}


# No scaling

alpha = delta
alphas = matrix(0, n, 2)

for (i in 1:n) {
    alpha = alpha %*% get.B(x[i])
    alphas[i,] = alpha
}


# Scaling

phi = delta
phis = matrix(0, n, 2)
l = 0

for (i in 1:n) {
    v = phi %*% get.B(x[i])
    u = sum(v)
    l = l + log(u)
    phi = v/u
    phis[i,] = phi
}
