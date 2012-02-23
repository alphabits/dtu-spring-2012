X = cbind(rep(1, 8), c(rep(0,4), rep(1,4)), c(rep(0,4), 1, rep(0,3)))
y = c(4.4, 3.4, 3.3, 2.5, 7.3, 4.9, 4.8, 4.4)

pseudo.inv = solve(t(X)%*%X)
a.hat = pseudo.inv%*%t(X)%*%y
