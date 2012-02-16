x = 1:10
y = c(3, 0, 4, 5, 6, 4, 9, 7, 4, 10)

l = function (theta) {
    -sum(dpois(y, exp(theta[1]*x+theta[2]), log=TRUE))
}
fit = optim(par=c(0,0), fn=l, hessian=TRUE)
