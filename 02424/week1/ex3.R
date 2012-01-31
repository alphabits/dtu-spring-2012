x = c(13, 5, 28, 28, 15, 4, 13, 4, 10, 17, 11, 13, 12, 17, 3)

minus.likelihood = function (lambda) {
    -sum(x*log(1-lambda) + log(lambda))
}

diff.likelihood = function (lambda) {
    length(x)/lambda - sum(x/(1-lambda))
}

estimate = optimize(f=minus.likelihood, interval=c(0.01, 0.99))
