plot.predicted.states = function (x, states, lambda, ...) {
    plot.ys = get.ys.for.state.plot(states, lambda)
    plot(x, type="l", ...)
    points(plot.ys, pch=20)
}

plot.HMM = function (lambda, gamma, max.x=25, plotfns=plot, ...) {
    delta = statdist(gamma)
    x = 0:max.x
    y = sapply(x, p.HMM, lambda=lambda, gamma=gamma, delta=delta)
    plotfns(x, y, ...)
}

plot.HMM.with.hist = function (x, lambda, gamma, ...) {
    hist(x, freq=F, ylim=c(0, 0.16), col="gray", border="white", ...)
    plot.HMM(lambda, gamma, plotfns=lines, lwd=3)
}

plot.data.with.means = function (x, means, ...) {
    plot(x, type="l", ...)
    for (m in means) {
        lines(c(-1000, 1000), rep(m, 2), lwd=3)
    }
}
