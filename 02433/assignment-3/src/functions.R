plot.and.save = function(filename, width, height, cex, plotfunction, ...) {
    save.the.plot = exists('SAVEPLOTS') && SAVEPLOTS
    if (save.the.plot) {pdf(sprintf('plots/%s', filename), width, height)}
    plotfunction(..., cex.lab=cex, cex.main=cex, cex.axis=cex)
    if (save.the.plot) {dev.off()}
}

save.result = function (filename, printfunction, ...) {
    sink(sprintf('results/%s', filename))
    printfunction(...)
    sink()
}


float.str = function (float) {
    return(sprintf('%.03f', float))
}

estimate.vector.to.list = function (v, m) {
    return(list(
        mu=v[1:m],
        sd=v[(m+1):(2*m)],
        gamma=matrix(v[(2*m+1):length(v)], nrow=m)
    ))
}

save.mle.res = function (mle.res, prefix) {
    m = length(mle.res$initial.values$mu)
    for (param.name in c("gamma", "sd", "mu")) {
        estimate = mle.res$mle$mle.list[[param.name]]
        initial = mle.res$initial.values[[param.name]]
        save.latex.table.content(sprintf('%s-%s-estimate.tex', prefix, param.name), 
                                 estimate, "\\num{%.02e}")
        save.latex.table.content(sprintf('%s-%s-initial.tex', prefix, param.name), 
                                 initial, "\\num{%.02e}")
    }
    save.result(sprintf('%s-timing.tex', prefix), print, mle.res$time.info)
    save.digit(sprintf('%s-mllk.tex', prefix), mle.res$mle$mllk)
    save.digit(sprintf('%s-aic.tex', prefix), mle.res$mle$AIC)
    save.digit(sprintf('%s-bic.tex', prefix), mle.res$mle$BIC)
    delta = statdist(mle.res$mle$mle.list$gamma)
    save.latex.table.content(sprintf('%s-statdist.tex', prefix), delta, "\\num{%.02e}")
    save.state.dependent.plot(mle.res, prefix)
    save.decoding(mle.res$decoding, prefix)
    if (!is.null(mle.res$var.cov)) {
        save.latex.table.content(sprintf('%s-varcov.tex', prefix), mle.res$var.cov, 
                                 "\\num{%.02e}")
        sds = sqrt(sapply(diag(mle.res$var.cov), max, 0))
        lower.bound = mle.res$mle$mle - 1.96*sds
        upper.bound = mle.res$mle$mle + 1.96*sds
        estimate.table = cbind(mle.res$mle$mle, lower.bound, upper.bound, sds)
        save.estimate.table(sprintf('%s-confidence.tex', prefix), estimate.table, m)
    }
}

get.estimate.table = function (mat, m) {
    row.output = c()
    gamma.idx = c(outer(1:m, 1:m, function (x, y) {paste(y, x, sep="")}))
    row.labels = c(sprintf("$\\mu_%s$", 1:m), sprintf("$\\sigma_{%s}$", 1:m),
                   sprintf("$\\gamma_{%s}$", gamma.idx))
    for (i in 1:dim(mat)[1]) {
        row.output[i] = paste(row.labels[i], "&", paste(sprintf("%.04f", mat[i,]), collapse=" & "))
    }
    return(paste(row.output, collapse="\\\\"))
}

save.estimate.table = function (filename, mat, m) {
    sink(sprintf('results/%s', filename))
    cat(get.estimate.table(mat, m))
    sink()
}

save.decoding = function (decoding, prefix) {
    pdf(sprintf('plots/%s-decoding.pdf', prefix), 12, 7)
    plot(power.production, type="l", col="grey", xlab="Time", ylab="Power production",
         main="Local decoding")
    lines(decoding)
    dev.off()
}

save.state.dependent.plot = function (mle.res, prefix) {
    mu.est = mle.res$mle$mle.list$mu
    sd.est = mle.res$mle$mle.list$sd
    m = length(mu.est)
    xs = seq(-2000, 22000, length.out=2000)
    distributions = matrix(0, 2000, m)
    for (i in 1:m) {
        distributions[,i] = dnorm(xs, mu.est[i], sd.est[i])
    }
    rng = range(distributions)
    ymax = rng[2] + 0.1*diff(rng)
    pdf(sprintf('plots/%s-state-dependent-plot.pdf', prefix), 12, 7)
    plot(xs, distributions[,1], type="l", ylim=c(0,ymax),
         main="State dependent normal distributions", xlab="x", ylab="p(x|c)",
         lwd=2)
    for (i in 2:m) {
        lines(xs, distributions[,i], lwd=2)
    }
    dev.off()

    pdf(sprintf('plots/%s-dist-plot.pdf', prefix), 12, 7)
    plot(power.production, type="l", col="gray", xlab="Time", ylab="Power production",
         main="State dependent distributions on data plot")
    for (i in 1:m) {
        lines(distributions[,i]*1e7, xs, lwd=2)
        abline(h=mu.est[i])
    }
    dev.off()
}

save.digit = function (filename, digit, digit.format='%.03f') {
    sink(sprintf('results/%s', filename))
    cat(sprintf(digit.format, digit))
    sink()
}


latex.table.content = function (mat, digit.format='%.03f') {
    if (class(mat) == "numeric") {
        mat = matrix(mat)
    }
    row.output = c()
    for (i in 1:dim(mat)[1]) {
        row.output[i] = paste(sprintf(digit.format, mat[i,]), collapse=" & ")
    }
    return(paste(row.output, collapse="\\\\"))
}

save.latex.table.content = function (filename, mat, digit.format='%.03f') {
    sink(sprintf('results/%s', filename))
    cat(latex.table.content(mat, digit.format))
    sink()
}
