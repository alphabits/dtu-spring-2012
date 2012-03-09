plot.and.save = function(filename, width, height, cex, plotfunction, ...) {
    save.the.plot = exists('SAVEPLOTS') && SAVEPLOTS
    if (save.the.plot) {pdf(sprintf('plots/%s', filename), width, height)}
    plotfunction(..., cex.lab=cex, cex.main=cex, cex.axis=cex)
    if (save.the.plot) {dev.off()}
}

get.initial.gamma = function (m, off.diag=0.1) {
    on.diag = 1 - (m-1)*off.diag
    dat = c()
    for (column in 1:m) {
        if (column > 1) {
            dat = c(dat, rep(off.diag, column-1))
        }
        dat = c(dat, on.diag)
        if (column < m) {
            dat = c(dat, rep(off.diag, m-column))
        }
    }
    return(matrix(dat, nrow=m))
}

float.str = function (float) {
    return(sprintf('%.03f', float))
}

statdist = function (gamma) {
    m = dim(gamma)[1]
    one.row = rep(1, m)
    U = matrix(rep(one.row, m), nrow=m)
    Id = diag(one.row)
    inv = solve(Id - gamma + U)
    return(one.row %*% inv)
}

save.fit.results = function (fit, label) {
    for (key in c('delta', 'gamma', 'lambda')) {
        save.latex.matrix(sprintf('%s-%s.tex', label, key), fit[[key]])
    }
}

repeated.fit.results.table = function (fits) {
    num.fits = length(fits)
    output.rows = rep("", num.fits)
    for (i in 1:num.fits) {
        fit = fits[[i]]
        output.rows[i] = sprintf('$%.02f$ & $%s$ & $%s$ & $%.02f$',
                                 fit$gamma.init[2,1],
                                 latex.matrix(fit$lambda.init, '%i'),
                                 latex.matrix(fit$lambda),
                                 fit$mllk)
    }
    return(paste(output.rows, collapse=" \\\\\n"))
}

save.repeated.fit.results = function (fits, label) {
    sink(sprintf('results/%s-fit-table.tex', label))
    cat(repeated.fit.results.table(fits))
    sink()
}

save.comparison.table = function (filename, sample, fits) {
    path = sprintf('results/%s', filename)
    sink(path)
    cat(get.comparison.table(sample, fits))
    sink()
}

get.comparison.table = function (sample, fits) {
    sample.corr = acf(sample, plot=F)[[1]]
    sample.mean = mean(sample)
    sample.var = var(sample)
    output.rows = c(get.comparison.table.row('Sample', sample.corr, sample.mean,
                                              sample.var))
    for (fit.num in 1:length(fits)) {
        label = paste(fit.num+1, '-states', sep="")
        fit = fits[[fit.num]]
        output.rows = c(output.rows, 
                        get.comparison.table.row(label, 
                                                 fit$moments$correlations,
                                                 fit$moments$expectation,
                                                 fit$moments$variance))
    }
    return(paste(output.rows, collapse="\n"))
}

get.comparison.table.row = function (label, corrs, mean, var) {
    return(paste(
        label, '&', 
        float.str(mean), '&',
        float.str(var), '&',
        paste(float.str(corrs[2:8]), collapse=" & "), "\\\\"
    ))
}

get.information.content.table = function (fits) {
    rows = c()
    for (fit.num in 1:length(fits)) {
        fit = fits[[fit.num]]
        label = paste(fit.num+1, '-states', sep='')
        rows = c(rows, paste(label, '&', 
                             float.str(fit$AIC), '&', 
                             float.str(fit$BIC), "\\\\"))
    }
    return(paste(rows, collapse="\n"))
}

save.information.content.table = function (filename, fits) {
    path = sprintf('results/%s', filename)
    sink(path)
    cat(get.information.content.table(fits))
    sink()
}

plot.fit.results = function (x, fit, label) {
    lambda = fit$lambda
    gamma = fit$gamma
    plot.and.save(sprintf('%s-%s.pdf', label, 'histogram'), 12, 9, 1.5,
                  plot.HMM.with.hist, x, lambda, gamma, main="", xlab="",
                  ylab="")
    plot.and.save(sprintf('%s-%s.pdf', label, 'means-on-data'), 12, 9, 1.5,
                  plot.data.with.means, x, lambda, main="", xlab="", 
                  ylab="")
}

latex.matrix = function (mat, digit.format='%.03f') {
    if (!class(mat) == "matrix") {
        mat = matrix(mat, ncol=length(mat))
    }
    return(paste("\\begin{pmatrix}", 
                 latex.table.content(mat, digit.format),
                 "\\end{pmatrix}",
                 sep="\n"))
}

latex.table.content = function (mat, digit.format='%.03f') {
    row.output = c()
    for (i in 1:dim(mat)[1]) {
        row.output[i] = paste(sprintf(digit.format, mat[i,]), collapse=" & ")
    }
    return(paste(row.output, collapse="\\\\"))
}

save.latex.matrix = function (filename, mat, digit.format='%.03f') {
    sink(sprintf('results/%s', filename))
    cat(latex.matrix(mat, digit.format))
    sink()
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

p.HMM = function (x, lambda, gamma, delta=FALSE) {
    if (delta == FALSE) {
        delta = statdist(gamma)
    }

    return(sum(dpois(x, lambda=lambda)*delta))
}

bootstrap.errors = function (x, initial.fit, num.samples, lambda.0, gamma.0) {
    sample.size = length(x)
    num.states = length(initial.fit$lambda)
    lambdas = matrix(0, num.samples, num.states)
    trans.probs = matrix(0, num.samples, num.states^2)
    for (sample.num in 1:num.samples) {
        samp = pois.HMM.generate_sample(sample.size, num.states, 
                                        initial.fit$lambda, 
                                        initial.fit$gamma)
        samp.fit = pois.HMM.mle(samp, num.states, lambda.0, gamma.0)
        lambdas[sample.num,] = samp.fit$lambda
        trans.probs[sample.num,] = c(samp.fit$gamma)
    }
    return(list(lambda=cov(lambdas), trans.prob=cov(trans.probs)))
}

# From exercise 9 in week 2
pois.HMM.moments = function (m, lambda, gamma, lag.max=10) {
    delta = statdist(gamma)
    lambda = matrix(lambda)
    Lambda = diag(c(lambda))
    expectation = delta %*% lambda
    expectation.sq = expectation^2
    M = delta %*% Lambda %*% lambda
    variance = expectation + M - expectation.sq
    correlations = matrix(c(1, rep(0, lag.max)), nrow=1)
    tmp.gamma = gamma
    for (i in 2:(lag.max+1)) {
        cov = delta %*% Lambda %*% tmp.gamma %*% lambda - expectation.sq
        correlations[i] = cov/variance
        tmp.gamma = tmp.gamma %*% gamma
    }

    return(list(
        expectation=expectation,
        variance=variance,
        correlations=correlations
    ))
}


calculate.fit.result = function (fit.func, lambdas, off.diags, m, ...) {
    fit.res = list()
    for (lambda in lambdas) {
        for (off.diag in off.diags) {
            gamma.init = get.initial.gamma(m, off.diag)
            fit = fit.func(x, m, lambda, gamma.init, ...)
            fit$lambda.init = lambda
            fit$gamma.init = gamma.init
            fit$moments = pois.HMM.moments(m, fit$lambda, fit$gamma)
            fit.res[[(length(fit.res)+1)]] = fit
        }
    }
    return(fit.res)
}

calculate.fit.result.ml = function (lambdas, off.diags, m) {
    return(calculate.fit.result(pois.HMM.mle, lambdas, off.diags, m))
}

calculate.fit.result.em = function (lambdas, off.diags, m, ...) {
    return(calculate.fit.result(pois.HMM.EM, lambdas, off.diags, m, ...))
}
