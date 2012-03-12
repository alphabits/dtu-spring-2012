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

get.ys.for.state.plot = function (states, lambda) {
    ys = rep(0, length(states))
    for (state.num in 1:length(lambda)) {
        active.points = states==state.num
        ys = ys + lambda[state.num]*active.points
    }
    return(ys)
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

get.state.prediction.table = function (preds) {
    output = c()
    for (i in 1:dim(preds)[1]) {
        output[i] = paste(i, " & ", 
                          paste(float.str(preds[i,]), collapse=" & "))
    }
    return(paste(output, collapse="\\\\"))
}
