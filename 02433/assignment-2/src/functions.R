plot.and.save = function(filename, width, height, cex, plotfunction, ...) {
    save.the.plot = exists('SAVEPLOTS') && SAVEPLOTS
    if (save.the.plot) {pdf(sprintf('plots/%s', filename), width, height)}
    plotfunction(..., cex.lab=cex, cex.main=cex, cex.axis=cex)
    if (save.the.plot) {dev.off()}
}

float.str = function (float) {
    return(sprintf('%.03f', float))
}

save.mle.res = function (mle.res, prefix) {
    save.latex.table.content(sprintf('%s-estimate.tex', prefix), matrix(mle.res$mle, nrow=1),
                             "\\num{%.02e}")
    save.digit(sprintf('%s-mllk.tex', prefix), mle.res$mllk)
    save.digit(sprintf('%s-aic.tex', prefix), mle.res$AIC)
    save.digit(sprintf('%s-bic.tex', prefix), mle.res$BIC)
    if (!is.null(mle.res$var.cov)) {
        save.latex.table.content(sprintf('%s-varcov.tex', prefix), mle.res$var.cov, 
                                 "\\num{%.02e}")
        sds = sqrt(diag(mle.res$var.cov))
        lower.bound = mle.res$mle - 1.96*sds
        upper.bound = mle.res$mle + 1.96*sds
        estimate.table = cbind(mle.res$mle, lower.bound, upper.bound, sds)
        save.estimate.table(sprintf('%s-confidence.tex', prefix), estimate.table)
    }
}

save.digit = function (filename, digit, digit.format='%.03f') {
    sink(sprintf('results/%s', filename))
    cat(sprintf(digit.format, digit))
    sink()
}

get.estimate.table = function (mat) {
    row.output = c()
    row.labels = c("$\\theta$", "$r_0$", "$K$", "$Q$", "$R$")
    for (i in 1:dim(mat)[1]) {
        row.output[i] = paste(row.labels[i], "&", paste(sprintf("%.04f", mat[i,]), collapse=" & "))
    }
    return(paste(row.output, collapse="\\\\"))
}

save.estimate.table = function (filename, mat) {
    sink(sprintf('results/%s', filename))
    cat(get.estimate.table(mat))
    sink()
}

latex.table.content = function (mat, digit.format='%.03f') {
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

save.decoding.plot = function (decoding, prefix, main) {
    plot.and.save(sprintf('%s-decoding.pdf', prefix), 12, 5, 1.5, 
                  plot.decoding, decoding, main)
}

plot.decoding = function (decoding, main, ...) {
    plot(decoding$model$data, type="l", xlab="Time", col="gray", 
         ylab="Log-population size", main=main, lwd=2, ...)
    lines(decoding$continous, lwd=2)
    legend("bottomright", c("Observations", "Decoding"), lty=c(1,1), 
           col=c("gray", "black"), lwd=c(2.5, 2.5), inset=c(0.02, 0.05))
}
