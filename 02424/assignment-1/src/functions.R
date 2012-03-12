plot.and.save = function(filename, width, height, cex, plotfunction, ...) {
    save.the.plot = exists('SAVEPLOTS') && SAVEPLOTS
    if (save.the.plot) {pdf(sprintf('plots/%s', filename), width, height)}
    plotfunction(..., cex.lab=cex, cex.main=cex, cex.axis=cex)
    if (save.the.plot) {dev.off()}
}

float.str = function (float) {
    return(sprintf('%.03f', float))
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

