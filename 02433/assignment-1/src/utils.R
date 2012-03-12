plot.and.save = function(filename, width, height, cex, plotfunction, ...) {
    save.the.plot = exists('SAVEPLOTS') && SAVEPLOTS
    if (save.the.plot) {pdf(sprintf('plots/%s', filename), width, height)}
    plotfunction(..., cex.lab=cex, cex.main=cex, cex.axis=cex)
    if (save.the.plot) {dev.off()}
}

float.str = function (float) {
    return(sprintf('%.03f', float))
}
