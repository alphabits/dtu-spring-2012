source('src/functions.R')
source('src/loaddata.R')

SAVEPLOTS = TRUE

for (i in 1:2) {
    filename = sprintf('dataset-%s.pdf', i)
    main = paste("Noisy observations of log-population size of population", i)
    plot.and.save(filename, 12, 5, 1.5, plot, dat[[i]], type="l", lwd=2,
                  main=main, xlab="Time", ylab="Log-population size")
}
