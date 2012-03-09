# All paths are relative to assignment root dir
source('src/loaddata.R')
source('src/functions.R')

SAVEPLOTS = TRUE

plot.and.save("sales-series.pdf", 12, 7, 1.5, plot,
    x, type="l", xlab="Week", ylab="Weekly soap sale",
    main="Sales of soap product for 242 consecutive weeks")

plot.and.save("acf-data.pdf", 8, 7, 1.5, acf, x, 
              main="ACF of the soap sales data")
