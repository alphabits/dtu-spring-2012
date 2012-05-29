source('src/loaddata.R')
source('src/functions.R')

SAVEPLOTS = TRUE


plot.trainingdata = function (idx, ...) {
    first.day.shift = day.shifts[which(day.shifts>=idx[1])[1]]
    date.idx = seq(first.day.shift, idx[length(idx)], by=2*measurements.pr.day)
    plot(dates[idx], power.production[idx], xlab="Date", ylab="Power production in kW", 
         type="l", main="Power production at Klim wind farm at 5 minute resolution", 
         xaxt="n", ...)
    axis(1, at=dates[date.idx], labels=format(dates[date.idx], "%d %b"), 
         cex.axis=list(...)$cex.axis)
}

plot.and.save('training-dataset.pdf', 14, 6, 1.5, plot.trainingdata, 1:8000)
plot.and.save('training-dataset-1.pdf', 14, 6, 1.5, plot.trainingdata, 1:4000)
plot.and.save('training-dataset-2.pdf', 14, 6, 1.5, plot.trainingdata, 4001:8000)

plot.and.save('acf-training-dataset.pdf', 7, 7, 1.5, acf, power.production,
              main="ACF for training dataset")
plot.and.save('pacf-training-dataset.pdf', 7, 7, 1.5, pacf, power.production,
              main="PACF for training dataset")

sink('results/data-summary-training-data.txt')
print(summary(power.production))
sink()

plot(power.production[1:1000], type="l", main="Testing my setup which is nice")
res = ets(power.production)

