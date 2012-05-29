source('src/loaddata.R')

model.1 = auto.arima(power.production)
model.2 = arima(power.production, order=c(6,0,0))
model.3 = arima(power.production, order=c(5,0,0))

calc.setup = list(
    list(model=model.1, predictions=rep(0,2000), errors=rep(0,2000)),
    list(model=model.2, predictions=rep(0,2000), errors=rep(0,2000)),
    list(model=model.3, predictions=rep(0,2000), errors=rep(0,2000))
)

for (k in 1:num.test.samples) {
    idx = num.training.samples + k
    for (model.num in 1:length(calc.setup)) {
        setup = calc.setup[[model.num]]
        tmp.prediction = predict(setup$model, n.ahead=1)
        calc.setup[[model.num]]$model = Arima(dat[1:idx,3], model=setup$model)
        calc.setup[[model.num]]$predictions[k] = tmp.prediction$pred[1]
        calc.setup[[model.num]]$errors[k] = tmp.prediction$se[1]
    }
}

qqplot.fns = function(x, ...) {
    qqnorm(x, ...)
    qqline(x)
}

# Forgot the labels and prediction error calculations in calc.setup
org.models = list(model.1, model.2, model.3)
model.labels = c("auto-arima", "arima600", "arima500")
model.descriptions = c("ARIMA(1,1,2)", "ARIMA(6,0,0)", "ARIMA(5,0,0)")
prediction.errors = rep(0, length(calc.setup))

for (i in 1:3) {
    lbl = model.labels[i]
    desc = model.descriptions[i]
    setup = calc.setup[[i]]
    tmp.R = sqrt(sum((setup$predictions - power.production.test)^2)/1999)
    prediction.errors[i] = tmp.R
    save.result(sprintf('prediction-error-%s.txt', lbl), cat, round(tmp.R, 2))

    save.result(sprintf('%s-model.txt', lbl), print, org.models[[i]])

    plot.and.save(sprintf('acf-residuals-%s.pdf', lbl), 8, 7, 1.5, acf, residuals(org.models[[i]]),
                  main=sprintf("ACF for residuals for %s model", desc))
    plot.and.save(sprintf('rstandard-%s.pdf', lbl), 12, 7, 1.5, plot, rstandard(org.models[[i]]),
                  xlab="Index", ylab="Standardized residual", type="l",
                  main=sprintf("Standardized residuals for the %s model", desc))


    plot.and.save(sprintf('qq-plot-%s.pdf', lbl), 7, 7, 1.5, qqplot.fns, residuals(org.models[[i]]), 
                  main=sprintf("QQ plot of residuals for %s model", desc))
}


plot.and.save('acf-diffed-data.pdf', 7, 7, 1.5, acf, diff(power.production),
              main="ACF for first difference of training data")
plot.and.save('pacf-diffed-data.pdf', 7, 7, 1.5, pacf, diff(power.production),
              main="PACF for first difference of training data")
