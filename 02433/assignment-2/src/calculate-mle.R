source('src/loaddata.R')
source('src/functions.R')
source('src/model.R')

# Define initial parameters to test and define
# a model object for each dataset
params = list(c(0.5, 0.2, 1000, 0.01, 0.05), 
              c(0.8, 0.8, 100, 1, 2))
models = list(get.HMM.model(dat1, 250, c(2.1, 8.4)),
              get.HMM.model(dat2, 250, c(2.1, 8.4)))
mles = list()

# Just a helper function
get.label = function (param.num, model.num) {
    return(sprintf('cont-dataset-%s-param-%s', model.num, param.num))
}

# Estimate parameters for each combination of initial parameters
# and models
for (param.num in 1:length(params)) {
    for (model.num in 1:length(models)) {
        label = get.label(param.num, model.num)
        mles[[label]] = models[[model.num]]$mle(params[[param.num]], 
                                                print.level=2)
    }
}

# Save the results in .tex files to include in report
for (k in names(mles)) {
    save.mle.res(mles[[k]], k)
}

# Calculate local decodings and save plots as pdf
decodings = list()
for (model.num in 1:length(models)) {
    params = mles[[get.label(1, model.num)]]$mle
    model = models[[model.num]]
    discrete.decoding = model$local.decoding(params)
    cont.decoding = model$p.i[discrete.decoding]
    decodings[[model.num]] = list(discrete=discrete.decoding,
                                  continous=cont.decoding,
                                  model=model)
}
for (model.num in 1:length(models)) {
    save.decoding.plot(decodings[[model.num]], get.label(1, model.num),
                       paste("Local decoding of population", model.num))
}

# Calculate the estimated correlations between parameter estimates
vc = mles[[get.label(1, 2)]]$var.cov
for (i in 1:4) {
    for (j in (i+1):5) {
        print(c(i, j))
        print(vc[i, j]/sqrt(vc[i, i]*vc[j, j]))
    }
}
