source('src/loaddata.R')
source('src/functions.R')
source('src/model.R')

params = c(0.5, 0.2, 1000, 0.01, 0.05)
model = get.HMM.model(dat1, 250, c(2.1, 8.4))
mllk = model$mllk(params)
