params = c(0.5, 0.2, 1000, 0.01, 0.05)

model = get.HMM.model(dat1, 250, c(2.1, 8.4))
mllk = model$mllk(params)
