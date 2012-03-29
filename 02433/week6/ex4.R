source('../assignment-1/src/A1.R')
source('../assignment-1/src/A2.R')

gamma = rbind(c(0.6, 0.4),
              c(0.4, 0.6))

lambda1 = c(2, 5)
lambda2 = c(2, 7)

samp = pois.HMM.generate_sample(100, 2, lambda1, gamma)
samp = c(samp, pois.HMM.generate_sample(20, 2, lambda2, gamma))
