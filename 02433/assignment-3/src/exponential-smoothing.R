source('src/loaddata.R')
source('src/functions.R')

library('forecast')


auto.fit = ets(power.production)
sink('results/exponential-smoothing.txt')
print(auto.fit)
sink()

R.random.walk = sqrt(sum(diff(power.production.test)^2)/1998)
sink('results/prediction-error-random-walk.txt')
cat(sprintf("%0.2f", R.random.walk))
sink()
