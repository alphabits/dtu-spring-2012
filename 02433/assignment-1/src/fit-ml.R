# All paths are relative to assignment root
source('src/loaddata.R')
source('src/functions.R')
source('src/A1.R')
source('src/A2.R')

SAVEPLOTS = TRUE

off.diag.to.test = c(0.01, 0.05, 0.1)


lambdas.2 = list(c(2,8), c(3,7), c(4,6))
fit.2.res = calculate.fit.result.ml(lambdas.2, off.diag.to.test, 2)
# One fit chosen for the report
fit.2 = fit.2.res[[3]]

lambdas.3 = list(c(1,5,9), c(2,5,8), c(3,5,7))
fit.3.res = calculate.fit.result.ml(lambdas.3, off.diag.to.test, 3)
# One fit chosen for the report
fit.3 = fit.3.res[[6]]

lambdas.4 = list(c(1, 4, 6, 9), c(2, 4, 6, 8), c(1, 4, 8, 12))
fit.4.res = calculate.fit.result.ml(lambdas.4, off.diag.to.test, 4)
# One fit chosen for the report
fit.4 = fit.4.res[[1]]
