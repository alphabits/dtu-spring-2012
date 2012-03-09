# All paths are relative to assignment root
source('src/loaddata.R')
source('src/functions.R')
source('src/A1.R')
source('src/A2.R')

SAVEPLOTS = TRUE

off.diag.to.test = c(0.01, 0.05, 0.1)


lambdas.2 = list(c(2,8), c(3,7), c(4,6))
fit.2.res = calculate.fit.result.em(lambdas.2, off.diag.to.test, 2, 
                                    delta=c(0.5,0.5))
# One fit chosen for the report

lambdas.3 = list(c(1,5,9), c(2,5,8), c(3,5,7))
fit.3.res = calculate.fit.result.em(lambdas.3, off.diag.to.test, 3, 
                                    delta=c(1/3,1/3,1/3))
# One fit chosen for the report

lambdas.4 = list(c(1, 4, 6, 9), c(2, 4, 6, 8), c(1, 4, 8, 12))
fit.4.res = calculate.fit.result.em(lambdas.4, off.diag.to.test, 4, 
                                    delta=c(0.25,0.25,0.25,0.25))
# One fit chosen for the report
