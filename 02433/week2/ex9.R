# Function taken from 
# http://134.76.173.220/hmm-with-r/appendix-a/A2.txt 
pois.HMM.generate_sample <-                               
 function(n,m,lambda,gamma,delta=NULL)                     
{                                                           
 if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))    
 mvect <- 1:m                                               
 state <- numeric(n)                                        
 state[1] <- sample(mvect,1,prob=delta)                     
 for (i in 2:n)                                             
   state[i]<-sample(mvect,1,prob=gamma[state[i-1],])        
 x <- rpois(n,lambda=lambda[state])                      
 x                                                          
}


# Start exercises 9 and 11

statdist = function (gamma) {
    m = dim(gamma)[1]
    one.row = rep(1, m)
    U = matrix(rep(one.row, m), nrow=m)
    Id = diag(one.row)
    inv = solve(Id - gamma + U)
    return(one.row %*% inv)
}

pois.HMM.moments = function (m, lambda, gamma, lag.max=10) {
    delta = statdist(gamma)
    Lambda = diag(c(lambda))
    expectation = delta %*% t(lambda)
    expectation.sq = expectation^2
    M = delta %*% Lambda %*% t(lambda)
    variance = expectation + M - expectation.sq
    correlations = matrix(c(1, rep(0, lag.max)), nrow=1)
    tmp.gamma = gamma
    for (i in 2:(lag.max+1)) {
        print(i)
        cov = delta %*% Lambda %*% tmp.gamma %*% t(lambda) - expectation.sq
        correlations[i] = cov/variance
        tmp.gamma = tmp.gamma %*% gamma
    }

    return(list(
        expectation=expectation,
        variance=variance,
        correlations=correlations
    ))
}

test.moment.function = function () {
    gamma = rbind(c(0.8,0.2,0),c(0.05,0.7,0.25),c(0.1,0.25,0.65))
    lambda = t(matrix(c(1,5,10)))
    m = length(lambda)
    lag.max = 15
    mom = pois.HMM.moments(m, lambda, gamma, lag.max)
    samp = pois.HMM.generate_sample(10000, m, lambda, gamma)
    mome = list(
        expectation=mean(samp),
        variance=var(samp),
        correlations=acf(samp)
    )
    return(list(mom=mom, mome=mome))
}
