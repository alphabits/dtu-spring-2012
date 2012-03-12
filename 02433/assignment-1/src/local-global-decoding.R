# This source file assumes that src/fit-em.R has
# been run.


local.decoding = pois.HMM.local_decoding(x, 3, fit.3$lambda, fit.3$gamma,
                                         fit.3$delta)

global.decoding = pois.HMM.viterbi(x, 3, fit.3$lambda, fit.3$gamma,
                                   fit.3$delta)

plot.and.save('local-decoding.pdf', 12, 7, 1.5, plot.predicted.states,
              x, local.decoding, fit.3$lambda, 
              main="Local decoding for EM fitted 3-state Poisson HMM",
              xlab="Week", ylab="Weekly soap sale")

plot.and.save('global-decoding.pdf', 12, 7, 1.5, plot.predicted.states,
              x, global.decoding, fit.3$lambda, 
              main="Global decoding for EM fitted 3-state Poisson HMM",
              xlab="Week", ylab="Weekly soap sale")

state.preds = pois.HMM.state_prediction(x, 3, fit.3$lambda, fit.3$gamma,
                                        H=10)
