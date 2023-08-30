# cfa analogon to jags_dyn.R model
model{

############################
#
# Global
#
############################
# Constants
  for(k in 1:10){
    ly[(k-1)*10+1] <- 1
    ty1[(k-1)*10+1] <- muy0}

  for(i in 1:N){
  # i = person

  ##############################################
  # within level
  ##############################################


    for(j in 1:100){
      # j = item
      y[i,j] ~ dnorm(muy[i,j],psiy[j])
      muy[i,j] <- ty1[j] + ly[j]*xi[i,index[j]]
      }

    # factors
    # psixi is precision matrix
    # muxi is mean vector
    xi[i,1:10] ~ dmnorm(muxi,psixi)

    }

  ##############################################
  # priors
  ##############################################

  for(j in 1:100){
    psiy[j] ~ dgamma(4,4) }


  for(j in 1:90){
    ly[index3[j]]  ~ dnorm(1,1)T(0,)
    ty1[index3[j]] ~ dnorm(muy0,1) }

  for(j in 1:10){muxi[j] = 0}
  psixi[1:10,1:10] ~ dwish(R0[1:10,1:10],10)

  ##############################################
  # LL
  ##############################################
  for(i in 1:N){
    for(j in 1:100){
      log_lik0[i,j] <- logdensity.norm(y[i,j],muy[i,j],psiy[j])
    }
    log_lik[i] <- sum(log_lik0[i,])
  }

  ##############################################
  # transformations
  ##############################################
  for(j in 1:100){sigmay[j]<- 1/psiy[j] } #zac: right?
  sigmaxi[1:10,1:10] <- inverse(psixi[1:10,1:10])

  ##############################################
  # STD loadings and covariances
  ##############################################
  for(j in 1:100){
    vx[j] <- ly[j]^2*sigmaxi[index[j],index[j]] + sigmay[j]
    lySt[j] <- ly[j]*sqrt(sigmaxi[index[j],index[j]] / vx[j])  }

  for(q in 1:90){
    Fcor[index2[q,1],index2[q,2]] <- sigmaxi[index2[q,1],index2[q,2]] / sqrt(sigmaxi[index2[q,1],index2[q,1]]*sigmaxi[index2[q,2],index2[q,2]]) }

}
