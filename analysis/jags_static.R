# static analogon to jags_dyn.R model
# jags_static2: p0 instead of logit construction for state prior
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
      y[i,j] ~ dnorm(muy[i,j,C[i]],psiy_ij[i,j,C[i]])

      # attention state
      muy[i,j,1] <- ty1[j] + ly[j]*xi[i,index[j]]
      psiy_ij[i,j,1] <- psiy[j]
      # non-attentive state
      muy[i,j,2] <- ty2
      psiy_ij[i,j,2] <- psiy_i[i]
    }

    # factors
    # psixi is precision matrix
    # muxi is mean vector
    xi[i,1:10] ~ dmnorm(muxi,psixi)

    C[i] ~ dcat(p0)
  }

  ##############################################
  # priors
  ##############################################

  p0[1] ~ dunif(0.5,1)
  p0[2] <- 1-p0[1]

  for(j in 1:100){
    psiy[j] ~ dgamma(4,4) }

  for(i in 1:N){
    psiy_i[i] ~ dgamma(4,4) }

  ty2 <- muy0 # Dropped person index, so it's just a single parameter

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
      log_lik0[i,j] <- logdensity.norm(y[i,j],muy[i,j,C[i]],psiy_ij[i,j,C[i]])
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
