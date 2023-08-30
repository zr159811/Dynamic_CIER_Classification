# jags_dyn2: omit variables to reduce RAM requirements
# 2b: fix PCb
# 4: adapt psi prior
model{

  # constants
  for(k in 1:10){
    ly[(k-1)*10+1] <- 1
    ty1[(k-1)*10+1] <- muy0}

  # switching probabilities in non-attention state
  PCb[2] <- 1-p.return
  PCb[1] <- p.return

  for(i in 1:N){
  # i = person

  ##############################################
  # within level
  ##############################################

    for(j in 1:100){
    # j = item
      y[i,j] ~ dnorm(
        ifelse(C[i,j]==1, ty1[j] + ly[j]*xi[i,index[j]], ty2),
        ifelse(C[i,j]==1, psiy[j],psiy_i[i]))
    }


    # factors
    # psixi is precision matrix
    # muxi is mean vector
    xi[i,1:10] ~ dmnorm(muxi,psixi)


    for(j in 2:100){
      # j = time

      # switching probabilities in attention state
      PCa[i,ord[i,j],1] <- 1/(exp(-(b0[1]+b0[2]*vy[i,ord[i,j]] + b0[3]*cd[i,ord[i,j]]))+1)
      PCa[i,ord[i,j],2] <- 1-PCa[i,ord[i,j],1]

      # markov states
      C[i,ord[i,j]] ~ dcat(ifelse(C[i,ord[i,j-1]]==1,PCa[i,ord[i,j],1:2],PCb[1:2]))
    }

    # first state distributed with p0
    C[i,ord[i,1]] ~ dcat(p0)
    PCa[i,ord[i,1],1:2] <- p0

  }

  ##############################################
  # priors
  ##############################################

  for(j in 1:100){
    psiy[j] ~ dgamma(.1,.1) }

  for(i in 1:N){
    psiy_i[i] ~ dgamma(.1,.1) }

  ty2 <- muy0 # Dropped person index, so it's just a single parameter

  for(j in 1:90){
    ly[index3[j]]  ~ dnorm(1,1)T(0,) # T instead of I for clarity
    ty1[index3[j]] ~ dnorm(muy0,1) }

  for(j in 1:10){muxi[j] = 0}
  psixi[1:10,1:10] ~ dwish(R0[1:10,1:10],10)

  b0[1]~dnorm(5,1/25) # increasing standard deviation of intercept.
  b0[2]~dnorm(0,1)
  b0[3]~dnorm(0,1)

  p0[1] ~ dunif(0.5,1) # 50% to 100% of cases begin in an attentive state.
  p0[2] <- 1-p0[1]

  ##############################################
  # Log-likelihood (for model comparison)
  ##############################################
  for(i in 1:N){
    for(j in 1:100){
      log_lik0[i,j] <- logdensity.norm(y[i,j],
                                       ifelse(C[i,j]==1, ty1[j] + ly[j]*xi[i,index[j]], ty2),
                                       ifelse(C[i,j]==1, psiy[j],psiy_i[i]))
    }
    log_lik[i] <- sum(log_lik0[i,])
  }


  ##############################################
  # STD loadings and covariances
  ##############################################

  for(j in 1:100){
    sigmay[j] <- 1/psiy[j] } #zac: right?
  sigmaxi[1:10,1:10] <- inverse(psixi[1:10,1:10])

  for(j in 1:100){
    vx[j] <- ly[j]^2*sigmaxi[index[j],index[j]] + sigmay[j]
    lySt[j] <- ly[j]*sqrt(sigmaxi[index[j],index[j]] / vx[j]) }

  for(q in 1:90){
    Fcor[index2[q,1],index2[q,2]] <- sigmaxi[index2[q,1],index2[q,2]] / sqrt(sigmaxi[index2[q,1],index2[q,1]]*sigmaxi[index2[q,2],index2[q,2]]) }

  # Compute if participant attentive at all time
  for(i in 1:N){
    Call[i]= prod(C[i,1:100]==1) }

  # Compute percentage attentive
  for(i in 1:N){
    Cratio[i]= mean(C[i,1:100]==1) }

}
