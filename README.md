## Data, models, and functions for the empirical example validation (Roman, Schmidt, Miller, and Brandt, In review)

Install:

`# install.packages("remotes")
remotes::install_github("zr159811/Dynamic_CIER_Classification")`

For demographic data

`cnrmturk::demographics`

For response to bogus items

`cnrmturk::bs`

Bogus responses in long format

`cnrmturk::bs.long`

For responses to items NOT in presentation order

`cnrmturk::ymat`

For raw presentation order of items

`cnrmturk::raw_ord`

For temporal presentation order of items

`cnrmturk::tord`

For response times

`cnrmturk::time`

For experimental condition

1.  Condition 1 = Control (survey instructions as usual)

2.  Condition 2 = Second Control, but instructed to SERIOUSLY pay
    attention because we have an algorithm that will know if they dont

3.  Condition 3 = Randomly told to Finish as fast as possible and not
    pay attention

4.  Condition 4 = Same instructions as 2, then randomly told to not pay
    attention but make it look like you did.

`cnrmturk::cond`

Person fit function overtime (Not in presentation order)

`cnrmturk::cd_mat_full`

A more convenient way of working with the objects is the pull data
function

You can choose conditions, and randomly sample cases with n or
standardize.

`d <- cnrmturk::pull_data(choice.cond = c(2,4), n = 50, standardize = F)`

Pulls all data in consistent order

`d <- cnrmturk::pull_data(choice.cond = c(1,2,3,4),standardize = F)`

This contains all necessary objects for the analysis

Response matrix NOT presentation order

`ymat <- d$ymat`

Presentation order

`iord <- d$iord`

Response time

`Rtime <- d$time`

Person fit function in item order

`CDfun <- d$cd`

Instructions

`d$instr`

Condition

`d$cond`

Raw order for verifying re-ordering by iord

`d$raw_ord`

To transform to presentation order:

`holder <- list()`

`for(i in 1:nrow(ymat)){holder[[i]] <- as.numeric(ymat[i,iord[i,]]) }`

`ymati <- do.call(rbind,holder)`

This is the response matrix in the order each participant was presented
the

items.

`ymati`

# Modeling

Regular CFA model:

`analysis/``Jags_CFA.R`

Static latent classification + CFA model:

`analysis/``jags_STA.R`

Dynamic latent classification + CFA model:

`analysis/jags_DYN.R`

To run these models the following pre-proc is necessary:

``` r
# Factor loading design matricies
index <- c(rep(1,10),
           rep(2,10),
           rep(3,10),
           rep(4,10),
           rep(5,10),
           rep(6,10),
           rep(7,10),
           rep(8,10),
           rep(9,10),
           rep(10,10))

xx <- expand.grid(1:10,1:10)
index2 <- xx[ xx[,1] != xx[,2],]

# Omits marker variables (first indicator of each factor)
index3 <- seq(1,100,1)[-c(1,11,21,31,41,51,61,71,81,91)]

# Creat input
input <- list(y = d$ymat, # raw ymatrix
              ord = d$iord, # Item order matrix
              R0 = diag(10), # number of factors = 10
              N = dim(d$ymat)[1], # N cases
              muy0 = 3.5, # Centerpoint of possible responses.
              cd = d$cd,
              vy = makevy2(ymat = d$ymat, # Creates variance function
                           iord = d$iord,
                           tord = d$tord,
                           startLag = 6), # Lagg of starting variance. 
              index = index,   # Design matricies
              index2 = index2,
              index3 = index3)


params <- c("lySt", # Standardized factor loadings
            "ly",   # Raw factor loadings
            "Fcor", # Factor correlations
            "sigmaxi", # residual variance of latent factors
            "sigmay", # Residual variance of items
            "ty1", # Item intercepts
            "b0", # Intercept probability of markov model
            "C", # Class assignment (CNR or Not CNR)
            "PC1", # Probability of CNR
            "xi", # Factor Scores
            "log_lik") # Log Liklihood

# Fit Dynamic Model
fit1 <-
  jags(
    data=input,
    parameters.to.save=params,
    n.iter=1000, n.chains=1,n.thin=1,n.burnin=500,
    model.file="analysis/jags_DYN.R")

# Creat input
input <- list(y = d$ymat, # raw ymatrix
              R0 = diag(10), # number of factors = 10
              N = dim(d$ymat)[1], # N cases
              muy0 = 3.5, # Centerpoint of possible responses.
              cd = d$cd[,100],
              vy = makevy2(ymat = d$ymat, # Creates variance function
                           iord = d$iord,
                           tord = d$tord,
                           startLag = 6)[,100], 
              index = index,   # Design matricies
              index2 = index2,
              index3 = index3)

params <- c("lySt", # Standardized factor loadings
            "ly",   # Raw factor loadings
            "Fcor", # Factor correlations
            "sigmaxi", # residual variance of latent factors
            "sigmay", # Residual variance of items
            "ty1", # Item intercepts
            "b0", # Intercept probability of markov model
            "C", # Class assignment (CNR or Not CNR) Item level
            "xi", # Factor Scores
            "log_lik") # Log Liklihood

# Fit static Model
fit1 <-
  jags(
    data=input,
    parameters.to.save=params,
    n.iter=1000, n.chains=1,n.thin=1,n.burnin=500,
    model.file="analysis/jags_STA.R")

# Creat input
input <- list(y = d$ymat, # raw ymatrix
              R0 = diag(10), # number of factors = 10
              N = dim(d$ymat)[1], # N cases
              muy0 = 3.5, # Centerpoint of possible responses.
              cd = d$cd[,100],
              vy = makevy2(ymat = d$ymat, # Creates variance function
                           iord = d$iord,
                           tord = d$tord,
                           startLag = 6)[,100], 
              index = index,   # Design matricies
              index2 = index2,
              index3 = index3)

params <- c("lySt", # Standardized factor loadings
            "ly",   # Raw factor loadings
            "Fcor", # Factor correlations
            "sigmaxi", # residual variance of latent factors
            "sigmay", # Residual variance of items
            "ty1", # Item intercepts
            "b0", # Intercept probability of markov model
            "C", # Class assignment (CNR or Not CNR) person level
            "xi", # Factor Scores
            "log_lik") # Log Liklihood

# Fit just CFA Model
fit1 <-
  jags(
    data=input,
    parameters.to.save=params,
    n.iter=1000, n.chains=1,n.thin=1,n.burnin=500,
    model.file="analysis/jags_CFA.R")
```
