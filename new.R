# Functions for the empirical example/ validation
# (Roman, Schmidt, Miller, and Brandt, 202X)

# For demographic data
cnrmturk::demographics

# For response to bogus items
cnrmturk::bs
# Bogus responses in long format
cnrmturk::bs.long

# For responses to items NOT in presentation order
cnrmturk::ymat

# For raw presentation order of items
cnrmturk::raw_ord

# For temporal presentation order of items
cnrmturk::tord

# For response times
cnrmturk::time

# For experimental condition
# Condition 1 = Control (survey instructions as usual)
# Condition 2 = Second Control, but instructed to SERIOUSLY pay attention because
#               we have an algorithm that will know if they dont
# Condition 3 = Randomly told to Finish as fast as possible and not pay attention
# Condition 4 = Same instructions as 2, then randomly told to not pay attention
#               but make it look like you did.
cnrmturk::cond

# Person fit function overtime (Not in presentation order)
cnrmturk::cd_mat_full


# A more convenient way of working with the objects is the pull data function
# You can choose conditions, and randomly sample cases with n or standardize.
d <- cnrmturk::pull_data(choice.cond = c(2,4),
                         n = 50,
                         standardize = F)

# Pulls all data in consistant order
d <- cnrmturk::pull_data(choice.cond = c(1,2,3,4),standardize = F)

# This contains all necessary objects for the analysis

# Response matrix NOT presentation order
ymat <- d$ymat

# Presentation order
iord <- d$iord

# Response time
Rtime <- d$time

# Person fit function in item order
CDfun <- d$cd

# Instructions
d$instr

# Condition
d$cond

# Raw order for verifying re-ordering by iord
d$raw_ord

# To transform to presentation order:
holder <- list()
for(i in 1:nrow(ymat)){
  holder[[i]] <-  as.numeric(ymat[i,iord[i,]])
}
ymati <- do.call(rbind,holder)

# This is the response matrix in the order each participant was presented the
# items.

ymati






