

#' Pull package data of online survey
#'
#' Returns relevant data in list. Contains the following elements:
#' \itemize{
#'   \item{ymat[person,item] with item values}
#'   \item{cond[person] with treatment condition}
#'   \item{iord[person,time] with item index}
#'   \item{tord[person,item] with time shown of item}
#'   \item{treatment.timing[person] with timing of treatment}
#'   \item{cd[person,item] with timing of treatment}
#'   \item{instr[person] with response to instruction question}
#' }
#'
#' @param bs use bs items (standard FALSE)
#' @param n number of persons to be sampled
#' @param choice.cond treatment condition to include (standard only 4)
#' @param standardize standardize item values (standard TRUE)
#' @param treatment.timing.min include only
#' @param treatment.timing.max
#'
#' @return list with data entries
#'
#' @examples
#' dat <- pull_data()
#' dim(dat$ymat)
#' all(colMeans(dat$ymat)<1e-10)
#' ls(dat)
#'
#' dim(dat$cd)
#'
#' ### check order
#' # show first person
#' dat$raw_ord[1]
#' # check first item shown (t=1)
#' which(dat$tord[1,]==1)
#' dat$iord[1,1]
#' # check last item shown (t=100)
#' colnames(ymat)[dat$iord[1,100]]
#'
#' #check place of first item (i=1)
#' dat$tord[1,1]
#'
#' @export
pull_data <- function(choice.cond = 4,
                      bs = FALSE,
                      n = NULL,
                      standardize = TRUE,
                      treatment.timing.min = -Inf,
                      treatment.timing.max = Inf,
                      choice.id = NULL
){

  # load data from package
  ymat <- cnrmturk::ymat
  tord <- cnrmturk::tord
  treatment.timing <- cnrmturk::treatment.timing
  time <- cnrmturk::time
  cond <- cnrmturk::cond
  raw_ord <- cnrmturk::raw_ord
  cd <- cnrmturk::cd_mat_full
  instr <- cnrmturk::instr

  # choose sample of persons based on choice.id
  if(!is.null(choice.id)){
    message("Choose given ids, ignore all other parameters")
    sample <- rownames(ymat) %in% choice.id
  } else {
    # choose sample of persons based on
    # treatment timing,treatment condition, and sample size n

    # safe indices of right conditions
    timing.right <-
      (treatment.timing >= treatment.timing.min &
         treatment.timing <= treatment.timing.max)
    timing.right[is.na(timing.right)]<-TRUE
    sample <- which(cond %in% choice.cond &timing.right)

    # subsample
    if(!is.null(n)){
      sample <- sample(sample,size = n)
    }
  }

  # if bullshit items
  if(bs){
    stop('not implemented')
  }

  ymat <- ymat[sample,]
  tord <- tord[sample,]
  treatment.timing <- treatment.timing[sample]
  time <- time[sample,]
  cond <- cond[sample]
  raw_ord <- raw_ord[sample]
  cd <- cd[sample,]
  instr <- instr[sample]

  if(standardize){
    if(is.logical(standardize)){
      message("Standardize mean and sd by person")
      for(j in 1:ncol(ymat))
        ymat[,j]<-scale(ymat[,j])
    } else if(standardize=="byitem"){
      message("Standardize mean and sd by item")
      for(i in 1:nrow(ymat))
        ymat[i,]<-scale(ymat[i,])
    }
  }

  ### gen iord
  # iord[j,t] returns item i at time t
  # tord[j,i] returns time t of item i
  iord <- matrix(NA, nrow(tord),ncol(tord))
  rownames(iord)<-rownames(tord)
  for(t in 1:ncol(iord)){
    for(j in 1:nrow(iord)){
      iord[j,t]<- which(tord[j,]==t)
    }
  }

  return(list(
    ymat = ymat,
    iord = iord,
    tord = tord,
    time = time,
    cond = cond,
    treatment.timing = treatment.timing,
    raw_ord = raw_ord,
    cd = cd,
    n_item = ncol(ymat),
    n_person = nrow(ymat),
    instr = instr
  ))

}


#' generate data.frame with survey data
#'
#' @param dat.list list of data, can be generated with pull_data, if no dat.list ist provided, the full data from the package is used
#'
#' @return data.frame
#' @export
gen_df <- function(dat.list=NULL){

  # if no data given, use from package
  if(is.null(dat.list)){
    dat.list <- pull_data(choice.cond = 1:4)
  }

  ### time
  time.df <- dat.list$time %>%
    as.data.frame()%>%
    rownames_to_column('id')%>%
    mutate(j=1:n())%>%
    pivot_longer(-c(id,j),names_to = 'iname',values_to="time")%>%
    group_by(iname)%>%
    mutate(timerank = rank(time)/n())%>%
    ungroup()

  ### item
  item.df <- dat.list$ymat %>%
    as.data.frame()%>%
    rownames_to_column('id')%>%
    mutate(j=1:n())%>%
    pivot_longer(-c(id,j),names_to = 'iname',values_to="value")%>%
    mutate(value = as.numeric(value),
           i = rep(1:n_distinct(iname),n_distinct(id)))


  ### order
  ord.df <- dat.list$tord%>%
    as.data.frame()%>%
    rename_with(function(x) gsub('V','',x))%>%
    rownames_to_column('id')%>%
    pivot_longer(-c(id),names_to = 'iname',values_to = 't')

  ### treatment.timing
  treatment.timing.df <- tibble(
    treatment.timing = dat.list$treatment.timing,
    id = names(dat.list$treatment.timing))


  ### instructions
  instr.df <- tibble(
    instr = dat.list$instr,
    id=names(dat.list$instr)
  )

  # condition
  cond.df <-
    tibble(
      cond=dat.list$cond,
      id=names(dat.list$cond)
    )

  ### joined info
  df.joined <-
    time.df %>%
    left_join(item.df)%>%
    left_join(ord.df)%>%
    left_join(treatment.timing.df)%>%
    left_join(instr.df)%>%
    left_join(cond.df)%>%
    mutate(across(time, as.numeric))


  df.joined <- df.joined%>%
    mutate(value.org = ifelse(iname %in% reverse.index, 7-value, value))

  return(df.joined)
}
