
#' Plots traceplots of worst mixing parameters
#'
#' @param mcmc mcmc output
#' @param summary summary output with Rhat column
#' @param names names of variables to consider, if NULL considers all
#' @param regex regex expression to find names
#' @param n number of parameters to plot
#'
#' @return ggplot
#' @export
plot_worst_trace <- function(mcmc=NULL,
                             summary=NULL,
                             names = NULL,
                             regex = NULL,
                             n = 12){

  # if not given, load from global environment
  if(is.null(mcmc)){
    mcmc <- res.stan
  }

  # ensure mcmc.list format for stanfit objects
  if("stanfit" %in% class(mcmc)){
    mcmc <- rstan::As.mcmc.list(mcmc)
  }

  # if not given laod from global environment
  if(is.null(summary)){
    summary <- mcmc.summary
  }

  # safe names that exist
  varnames.exist <- unique(summary$name)

  # choose names
  if(is.null(regex)){

    # use all names if none given
    if(is.null(names)){
      names <- unique(summary$name)}

    summary <- summary%>%
      filter(name %in% names)
  } else{

    summary <- summary%>%
      filter(grepl(pattern = regex,x = name))

    if(!is.null(names))
      warning('name paramter is ignored as regex is given.')
  }

  if(nrow(summary)==0){
    stop('No variable names selected with choice of regex/names. The following exist:',paste(varnames.exist,collapse = ', '))
  }

  # extract varnames with highest Rhat
  varnames <- summary%>%
    slice_max(Rhat,n = n)%>%
    pull(varname)


  # plot traceplot
  if("mcmc.list" %in% class(mcmc))
    p <- mcmcplots::traplot(mcmc,varnames)
  else if("draws_array" %in% class(mcmc))
   p <- bayesplot::mcmc_trace(mcmc,varnames)
  else
    stop('mcmc does not have right class')

  return(p)
}





#' Plot shifted states plot with treatment timing at 0
#'
#' @param df data.frame that contains
#' "treatment.timing",
#' "id" (person identifier),
#' ".value" (value to plot, normally posterior mean),
#' "t" (item timing)
#' @param summary should summary stats be plotted
#'
#' @return
#' @export
#'
#' @examples
plot_states_shifted <- function(df, summary=TRUE){

  if(any(!c("treatment.timing","id",".value","t") %in% colnames(df)))
    stop('Missing info in data.frame')

  if(summary)
    return(df %>%
             mutate(t.transform=t-treatment.timing)%>%
             group_by(t.transform)%>%
             mutate(
               mean = mean(.value),
               qlow = quantile(.value,probs=.1),
               qhigh = quantile(.value,probs=.9)
             )%>%
             ggplot(aes(x=t.transform,y=.value,group=id))+
             geom_line(alpha=.1)+
             geom_line(aes(y=mean),color="blue")+
             geom_line(aes(y=qlow),color="darkblue")+
             geom_line(aes(y=qhigh),color="darkblue")+
             geom_vline(xintercept = 0,color="red"))
  else
    return(df %>%
             mutate(t.transform=t-treatment.timing)%>%
             group_by(t.transform)%>%
             ggplot(aes(x=t.transform,y=.value,group=id))+
             geom_line(alpha=.1)+
             geom_vline(xintercept = 0,color="red"))
}
