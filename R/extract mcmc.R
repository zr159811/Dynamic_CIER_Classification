
#' Subset mcmc.summary and add dimensions
#'
#' Can parse columns with entries like "variable[1,2]" into multiple columns containing the dimensions
#'
#' @param summary mcmc.summary data.frame (but can be any df)
#' @param variable character to subset variable names in iden.column
#' @param dim vector with dimension names. optional.
#' @param iden.column name of column that contains variable names
#'
#' @return
#' @export
#'
#' @examples
subset_summary <- function(summary, variable, dim=NULL, iden.column="variable"){

  # safe column for convenience
  summary$temporary.iden <- dplyr::pull(summary,iden.column)

  ### subset variable
  summary$temporary.name <- gsub(" ", "",gsub("\\[.*"," ",summary$temporary.iden))
  summary <- summary[summary$temporary.name==variable,]

  # extract dimensions: split by ",", drop before brackets, drop all but numbers
  dims <- gsub("[^0-9.-]", "",
               gsub(".*\\["," ",
                    stringi::stri_split_regex(summary$temporary.iden,pattern = ",",simplify = TRUE)))

  # generic dimnames
  if(is.null(dim)){
    dim <- paste0("dim",1:ncol(dims))
  }

  # add dimension columns to df
  for(dim.cur in dim){
    summary[,dim.cur]<- as.numeric(dims[,which(dim==dim.cur)])
  }

  # drop tmeporary columns
  summary$temporary.name<-NULL
  summary$temporary.iden<-NULL

  return(summary)
}
