#' logistic function
#'
#' @param x numeric
#'
#' @return numeric
#' @export
logistic <- function(x) {
  return(1/(1+exp(-x)))
  }
