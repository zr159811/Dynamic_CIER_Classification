
#' Generates post-proc report for jags model results
#' @param path location of jags output object as defined
#' @return Returns pdf document report
#' @export
JagsReport<- function(path){
  require(knitr)
  require(stringr)
  require(kableExtra)
  require(ggplot2)
  # https://stackoverflow.com/questions/30377213/how-to-include-rmarkdown-file-in-r-package/30377598
  # system.file("rmd", "file.Rmd", package = "packagename")
  template <- system.file("rmd", "JagsReportTemplate.Rmd", package = "cnrmturk")

  MD <- read_lines(template)
  MD[21] <- str_replace(string = MD[21],
                        pattern = "JagsOut",
                        replacement = path)

  writeLines(MD,"JagsReport.Rmd")
  knit(input = "JagsReport.Rmd")
  pandoc(input = "JagsReport.md",
         format = "pdf")

  file.rename("JagsReport.pdf",paste0("JagsReport",format(Sys.time(),"%d_%m_%y"),".pdf"))

  unlink("JagsReportTemplate.md")
  unlink("JagsReport.Rmd")
  unlink("JagsReport.md")
  unlink("figure",
         recursive = T)
}


#' Generates post-proc report for jags model results with multiple conditions
#' @param path location of jags output object as defined
#' @return Returns pdf document report
#' @export
JagsReportMultiCon<- function(path){
  require(knitr)
  require(stringr)
  require(kableExtra)
  require(ggplot2)
  # https://stackoverflow.com/questions/30377213/how-to-include-rmarkdown-file-in-r-package/30377598
  # system.file("rmd", "file.Rmd", package = "packagename")
  template <- system.file("rmd", "JagsReportTemplateMultiCon.Rmd", package = "cnrmturk")

  MD <- read_lines(template)
  MD[21] <- str_replace(string = MD[21],
                        pattern = "JagsOut",
                        replacement = path)

  writeLines(MD,"JagsReport.Rmd")
  knit(input = "JagsReport.Rmd")
  pandoc(input = "JagsReport.md",
         format = "pdf")

 file.rename("JagsReport.pdf",paste0("JagsReportMC",format(Sys.time(),"%d_%m_%y"),".pdf"))

  unlink("JagsReportTemplate.md")
  unlink("JagsReport.Rmd")
  unlink("JagsReport.md")
  unlink("figure",
         recursive = T)
}
