
#' Creates mplus object
#' Gather information for creating Mplus input template files
#' @param data_file location of wide/long .csv file with outcomes, exposures, covariates, and weights for modeling
#' @param covariates data frame specifying the names of outputs as well as any technical covariates for modeling
#' @return object

#' @examples mplusObject(data_file, covariates, reference, hi_cutoff=0.75, lo_cutoff=0.25)
mplusObject <- function(data_file, additional_covariates){

  charOrNull <- function(x) {
    is.character(x) || is.null(x)
  }


  stopifnot(charOrNull(data_file))

  object<-list(data_file=data_file, additional_covariates=additional_covariates)

  class(object)<-c("list","mplusObject")

  return(object)



}
