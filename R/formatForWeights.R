#' Format data for creating balancing weights
#'
#' This code formats the imputed datasets for calculating balancing weights
#' @param msm_object msm object that contains all relevant user inputs
#' @param data original dataset
#' @param imputed_datasets output from imputeData
#' @return wide_long_datasets
#' @export
#' @importFrom readr read_csv
#' @importFrom stats reshape
#' @importFrom mice complete
#' @importFrom dplyr select
#' @importFrom purrr map
#' @seealso [formatDataStruct()], [imputeData()]
#' @examples formatForWeights(msm_object, data, imputed_datasets)
#'
formatForWeights <- function(msm_object, data, imputed_datasets) {
  ID <- msm_object$ID
  home_dir <- msm_object$home_dir
  m <- msm_object$m
  time_pts <- msm_object$time_pts
  time_var_exclude <- msm_object$time_var_exclude
  time_varying_covariates <- msm_object$time_varying_variables
  exposure <- msm_object$exposure
  outcome <- msm_object$outcome

  # Load necessary packages
  library(readr)
  library(stats)
  library(mice)
  library(dplyr)
  library(purrr)

  options(readr.num_columns = 0)

  # Cycles through imputed datasets and puts them in wide dataset
  wide_long_datasets <- map(1:m, function(k) {
    imp <- as.data.frame(mice::complete(imputed_datasets, k))

    # Finds all time-varying variables in wide format
    time_varying_wide <- apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep = "", collapse = ".")
    time_varying_wide <- sort(time_varying_wide)
    time_varying_wide <- c(ID, time_varying_wide)
    time_varying_wide <- time_varying_wide[!time_varying_wide %in% time_var_exclude] # Removes time-varying time pts that should not be there

    # Creates wide dataset
    imp_wide <- suppressWarnings(stats::reshape(data = imp,
                                                idvar = ID,
                                                v.names = time_varying_covariates,
                                                timevar = "WAVE",
                                                times = c(time_pts),
                                                direction = "wide"))
    imp_wide <- imp_wide[, !colnames(imp_wide) %in% time_var_exclude] # Confirm: only include what should be there

    return(imp_wide)
  })
  names(wide_long_datasets) <- 1:m

  return(wide_long_datasets)
}
