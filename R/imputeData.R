#' Imputes dataset so there is no missing at each time point using parallel processing to speed up
#'
#' Creates m imputed datasets from original datasets using mice::mice
#' @param msm_object msm object that contains all relevant user inputs
#' @param data_to_impute output from dataToImpute
#' @return imputed_datasets imputation results
#' @export
#' @importFrom mice mice
#' @importFrom mice ibind
#' @importFrom mice complete
#' @importFrom tibble tibble
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom tidyr complete
#' @importFrom knitr kable
#' @importFrom parallel detectCores
#' @importFrom doRNG %dorng%
#' @importFrom purrr map_dfr
#' @importFrom foreach getDoParWorkers
#' @importFrom foreach getDoParName
#' @examples imputeData(msm_object, data_to_impute, read_imps_from_file = "no")
#'
imputeData <- function(data, m, method, home_dir, exposure, outcome, tv_confounders, ti_confounders, read_imps_from_file="no") {
  # home_dir <- msm_object$home_dir
  ID <- "ID"
  # continuous_variables <- msm_object$continuous_variables
  # m <- msm_object$m
  imp_method <- method
  time_varying_covariates <- tv_confounders
  time_invar <- ti_confounders
  # factor_covariates <- msm_object$factor_covariates
  # exposure <- msm_object$exposure
  # outcome <- msm_object$outcome

  if (read_imps_from_file == "yes") {
    imputed_datasets <- list()
    # Error checking
    if (!file.exists(glue::glue("{home_dir}imputations/{exposure}-{outcome}_all_imp.rds"))) {
      stop("Imputations have not been created and saved locally. Please set 'read_imps_from_file' == 'no' and re-run.")
    }
    imp <- readRDS(glue::glue("{home_dir}imputations/{exposure}-{outcome}_all_imp.rds"))
    imputed_datasets <- imp
    cat("\n")
    cat(glue::glue("Reading in {imputed_datasets$m} imputations from the local folder."))
    cat("\n")
    return(imputed_datasets)
  } else {
    library(mice)
    library(doParallel)
    library(doRNG)
    library(purrr)
    library(tibble)

    # Set seed for reproducibility
    set.seed(123)

    cat(glue::glue("Creating {m} imputed datasets using the {imp_method} imputation method in mice. This may take some time to run."))
    cat("\n")

    # Configure parallelization
    nCores <- min(parallel::detectCores(), 8)
    options(mc.cores = nCores)
    options(cores = nCores)
    doParallel::registerDoParallel(cores = nCores)
    cat("### Using", foreach::getDoParWorkers(), "cores\n")
    cat("### Using", foreach::getDoParName(), "as the backend\n")

    data_to_impute <- tibble(data_to_impute)
    # data_to_impute[, ID] <- as.factor(data_to_impute[, ID])

    # Adding missing longitudinal time points in long format to create a full dataset
    # data_to_impute_full <- data_to_impute %>%
    #   complete(.data[[ID]], WAVE)

    # Filling in all person-level variables (i.e., those not designated as time-varying)
    # time_invar <- colnames(data_to_impute)[!colnames(data_to_impute) %in% time_varying_covariates]
    # time_invar <- time_invar[!time_invar %in% c(ID, "WAVE")]
    person_info <- data_to_impute[!duplicated(data_to_impute$ID),] %>%
      dplyr::select(c(all_of(ID), all_of(time_invar)))
    for (x in 1:nrow(person_info)) {
      id <- as.character(person_info[x, ID])
      data_to_impute_full[data_to_impute_full[colnames(data_to_impute) == ID] == id,
                          colnames(data_to_impute_full) %in% time_invar] <- person_info[x, time_invar]
    }

    # Making variable types
    # data_to_impute_full[, colnames(data_to_impute_full) %in% factor_covariates] <- lapply(data_to_impute_full[, colnames(data_to_impute_full) %in% factor_covariates], factor)

    # Conducts imputations using parallelized execution cycling through m
    imputed_datasets <- map_dfr(1:m, ~{
      cat(glue::glue("### Started iteration {i}\n"))
      miceout <- mice::mice(data_to_impute_full, m = 1, method = imp_method, maxit = 5, print = FALSE)
      cat(glue::glue("### Completed iteration {i}\n"))
      miceout
    })

    saveRDS(imputed_datasets, glue::glue("{home_dir}imputations/{exposure}-{outcome}_all_imp.rds"))

    # Print warnings
    cat("USER ALERT: Please view any logged events from the imputation below:", "\n")
    cat(knitr::kable(imputed_datasets$loggedEvents, caption = "Logged Events from mice", format = 'pipe'), sep = "\n")
    cat("\n")

    # Save out individual imputed datasets
    for (k in 1:m) {
      write.csv(mice::complete(imputed_datasets, k), file = glue::glue("{home_dir}imputations/{exposure}-{outcome}_imp{k}.csv"))
    }
    cat("See the 'imputations/' folder for a .csv file of each imputed dataset and an .rds file of all imputed datasets", "\n")

    imputed_datasets
  }
}
