#' Imputes dataset so there is no missing at each time point using parallel processing to speed up
#'
#' Creates m imputed datasets from original datasets using mice::mice
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
#'
imputeData <- function(data, m, method, home_dir, exposure, outcome, tv_confounders, ti_confounders, read_imps_from_file="no") {
  if (!dir.exists(home_dir)) {
    stop('Please provide a valid home directory.')
  }


  if (read_imps_from_file == "yes") {
    imputed_datasets <- list()

    if (!file.exists(glue::glue("{home_dir}/imputations/{exposure}-{outcome}_all_imp.rds"))) {
      stop("Imputations have not been created and saved locally. Please set 'read_imps_from_file' == 'no' and re-run.")
    }

    imp <- readRDS(glue::glue("{home_dir}/imputations/{exposure}-{outcome}_all_imp.rds"))
    imputed_datasets <- imp

    cat("\n")
    cat(glue::glue("Reading in {imputed_datasets$m} imputations from the local folder."))
    cat("\n")
    return(imputed_datasets)

  } else {

    if (sum(duplicated(data$"ID")) > 0){
      stop("Please provide a wide dataset with a single row per ID.")
    }

    # library(mice)
    # library(doParallel)
    # library(doRNG)
    # library(purrr)
    # library(tibble)

    imp_method <- method
    time_varying_covariates <- tv_confounders
    time_invar <- ti_confounders

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

    # Conducts imputations using parallelized execution cycling through m
    imputed_datasets <- foreach(i = seq_len(m), .combine = mice::ibind) %dorng% {
      cat("### Started iteration", i, "\n")
      miceout <- mice::mice(data_to_impute_full, m = 1, method = imp_method, maxit = 0,
                            print = F)
      cat("### Completed iteration", i, "\n")
      miceout
    }

    saveRDS(imputed_datasets, glue::glue("{home_dir}/imputations/{exposure}-{outcome}_all_imp.rds"))

    # Print warnings
    cat("USER ALERT: Please view any logged events from the imputation below:", "\n")
    cat(knitr::kable(imputed_datasets$loggedEvents, caption = "Logged Events from mice", format = 'pipe'), sep = "\n")
    cat("\n")

    # Save out individual imputed datasets
    for (k in 1:m) {
      write.csv(mice::complete(imputed_datasets, k),
                file = glue::glue("{home_dir}/imputations/{exposure}-{outcome}_imp{k}.csv"))
    }
    cat("See the 'imputations/' folder for a .csv file of each imputed dataset and an .rds file of all imputed datasets", "\n")

    imputed_datasets
  }
}
