#' Imputes dataset so there is no missing at each time point using parallel
#' processing to speed up
#'
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
#' @importFrom missMethods delete_MAR_1_to_x
#' @seealso {[mice::mice()],
#'   <https://cran.r-project.org/web/packages/mice/index.html>}
#' @param data data in wide format
#' @param m (optional) integer number of imputed datasets (default is 5)
#' @param method (optional) character string of imputation method from mice()
#'   (default is random forest "rf")
#' @param home_dir path to home directory
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param para_proc (optional) TRUE/FALSE whether to do parallel processing
#'   using multiple cores to speed up process (default = TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param read_imps_from_file (optional) "yes" or "no" indicatorto read in weights
#'   that have been previously run and saved locally (default is "no")
#' @return mice object of m imputed datasets
#' @examples
#' test <- data.frame(ID = 1:50,
#'                    A.1 = rnorm(n = 50),
#'                    A.2 = rnorm(n = 50),
#'                    A.3 = rnorm(n = 50),
#'                    B.1 = rnorm(n = 50),
#'                    B.2 = rnorm(n = 50),
#'                    B.3 = rnorm(n = 50),
#'                    C = rnorm(n = 50),
#'                    D.3 = rnorm(n = 50))
#' test[, c("A.1", "A.2", "A.3")] <- lapply(test[, c("A.1", "A.2", "A.3")], as.numeric)
#'
#' test_miss <- missMethods::delete_MAR_1_to_x(as.data.frame(test), p = 0.20,
#'                                             cols_mis = c("A.1", "B.2", "C"),
#'                                             cols_ctrl = c("B.1", "B.1", "B.1"), 3)
#' test_i <- imputeData(data = test_miss,
#'                      m = 3,
#'                      method = "rf",
#'                      exposure = "A",
#'                      outcome = "D.3",
#'                      para_proc = TRUE,
#'                      read_imps_from_file = "no",
#'                      save.out = FALSE)


imputeData <- function(data, m = 5, method = "rf", home_dir = NA, exposure, outcome, para_proc = TRUE,
                       read_imps_from_file = "no", save.out = TRUE) {

  if (save.out | read_imps_from_file == "yes"){
    if (missing(home_dir)){
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if (!dir.exists(home_dir)) {
      stop('Please provide a valid home directory.', call. = FALSE)
    }
  }
  if (missing(data)){
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  if (missing(exposure)){
    stop("Please supply a single exposure.", call. = FALSE)
  }
  if (missing(outcome)){
    stop("Please supply a single outcome.", call. = FALSE)
  }


  if(!is.character(method)){
    stop("Please provide as a character a valid imputation method abbreviation.", call. = FALSE)
  }
  if(!is.numeric(m)){
    stop("Please provide an integer value number of imputations.", call. = FALSE)
  }

  if (save.out | read_imps_from_file == "yes"){
    imp_dir <- file.path(home_dir, "imputations")
    if (!dir.exists(imp_dir)) {
      dir.create(imp_dir)
    }
  }


  if (read_imps_from_file == "yes") {
    imputed_datasets <- list()

    if (!file.exists(glue::glue("{home_dir}/imputations/{exposure}-{outcome}_all_imp.rds"))) {
      stop("Imputations have not been created and saved locally. Please set 'read_imps_from_file' == 'no' and re-run.", call. = FALSE)
    }

    imp <- readRDS(glue::glue("{home_dir}/imputations/{exposure}-{outcome}_all_imp.rds"))
    imputed_datasets <- imp

    cat("\n")
    cat(glue::glue("Reading in {imputed_datasets$m} imputations from the local folder."))
    cat("\n")
    return(imputed_datasets)

  }
  else {

    if (sum(duplicated(data$"ID")) > 0){
      stop("Please provide a wide dataset with a single row per ID.", call. = FALSE)
    }

    imp_method <- method
    data_to_impute <- tibble::tibble(data)

    cat(glue::glue("Creating {m} imputed datasets using the {imp_method} imputation method in mice. This may take some time to run."))
    cat("\n")

    if (para_proc){
      # Configure parallelization
      nCores <- min(parallel::detectCores(), 8)
      options(mc.cores = nCores)
      options(cores = nCores)
      doParallel::registerDoParallel(cores = nCores)

      cat("### Using", foreach::getDoParWorkers(), "cores\n")
      cat("### Using", foreach::getDoParName(), "as the backend\n")

      # Conducts imputations using parallelized execution cycling through m
      imputed_datasets <- foreach::foreach(i = seq_len(m), .combine = mice::ibind) %dorng% {
        cat("### Started iteration", i, "\n")
        miceout <- mice::mice(data_to_impute, m = 1, method = imp_method, maxit = 0, #change maxit to default 5 after testing!!!
                              print = F)
        cat("### Completed iteration", i, "\n")
        miceout
      }
    }
    else{
      imputed_datasets <- mice::mice(data_to_impute, m = m, method = imp_method, maxit = 0, #change maxit to default 5 after testing!!!
                                     print = F)
    }

    if(save.out){
      saveRDS(imputed_datasets, glue::glue("{home_dir}/imputations/{exposure}-{outcome}_all_imp.rds"))
    }

    # Print warnings
    cat("USER ALERT: Please view any logged events from the imputation below:", "\n")
    cat(knitr::kable(imputed_datasets$loggedEvents, caption = "Logged Events from mice", format = 'pipe'), sep = "\n")
    cat("\n")

    if(save.out){
      # Save out individual imputed datasets
      for (k in seq_len(m)) {
        write.csv(mice::complete(imputed_datasets, k),
                  file = glue::glue("{home_dir}/imputations/{exposure}-{outcome}_imp{k}.csv"))
      }
      cat("See the 'imputations/' folder for a .csv file of each imputed dataset and an .rds file of all imputed datasets", "\n")
    }

    imputed_datasets
  }
}
