
#' Formats long data
#'
#' @param home_dir path to home directory
#' @param data data in wide format as: a data frame, path to folder of imputed .csv files, or mids object
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint" suffix
#' @param time_var (optional) variable name in original dataset demarcating time
#' @param id_var (optional) variable name in original dataset demarcating ID
#' @param missing (optional) indicator for missing data in original dataset
#' @param factor_confounders (optional) list of variable names that are factors (default is numeric)
#' @return formatted long dataset
#' @export
#'
#' @examples
formatLongData <- function(home_dir, data, exposure, exposure_time_pts, outcome, tv_confounders, time_var = NA, id_var = NA, missing = NA, factor_confounders = NULL){

  if (missing(home_dir)){
    stop("Please supply a home directory.", call. = FALSE)
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
  if (missing(exposure_time_pts)){
    stop("Please supply the exposure time points at which you wish to create weights.", call. = FALSE)
  }
  if (missing(tv_confounders)){
    stop("Please supply a list of time-varying confounders.", call. = FALSE)
  }

  time_varying_covariates <- tv_confounders
  options(readr.num_columns = 0)

  if (!dir.exists(home_dir)) {
    stop('Please provide a valid home directory.', call. = FALSE)
  }


  # Reading and formatting LONG dataset
  if (!is.na(time_var)){
    colnames(data)[colnames(data) == time_var] <- "WAVE" # Assigning time variable
  }

  if(!is.na(id_var)){
    colnames(data)[colnames(data) == id_var] <- "ID" # Assigning time variable
  }

  if(!is.na(missing)){
    # data[data == missing] <- NA # Makes NA the missingness indicator
    is.na(data[data == missing]) <- TRUE
  }

  if (which(colnames(data) == "ID") != 1){
    data <- data[,which(colnames(data) == "ID"):ncol(data)]
  }

  # Exposure summary
  exposure_summary <- data %>%
    dplyr::filter(WAVE %in% exposure_time_pts) %>%
    dplyr::group_by(WAVE) %>%
    dplyr::summarize_at(dplyr::vars(all_of(exposure)),
                        list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)

  cat(knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure,
                                                      " Exposure Information"), format = 'pipe'), sep = "\n")
  knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"), format = 'html') %>%
    kableExtra::kable_styling() %>%
    kableExtra::save_kable(file = file.path(home_dir, paste0(exposure, "_exposure_info.html")))

  cat(paste0(exposure, " exposure descriptive statistics have now been saved in the home directory"), "\n")
  cat("\n")


  # Outcome summary
  outcome <- sapply(strsplit(outcome, "\\."), "[",1)
  outcome_summary <- data %>%
    dplyr::group_by(WAVE) %>%
    dplyr::summarize_at(dplyr::vars(all_of(outcome)),
                        list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)

  cat(knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ", outcome, " Information"), format = 'pipe'), sep = "\n")

  knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ", outcome, " Information"), format = 'html') %>%
    kableExtra::kable_styling() %>%
    kableExtra::save_kable(file = file.path(home_dir, paste0(outcome, "_outcome_info.html")))

  cat(paste0(outcome, " outcome descriptive statistics have now been saved in the home directory"), "\n")


  data$ID <- as.factor(data$ID)

  if(!is.null(factor_confounders)){
    if (sum(factor_confounders %in% colnames(data)) < length(factor_confounders)) {
      stop('Please provide factor covariates that correspond to columns in your data when creating the msm object', call. = FALSE)
    }
    # Formatting factor covariates
    data[, factor_confounders] <- lapply(data[, factor_confounders], factor)
    # Formatting numeric covariates
    numeric_vars <- colnames(data)[!colnames(data) %in% c(factor_covariates, id)]
    data[, numeric_vars] <- lapply(data[, numeric_vars], as.numeric)
  }

  as.data.frame(data)
}
