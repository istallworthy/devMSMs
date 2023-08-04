#' Formats dataset and creates required directory
#'
#' Creates the required directories, assesses the data structure, and makes it uniform for future functions and returns data
#' This function requires a clean dataset in long format that contains columns for ID, time, exposure, outcome, and potential covariate confounds
#'
#' @param msm_object msm object that contains all relevant user inputs
#' @return data formatted dataset
#' @export
#' @importFrom dplyr filter group_by summarize mutate n
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling save_kable
#' @importFrom gtools permutations quantile
#' @importFrom readr read_csv
#' @examples formatDataStruct(object)

formatDataStruct <- function(msm_object){
  data_path <- msm_object$data_path
  home_dir <- msm_object$home_dir
  missing <- msm_object$missing
  time_var <- msm_object$time_var
  time_varying_covariates <- msm_object$time_varying_variables
  factor_covariates <- msm_object$factor_covariates
  id <- msm_object$ID
  exposure <- msm_object$exposure
  outcome <- msm_object$outcome
  exposure_time_pts <- msm_object$exposure_time_pts
  exposure_epochs <- msm_object$exposure_epochs
  outcome_time_pt <- msm_object$outcome_time_pt
  time_varying_covariates <- msm_object$time_varying_variables
  time_pts <- msm_object$time_pts
  time_var_exclude <- msm_object$time_var_exclude
  hi_cutoff <- msm_object$hi_cutoff
  lo_cutoff <- msm_object$lo_cutoff

  cat("Formatting data")

  options(readr.num_columns = 0)

  # Error checking
  if (!file.exists(data_path)) {
    stop('Please provide a valid directory for your data in data_path when creating the msm object')
  }
  if (!dir.exists(home_dir)) {
    stop('Please provide a valid home directory in home_dir when creating the msm object')
  }

  # Creating all necessary directories within the home directory
  required_dirs <- c(
    "imputations", "original weights", "original weights/values", "original weights/histograms",
    "final weights", "final weights/values", "final weights/histograms",
    "forms", "pre balance", "pre balance/plots",
    "balance", "balance/plots", "balance/comparison values",
    "msms", "msms/estimated means", "msms/contrasts",
    "results figures"
  )

  for (dir in required_dirs) {
    if (!dir.exists(file.path(home_dir, dir))) {
      dir.create(file.path(home_dir, dir))
    }
  }

  cat("\n")
  # Reading and formatting LONG dataset
  data <- suppressWarnings(readr::read_csv(data_path, show_col_types = FALSE))

  colnames(data)[colnames(data) == time_var] <- "WAVE" # Assigning time variable
  data[data == missing] <- NA # Makes NA the missingness indicator

  # Exposure summary
  exposure_summary <- data %>%
    filter(WAVE %in% exposure_time_pts) %>%
    group_by(WAVE) %>%
    summarize_at(vars(all_of(exposure)), list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)

  cat(knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"), format = 'pipe'), sep = "\n")
  knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"), format = 'html') %>%
    kableExtra::kable_styling() %>%
    kableExtra::save_kable(file = file.path(home_dir, paste0(exposure, "_exposure_info.html")))

  cat(paste0(exposure, " exposure descriptive statistics have now been saved in the home directory"), "\n")
  cat("\n")

  # Exposure history summary
  exposure_epochs$epochs <- as.character(exposure_epochs$epochs)
  time_varying_wide <- apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep = "", collapse = ".")
  time_varying_wide <- sort(time_varying_wide)
  time_varying_wide <- c(id, time_varying_wide)
  time_varying_wide <- time_varying_wide[!time_varying_wide %in% time_var_exclude]
  data_wide <- suppressWarnings(stats::reshape(data = data, idvar = id, v.names = time_varying_covariates, timevar = "WAVE", times = c(time_pts), direction = "wide"))
  new <- data.frame(ID = data_wide[, id])
  colnames(new) <- id
  # Averages exposure across time points that constitute the exposure epochs (e.g., infancy = 6 & 15)
  for (e in 1:nrow(exposure_epochs)) {
    epoch <- exposure_epochs[e, 1]
    temp <- data.frame(row.names = 1:nrow(data_wide))
    new_var <- paste0(exposure, "_", epoch)
    # Finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
    for (l in 1:length(as.numeric(unlist(exposure_epochs[e, 2])))) {
      level <- as.numeric(unlist(exposure_epochs[e, 2]))[l]
      z <- data_wide[, which(grepl(paste0(exposure, ".", level), names(data_wide)))]
      temp <- cbind(temp, z)
    }
    new <- new %>%
      mutate(!!new_var := rowMeans(temp, na.rm = TRUE)) %>%
      mutate_at(vars(all_of(new_var)), factor)
  }
  # Assigning history (e.g., h-h-h) based on user-specified hi/lo cutoffs
  tot_hist <- apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                    paste, sep = "", collapse = "-")
  new$history <- apply(new, 1, function(x) {
    paste(lapply(1:nrow(exposure_epochs), function(y) {
      if (is.na(new[x, y + 1])) {
        return(NA)
      }
      if (new[x, y + 1] >= as.numeric(quantile(new[, y + 1], probs = hi_cutoff, na.rm = TRUE))) {
        return("h")
      }
      if (new[x, y + 1] <= as.numeric(quantile(new[, y + 1], probs = lo_cutoff, na.rm = TRUE))) {
        return("l")
      }
    }), collapse = "-")
  })

    # Summarizing n's by history
    his_summ <- new %>%
      group_by(history) %>%
      summarize(n = n()) %>%
      filter(!is.na(history)) %>%
      filter(history != "NULL")

    cat(paste0("USER ALERT: Out of the total of ", nrow(data_wide), " individuals in the sample, below is the distribution of the ", sum(his_summ$n), " (",
               round((sum(his_summ$n) / nrow(data_wide)) * 100, 2), "%) that fall into ", nrow(his_summ), " out of the ", length(tot_hist),
               " the total user-defined exposure histories created from ",
               lo_cutoff * 100, "th and ", hi_cutoff * 100, "th percentile values for low and high levels of exposure ", exposure,
               ", respectively, across ", paste(exposure_epochs$epochs, collapse = ", "),
               ". Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision :"), "\n")

    if (nrow(his_summ) != length(tot_hist)) {
      cat(paste0("USER ALERT: There are no individuals in your sample that fall into ", tot_hist[!tot_hist %in% his_summ$history],
                 " exposure history/histories. You may wish to consider different high/low cutoffs or choose a different measure to avoid extrapolation."), "\n")
    }

    cat(knitr::kable(his_summ, caption = paste0("Summary of User-Specified Exposure ", exposure, " Histories Based on Exposure Epochs"), format = 'pipe', row.names = F), sep = "\n")
    cat("\n")

    # Outcome summary
    outcome_summary <- data %>%
      # filter(WAVE %in% outcome_time_pt) %>%
      group_by(WAVE) %>%
      summarize_at(vars(all_of(outcome)), list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)

    cat(knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ", outcome, " Information"), format = 'pipe'), sep = "\n")
    knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ", outcome, " Information"), format = 'html') %>%
      kableExtra::kable_styling() %>%
      kableExtra::save_kable(file = file.path(home_dir, paste0(outcome, "_outcome_info.html")))

    cat(paste0(outcome, " outcome descriptive statistics have now been saved in the home directory"), "\n")

    if (sum(factor_covariates %in% colnames(data)) < length(factor_covariates)) {
      stop('Please provide factor covariates that correspond to columns in your data when creating the msm object')
    }

    # Formatting factor covariates
    data[, factor_covariates] <- lapply(data[, factor_covariates], factor)
    data[, id] <- factor(data[, id])

    # Formatting numeric covariates
    numeric_vars <- colnames(data)[!colnames(data) %in% c(factor_covariates, id)]
    data[, numeric_vars] <- lapply(data[, numeric_vars], as.numeric)

    return(data)
}
