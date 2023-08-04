#' Truncates weights
#'
#' Code to 'cut off' weights at the 90th percentile and populate all of those above the 90th percentile value to avoid heavy tails that can bias results.
#'
#' @param data_for_model_with_weights output from condenseWeights
#' @param msm_object msm object that contains all relevant user inputs
#' @return dat_w_t list of truncated weights for each cutoff value
#' @export
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 ggplot ggsave geom_histogram
#' @importFrom WeightIt trim
#' @examples truncateWeights(object, data_for_model_with_weights)
#'
truncateWeights <- function(data_for_model_with_weights, msm_object) {
  home_dir <- msm_object$home_dir
  exposure <- msm_object$exposure
  exposure_epochs <- msm_object$exposure_epochs
  outcome <- msm_object$outcome
  weights_percentile_cutoff <- msm_object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity <- as.numeric(unlist(strsplit(msm_object$weights_percentile_cutoffs_sensitivity, " ")))
  m <- msm_object$m
  ID <- msm_object$ID
  weights_method <- msm_object$weights_method

  if (weights_percentile_cutoff > 1 || weights_percentile_cutoff < 0) {
    stop('Please select a percentile cutoff between 0 and 1 in the msmObject')
  }

  dat_w <- data_for_model_with_weights

  # Creating list of all cutoff values to cycle through (user-specified plus 2 others for sensitivity checks)
  all_cutoffs <- c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  # Function to truncate weights and save data
  truncate_and_save <- function(cutoff) {
    data_to_return <- list()

    lapply(1:length(dat_w), function(t) { # Cycling through imputations
      k <- t
      data <- dat_w[[t]]
      name <- paste0(cutoff, "_weight_cutoff") # Labeling new weight column
      cutoff_weights <- WeightIt::trim(data[, "weights"], cutoff, lower = FALSE)

      # Create a data frame to store truncated weights
      truncated_weights <- data.frame(ID = data[, ID], .imp = k, truncated_weight = cutoff_weights)

      # Add exposure epochs
      for (e in 1:nrow(exposure_epochs)) {
        epoch <- exposure_epochs[e, 1]
        temp <- data.frame(row.names = 1:nrow(data))
        new_var <- paste0(exposure, "_", epoch)
        for (l in 1:length(as.numeric(unlist(exposure_epochs[e, 2])))) {
          level <- as.numeric(unlist(exposure_epochs[e, 2]))[l]
          z <- data[, which(grepl(paste0(exposure, ".", level), names(data)))]
          temp <- cbind(temp, z)
        }
        new <- data %>%
          mutate(.imp = k, !!new_var := rowMeans(temp, na.rm = TRUE)) %>%
          select(ID, .imp, all_of(new_var))
        new[[new_var]] <- as.numeric(new[[new_var]])
      }

      data <- data %>%
        bind_cols(truncated_weights) %>%
        select(-weights) # Drop the original weight column

      # Save histogram of new weights
      ggplot(data, aes(x = truncated_weight)) +
        geom_histogram(color = 'black', bins = 15) +
        ggsave(paste("Hist_", exposure, "-", outcome, "_weights_cutoff_", cutoff, "_", weights_method, "_imp_", k, ".png", sep = ""),
               path = paste0(home_dir, "final weights/histograms/"), height = 8, width = 14)

      # Save truncated weight data
      write_csv(data, paste0(home_dir, "final weights/", exposure, "-", outcome, "_weights_cutoff_", cutoff, "_", weights_method, "_imp_", k, ".csv"))

      # Add data to the list
      data_to_return[[as.character(cutoff)]] <- data
    })

    return(data_to_return)
  }

  # Truncate weights for each cutoff value
  dat_w_t <- lapply(all_cutoffs, truncate_and_save)
  names(dat_w_t) <- all_cutoffs

  saveRDS(dat_w_t, paste0(home_dir, "final weights/values/", exposure, "-", outcome, "_", weights_method, "_dat_w_t.rds"))

  cat("\n")
  cat("USER ALERT: final truncated weights (using the user-specified cutoff values and 2 other values for subsequent sensitivity analyses) have now each been saved as a dataset in 'final weights' folder","\n")
  cat(paste0("A histogram of weights for exposure ", exposure, " on ", outcome, ", all imputations has been saved to the 'final weights/histograms/' folder"),"\n")
  cat("\n")
  cat("\n")

  return(dat_w_t)
}
