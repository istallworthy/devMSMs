#' Creates balancing weights for potential confounding covariates in relation to exposure at each time point
#'
#' Returns a list of weights_models
#' @param object msm object that contains all relevant user inputs
#' @param wide_long_datasets from formatForWeights
#' @param short_forms from createShortForms
#' @param read_in_from_file optional binary indicator of whether weights should be read in from a local file
#' @return weights_models
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 ggsave
#' @importFrom WeightIt weightitMSM
#' @importFrom cobalt bal.tab
#' @seealso [CBPS::CBPS()], [formatForWeights()], [createShortForms()]
#' @examples createWeights(object, wide_long_datasets, short_forms, read_in_from_file="no")
createWeights <- function(data, exposure, outcome, tv_confounders, formulas, method = "cbps", read_in_from_file = "no") {
  # Load required libraries
  library(cobalt)
  library(ggplot2)
  library(WeightIt)

  # Extracting arguments from the msm object
  ID <- "S_ID"
  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))
  time_varying_covariates <- tv_confounders
  weights_method <- method

  # Determining exposure type
  exposure_type <- if (length(unique(d[, paste0(exposure, ".", exposure_time_pts[1])])) < 3) "binary" else "continuous"

  # Getting form type
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[",1)


  if(! method %in% c("ps", "glm", "gbm", "bart", "super", "cbps")){
    stop("Please provide a weights method from this list: 'ps', 'glm', 'gbm', 'bart', 'super', 'cbps'.")
  }

  weights_dir <- file.path(home_dir, "weights")
  if (!dir.exists(weights_dir)) {
    dir.create(weights_dir)
  }

  values_dir <- file.path(home_dir, "weights", "values")
  if (!dir.exists(values_dir)) {
    dir.create(values_dir)
  }

  hist_dir <- file.path(home_dir, "weights", "histograms")
  if (!dir.exists(hist_dir)) {
    dir.create(hist_dir)
  }

  if (read_in_from_file == "yes") {
    tryCatch({
      dat_w <- readRDS(paste0(home_dir, "weights/", exposure, "-", outcome, "_", form_name, "_", weights_method, "_fits.rds"))
      message("Reading in balancing weights from the local folder.")
    }, error = function(x) {
      stop("These weights have not previously been saved locally. Please re-run with read_in_from_file='no'")
    })
    return(dat_w)
  }

  message(paste0("Creating longitudinal balancing weights using the ", weights_method, " weighting method and the ", form_name, " formulas"))

  # List of formulas for each time point
  form <- formulas[grepl(paste("form_", exposure, "-", outcome, sep = ""), names(formulas))]
  form <- unname(form)

  # Helper function to calculate weights
  calculate_weights <- function(data, form, weights_method) {
    set.seed(1234)
    fit <- weightitMSM(form, data = data,
                       method = weights_method,
                       stabilize = TRUE,
                       density = "dt_2",
                       use.kernel = TRUE,
                       weightit.force = TRUE,
                       over = FALSE)
    return(fit)
  }


  if(class(data) == "mids"){
    # Cycling through imputed datasets
    weights <- lapply(1:data$m, function(i) {
      d <-as.data.frame(mice::complete(data,i))
      v <- sapply(strsplit(tv_confounders, "\\."), "[",1)
      v <- v[!duplicated(v)]
      d_w <- suppressWarnings(stats::reshape(data = d,
                                             idvar = ID,
                                             v.names = v,
                                             timevar = "WAVE",
                                             times = c(exposure_time_pts),
                                             direction = "wide"))
      # Creating weights
      fit <- calculate_weights(d_w, form, weights_method)

      d$weights <- fit$weights

      # Save weights merged with ID variable
      write.csv(x = d, file = paste0(home_dir, "/weights/values/", exposure, "-", outcome, "_", form_name, "_", weights_method, "_", i, ".csv"))

      # Writes image of the histogram of weights to assess heavy tails
      ggplot(data = as.data.frame(fit$weight), aes(x = fit$weight)) +
        geom_histogram(color = 'black', bins = 15)
      ggsave(paste0(home_dir, "/weights/histograms/", "Hist_", exposure, "-", outcome, "_", form_name, "_", weights_method, "_", i, ".png"),
             height = 8, width = 14)

      fit
    })
    message("Weights for each imputation have now been saved into the 'weights/values/' folder.")
    message("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
  }

  if (class(data) == "data.frame"){
    # Creating weights
    weights <- calculate_weights(data, form, weights_method)

    data$weights <- fit$weights

    # Save weights merged with ID variable
    write.csv(x = data, file = paste0(home_dir, "/weights/values/", exposure, "-", outcome, "_", form_name, "_", weights_method, ".csv"))

    # Writes image of the histogram of weights to assess heavy tails
    ggplot(data = as.data.frame(fit$weight), aes(x = fit$weight)) +
      geom_histogram(color = 'black', bins = 15)
    ggsave(paste0(home_dir, "/weights/histograms/", "Hist_", exposure, "-", outcome, "_", form_name, "_", weights_method, ".png"),
           height = 8, width = 14)

    message("Weights have now been saved into the 'weights/values/' folder.")
    message("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
  }

  saveRDS(weights, file = paste0(home_dir, "/weights/", exposure, "-", outcome, "_", form_name, "_", weights_method, "_fit.rds"))
  message("Weights models have been saved as an .rds object in the 'weights' folder.")

  fit
}


