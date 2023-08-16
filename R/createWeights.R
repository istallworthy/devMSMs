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
createWeights <- function(home_dir, data, exposure, outcome, tv_confounders, formulas, method = "cbps", read_in_from_file = "no") {

  #error checking
  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.")
  }

  if(! method %in% c("ps", "glm", "gbm", "bart", "super", "cbps")){
    stop("Please provide a weights method from this list: 'ps', 'glm', 'gbm', 'bart', 'super', 'cbps'.")
  }

  if (!class(data) %in% c("mids", "data.frame", "character")) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.")
  }

  # Load required libraries
  library(cobalt)
  library(ggplot2)
  library(WeightIt)

  # Extracting arguments from the msm object

  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)], "\\."), "[",2))
  time_varying_covariates <- tv_confounders
  weights_method <- method
  exposure_type <- if (length(unique(d[, paste0(exposure, ".", exposure_time_pts[1])])) < 3) "binary" else "continuous"
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[",1)

  # creating directories
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
      weights = paste0(home_dir, "/weights/", exposure, "-", outcome, "_", form_name, "_", weights_method, "_fit.rds")
      message("Reading in balancing weights from the local folder.")
    }, error = function(x) {
      stop("These weights have not previously been saved locally. Please re-run with read_in_from_file='no'")
    })
    weights

  } else {

    message(paste0("Creating longitudinal balancing weights using the ", weights_method, "  weighting method and the ", form_name, " formulas"), "\n")
    cat("\n")

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

        if (sum(duplicated(d$"ID")) > 0){
          stop("Please provide wide imputed datasets with a single row per ID.")
        }

        fit <- calculate_weights(d_w, form, weights_method)

        d$weights <- fit$weights

        cat(paste0("USER ALERT: For imputation ", i, " and the ", exposure, "-", outcome, " relation, the median weight value is ",
                       round(median(fit$weights),2) ," (SD= ", round(sd(fit$weights),2), "; range= ", round(min(fit$weights),2), "-",
                       round(max(fit$weights),2), ")."), "\n")
        cat('\n')

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


    if(class(data) == "character"){
      if (!dir.exists(data)) {
        stop("Please provide a valid directory path with imputed datasets, a data frame, or a 'mids' object for the 'data' field.")
      }
      if (length(dir(data)) < 2) {
        stop("If you specify data as a directory, please supply more than 1 imputed dataset.")
      }

      files <- list.files(data, full.names = TRUE, pattern = "\\.csv")

      # Read and process imputed datasets
      data <- lapply(files, function(file) {
        imp_data <- read.csv(file)
        imp_data
      })

      # Cycling through list of imputed datasets
      weights <- lapply(1:length(data), function(i) {
        d <- data[[i]]

        if (sum(duplicated(d$"ID")) > 0){
          stop("Please provide wide imputed datasets with a single row per ID.")
        }

        fit <- calculate_weights(d, form, weights_method)

        d$weights <- fit$weights

        message(paste0("USER ALERT: For imputation", i, " and the ", exposure, "-", outcome, " relation, the median weight value is ", round(median(fit$weights),2) ,
                       " (SD= ", round(sd(fit$weights),2), "; range= ", round(min(fit$weights),2), "-", round(max(fit$weights),2), ")."), "\n")
        cat('\n')

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

      if (sum(duplicated(data$"ID")) > 0){
        stop("Please provide wide dataset with a single row per ID.")
      }

      # Creating weights
      weights <-  lapply(1, function(i) {
        calculate_weights(data, form, weights_method)
      })
      names(weights) <- "0"

      data$weights <- weights$weights

      message(paste0("USER ALERT: For  the ", exposure, "-", outcome, " relation, the median weight value is ", round(median(weights$weights),2) ,
                     " (SD= ", round(sd(weights$weights),2), "; range= ", round(min(weights$weights),2), "-", round(max(weights$weights),2), ")."), "\n")
      cat('\n')

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

    weights
  }
}


