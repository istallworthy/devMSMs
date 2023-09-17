#' Creates IPTW balancing weights
#'
#' Creates IPTW balancing weights at each user-specified exposure time point
#' using balancing formulas that relate exposure at each time point to all
#' relevant confounders.
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 ggsave
#' @importFrom WeightIt weightitMSM
#' @seealso {[WeightIt::WeightItMSM()], <url1>}
#' @param home_dir path to home directory
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix
#' @param formulas list of balancing formulas at each time point output from
#'   createFormulas()
#' @param method (optional) character string of WeightItMSM() balancing method
#'   abbreviation (default is Covariate Balancing Propensity Score "cbps")
#' @param SL.library required for superLearner weighting method; see
#'   SuperLearner::listWrappers() for options
#' @param read_in_from_file (optional) "yes" or "no" indicator to read in
#'   weights that have been previously run and saved locally (default is "no")
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return list of IPTW balancing weights
#' @export
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
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),,
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     save.out = FALSE)
#'
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                    formulas = f,
#'                    method = "cbps",
#'                    save.out = FALSE)


createWeights <- function(home_dir, data, exposure, outcome, tv_confounders, formulas, method = "cbps",
                          SL.library = NA, read_in_from_file = "no", verbose = TRUE, save.out = TRUE) {

  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if(!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.", call. = FALSE)
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
  if (missing(tv_confounders)){
    stop("Please supply a list of time-varying confounders.", call. = FALSE)
  }
  if (missing(formulas)){
    stop("Please supply a list of balancing formulas.", call. = FALSE)
  }

  if(!inherits(method, "character")){
    stop("Please provide as a character string a weights method from this list: 'ps', 'glm', 'gbm', 'bart', 'super', 'cbps'.", call. = FALSE)
  }
  if(! method %in% c("ps", "glm", "gbm", "bart", "super", "cbps")){
    stop("Please provide a weights method from this list: 'ps', 'glm', 'gbm', 'bart', 'super', 'cbps'.", call. = FALSE)
  }

  if (!mice::is.mids(data) & !is.data.frame(data) & !inherits(data, "list")) {
    stop("Please provide either a 'mids' object, a data frame, or a list of imputed data frames in the 'data' field.", call. = FALSE)
  }

  if(!inherits(formulas, "list")){
    stop("Please provide a list of formulas for each exposure time point", call. = FALSE)
  }

  # Load required libraries
  # library(cobalt)

  # exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)], "\\."), "[", 2))
  time_varying_covariates <- tv_confounders
  weights_method <- method
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[", 1)

  if(save.out){
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
  }

  if (read_in_from_file == "yes") {
    tryCatch({
      weights <- readRDS(paste0(home_dir, "/weights/", exposure, "-", outcome, "_", form_name, "_",
                                weights_method, "_fit.rds"))

      if (verbose){
        message("Reading in balancing weights from the local folder.")
      }
    }, error = function(x) {
      stop("These weights have not previously been saved locally. Please re-run with read_in_from_file='no'",
           call. = FALSE)
    })
    weights
  }
  else {

    # List of formulas for each time point
    # form <- formulas[grepl(paste("form_", exposure, "-", outcome, sep = ""), names(formulas))]
    form <- formulas
    form <- unname(form)

    # Helper function to calculate weights
    calculate_weights <- function(data, form, weights_method, SL.library) {

      if(weights_method == "super"){
        fit <- weightitMSM(form,
                           data = data,
                           method = weights_method,
                           stabilize = TRUE,
                           density = "dt_2", #do we want this?
                           use.kernel = TRUE,
                           include.obj = TRUE,
                           SL.library = SL.library,
                           over = FALSE)
      }
      else{
        fit <- weightitMSM(form,
                           data = data,
                           method = weights_method,
                           stabilize = TRUE,
                           density = "dt_2", #do we want this?
                           use.kernel = TRUE,
                           include.obj = TRUE,
                           over = FALSE)
      }
      fit
    }


    if(mice::is.mids(data)){
      # Cycling through imputed datasets
      weights <- lapply(seq_len(data$m), function(i) {
        d <- as.data.frame(mice::complete(data, i))

        if (length(which(is.na(d))) > 0){
          stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
               call. = FALSE)
        }
        if (sum(duplicated(d$"ID")) > 0){
          stop("Please provide wide imputed datasets with a single row per ID.", call. = FALSE)
        }

        fit <- calculate_weights(d, form, weights_method, SL.library)

        d$weights <- fit$weights

        cat(paste0("USER ALERT: For imputation ", i, " and the ", weights_method, " weighting method, the median weight value is ",
                   round(median(fit$weights), 2) ," (SD= ", round(sd(fit$weights), 2), "; range= ", round(min(fit$weights), 2), "-",
                   round(max(fit$weights), 2), ")."), "\n")
        cat('\n')

        if(save.out){
          # Save weights merged with ID variable
          write.csv(x = d, file = paste0(home_dir, "/weights/values/", exposure, "-", outcome, "_", form_name,
                                         "_", weights_method, "_", i, ".csv"))
        }

        # Writes image of the histogram of weights to assess heavy tails
        p <- ggplot2::ggplot(data = as.data.frame(fit$weight), aes(x = fit$weight)) +
          ggplot2::geom_histogram(color = 'black', bins = 15) +
          ggplot2::xlab("Weights") +
          ggplot2::ggtitle(paste0("Distribution of ", weights_method, " weights"))

        if(verbose){
          print(p)
        }

        if(save.out){
          ggsave(paste0(home_dir, "/weights/histograms/", "Hist_", exposure, "-", outcome, "_", form_name,
                        "_", weights_method, "_", i, ".png"), plot = p, height = 8, width = 14)
        }

        fit
      })

      if (save.out & verbose){
        cat("Weights for each imputation have now been saved into the 'weights/values/' folder.")
        cat("\n")
        cat("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely have heavy tails.")
      }
    }

    else if(inherits(data, "list")){
      # Cycling through list of imputed datasets
      weights <- lapply(seq_len(length(data)), function(i) {
        d <- data[[i]]

        if (length(which(is.na(d))) > 0){
          stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
               call. = FALSE)
        }
        if (sum(duplicated(d$"ID")) > 0){
          stop("Please provide wide imputed datasets with a single row per ID.", call. = FALSE)
        }

        fit <- calculate_weights(d, form, weights_method, SL.library)

        d$weights <- fit$weights

        if (verbose){
          message(paste0("USER ALERT: For imputation", i, " and ", weights_method,
                         ", weighting method, the median weight value is ", round(median(fit$weights), 2) ,
                         " (SD= ", round(sd(fit$weights), 2), "; range= ", round(min(fit$weights), 2), "-",
                         round(max(fit$weights), 2), ")."), "\n")
          cat('\n')
        }

        if(save.out){
          # Save weights merged with ID variable
          write.csv(x = d, file = paste0(home_dir, "/weights/values/", exposure, "-", outcome,
                                         "_", form_name, "_", weights_method, "_", i, ".csv"))
        }

        # Writes image of the histogram of weights to assess heavy tails
        p <- ggplot2::ggplot(data = as.data.frame(fit$weight), aes(x = fit$weight)) +
          ggplot2::geom_histogram(color = 'black', bins = 15) +
          ggplot2::xlab("Weights") +
          ggplot2::ggtitle(paste0("Distribution of ", weights_method, " weights"))

        if(verbose){
          print(p)
        }

        if(save.out){
          ggplot2::ggsave(paste0(home_dir, "/weights/histograms/", "Hist_", exposure, "-", outcome,
                                 "_", form_name, "_", weights_method, "_", i, ".png"), plot = p,
                          height = 8, width = 14)
        }

        fit
      })

      if (save.out & verbose){
        cat("Weights for each imputation have now been saved into the 'weights/values/' folder.")
        cat("\n")
        cat("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
      }
    }

    else if (is.data.frame(data)){
      if (sum(duplicated(data$"ID")) > 0){
        stop("Please provide wide dataset with a single row per ID.", call. = FALSE)
      }
      if (length(which(is.na(data))) > 0){
        stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
             call. = FALSE)
      }

      # Creating weights
      weights <-  lapply(1, function(i) {
        calculate_weights(data, form, weights_method, SL.library)
      })

      data$weights <- weights[[1]]$weights

      if (verbose){
        message(paste0("USER ALERT: for the ", weights_method, " weighting method, the median weight value is ",
                       round(median(data$weights), 2) , " (SD = ", round(sd(data$weights), 2), "; range = ",
                       round(min(data$weights), 2), "-", round(max(data$weights), 2), ")."), "\n")
        cat('\n')
      }

      names(weights) <- "0"

      if(save.out){
        # Save weights merged with ID variable
        write.csv(x = data, file = paste0(home_dir, "/weights/values/", exposure, "-", outcome,
                                          "_", form_name, "_", weights_method, ".csv"))
      }

      # Writes image of the histogram of weights to assess heavy tails
      p <- ggplot2::ggplot(data = as.data.frame(weights[[1]]$weights), aes(x = weights[[1]]$weights)) +
        ggplot2::geom_histogram(color = 'black', bins = 15) +
        ggplot2::xlab("Weights") +
        ggplot2::ggtitle(paste0("Distribution of ", weights_method, " weights"))

      if(verbose){
        print(p)
      }

      if(save.out){
        ggplot2::ggsave(paste0(home_dir, "/weights/histograms/", "Hist_", exposure, "-", outcome,
                               "_", form_name, "_", weights_method, ".png"), plot = p, height = 8, width = 14)
        if (verbose){
          cat("Weights have now been saved into the 'weights/values/' folder.")
          cat("\n")
          cat("Weights histograms have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
        }
      }
    }

    if(save.out){
      saveRDS(weights, file = paste0(home_dir, "/weights/", exposure, "-", outcome,
                                     "_", form_name, "_", weights_method, "_fit.rds"))
      if (verbose){
        message("Weights models have been saved as an .rds object in the 'weights' folder.")
      }
    }

    weights
  }
}


