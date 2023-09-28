#' Creates IPTW balancing weights
#'
#' Creates IPTW balancing weights at each user-specified exposure time point
#' using balancing formulas that relate exposure at each time point to all
#' relevant confounders.
#'
#' @export
#' @importFrom SuperLearner SuperLearner
#' @seealso {[WeightIt::weightitMSM()],
#'   <https://ngreifer.github.io/WeightIt/reference/weightitMSM.html>}
#' @param home_dir path to home directory (required if 'save.out' = TRUE)
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
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
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    method = "cbps",
#'                    save.out = FALSE)


createWeights <- function(home_dir, data, exposure, outcome, formulas, method = "cbps", SL.library = "SL.glm", criterion = "p.mean",
                          density = "dt_2", read_in_from_file = "no", verbose = TRUE, save.out = TRUE) {

  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if(!is.character(home_dir)){
      stop("Please provide a valid home directory path as a string if you wish to save output locally.", call. = FALSE)
    }
    else if(!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.", call. = FALSE)
    }
  }

  if (missing(data)){
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  else if (!inherits(data, "mids") && !is.data.frame(data) && !is.list(data)) {
    stop("Please provide either a 'mids' object, a data frame, or a list of imputed data frames in the 'data' field.",
         call. = FALSE)
  }
  else if(is.list(data) && !is.data.frame(data)){
    if (sum(sapply(data, is.data.frame)) != length(data)){
      stop("Please supply a list of data frames that have been imputed.", call. = FALSE)
    }
  }

  if (missing(exposure)){
    stop("Please supply a single exposure.", call. = FALSE)
  }
  else if(!is.character(exposure) || length(exposure) != 1){
    stop("Please supply a single exposure as a character.", call. = FALSE)
  }

  if (missing(outcome)){
    stop("Please supply a single outcome.", call. = FALSE)
  }
  else if(!is.character(outcome) || length(outcome) != 1){
    stop("Please supply a single outcome as a character.", call. = FALSE)
  }

  if (missing(formulas)){
    stop("Please supply a list of balancing formulas.", call. = FALSE)
  }
  else if(!is.list(formulas) | is.data.frame(formulas)){
    stop("Please provide a list of formulas for each exposure time point", call. = FALSE)
  }

  if(!is.character(method)){
    stop("Please provide as a character string a weights method from this list: 'ps', 'glm', 'gbm', 'bart', 'super', 'cbps'.",
         call. = FALSE)
  }
  else if(! method %in% c("ps", "glm", "gbm", "bart", "super", "cbps")){
    stop("Please provide a weights method from this list: 'ps', 'glm', 'gbm', 'bart', 'super', 'cbps'.",
         call. = FALSE)
  }


  if(!is.logical(verbose)){
    stop("Please set verbose to either TRUE or FALSE.", call. = FALSE)
  }
  else if(length(verbose) != 1){
    stop("Please provide a single TRUE or FALSE value to verbose.", call. = FALSE)
  }

  if(!is.logical(save.out)){
    stop("Please set save.out to either TRUE or FALSE.", call. = FALSE)
  }
  else if(length(save.out) != 1){
    stop("Please provide a single TRUE or FALSE value to save.out.", call. = FALSE)
  }



  weights_method <- method
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[", 1)

  if(save.out){
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
      # weights <- readRDS(paste0(home_dir, "/weights/", exposure, "-", outcome, "_", form_name, "_",
      #                           weights_method, "_fit.rds"))
      weights <- readRDS(sprintf("%s/weights/%s-%s_%s_%s_fit.rds",
                                 home_dir, exposure, outcome, form_name, weights_method))

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
    form <- formulas
    form <- unname(form)

    # Helper function to calculate weights
    calculate_weights <- function(data, form, weights_method, SL.library, criterion, density, verbose) {

      if(weights_method == "super"){
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     density = density, #continuous exposures
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     SL.library = SL.library,
                                     verbose = verbose,
                                     over = FALSE)
      }
      else if (weights_method == "glm"){
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     density = density, #continuous exposures
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     verbose = verbose,
                                     over = FALSE)

      }
      else if (weights_method == "gbm"){
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     density = density ,#continuous exposures
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     criterion = criterion,
                                     verbose = verbose,
                                     over = FALSE)
      }
      else{
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     density = density, #continuous exposures
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     verbose = verbose,
                                     over = FALSE)
      }
      fit
    }


    if (inherits(data, "mids")) {

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

        if (verbose){

          cat(sprintf("For imputation %s and the %s weighting method, the median weight value is %s
                     (SD= %s; range= %s-%s ). \n",
                      i,
                      weights_method,
                      round(median(fit$weights), 2),
                      round(sd(fit$weights), 2),
                      round(min(fit$weights), 2),
                      round(max(fit$weights), 2)))

          cat('\n')
        }

        if(save.out){
          # Save weights merged with ID variable
          write.csv(x = d, file = sprintf("%s/weights/values/%s-%s_%s_%s_%s.csv",
                                          home_dir, exposure, outcome, form_name, weights_method, i))
        }

        # Writes image of the histogram of weights to assess heavy tails
        p <- ggplot2::ggplot(data = as.data.frame(fit$weight), ggplot2::aes(x = fit$weight)) +
          ggplot2::geom_histogram(color = 'black', bins = 15) +
          ggplot2::xlab("Weights") +
          ggplot2::ggtitle(sprintf("Distribution of %s weights",
                                   weights_method))

        if(verbose){
          print(p)
        }

        if(save.out){
          ggsave(sprintf("%s/weights/histograms/Hist_%s-%s_%s_%s_%s.png",
                         home_dir, exposure, outcome, form_name, weights_method, i),
                 plot = p, height = 8, width = 14)
        }

        fit
      })

      if (save.out & verbose){
        cat("Weights for each imputation have now been saved into the 'weights/values/' folder.")
        cat("\n")
        cat("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely have heavy tails.")
      }
    }

    else if(is.list(data) && !is.data.frame(data)){
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

        fit <- calculate_weights(d, form, weights_method, SL.library, criterion, density, verbose)

        d$weights <- fit$weights

        if (verbose){
          # cat(paste0("For imputation", i, " and ", weights_method,
          #            ", weighting method, the median weight value is ", round(median(fit$weights), 2) ,
          #            " (SD= ", round(sd(fit$weights), 2), "; range= ", round(min(fit$weights), 2), "-",
          #            round(max(fit$weights), 2), ")."), "\n")
          cat(sprintf("For imputation %s and the %s weighting method, the median weight value is %s
                     (SD= %s; range= %s-%s ). \n",
                      i,
                      weights_method,
                      round(median(fit$weights), 2),
                      round(sd(fit$weights), 2),
                      round(min(fit$weights), 2),
                      round(max(fit$weights), 2)))
          cat('\n')
        }

        if(save.out){
          # Save weights merged with ID variable
          write.csv(x = d,
                    file = sprintf("%s/weights/values/%s-%s_%s_%s_%s.csv",
                                   home_dir, exposure, outcome, form_name, weights_method, i))
        }

        # Writes image of the histogram of weights to assess heavy tails
        p <- ggplot2::ggplot(data = as.data.frame(fit$weight), ggplot2::aes(x = fit$weight)) +
          ggplot2::geom_histogram(color = 'black', bins = 15) +
          ggplot2::xlab("Weights") +
          ggplot2::ggtitle(paste0("Distribution of ", weights_method, " weights"))

        if(verbose){
          print(p)
        }

        if(save.out){
          ggplot2::ggsave(sprintf("%s/weights/histograms/Hist_%s-%s_%s_%s_%s.png",
                                  home_dir, exposure, outcome, form_name, weights_method, i),
                          plot = p, height = 8, width = 14)
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
        calculate_weights(data, form, weights_method, SL.library, criterion, density, verbose)
      })

      data$weights <- weights[[1]]$weights

      if (verbose){
        # cat(paste0("For the ", weights_method, " weighting method, the median weight value is ",
        #            round(median(data$weights), 2) , " (SD = ", round(sd(data$weights), 2), "; range = ",
        #            round(min(data$weights), 2), "-", round(max(data$weights), 2), ")."), "\n")
        cat(sprintf("For the %s weighting method, the median weight value is %s
                    (SD = %s; range = %s-%s). \n",
                    weights_method,
                    round(median(data$weights), 2),
                    round(sd(data$weights), 2),
                    round(min(data$weights), 2),
                    round(max(data$weights))))

        cat('\n')
      }

      names(weights) <- "0"

      if(save.out){
        # Save weights merged with ID variable
        write.csv(x = data, file =  sprintf("%s/weights/values/%s-%s_%s_%s.csv",
                                            home_dir, exposure, outcome, form_name, weights_method))
      }

      # Writes image of the histogram of weights to assess heavy tails
      p <- ggplot2::ggplot(data = as.data.frame(weights[[1]]$weights),
                           ggplot2::aes(x = weights[[1]]$weights)) +
        ggplot2::geom_histogram(color = 'black', bins = 15) +
        ggplot2::xlab("Weights") +
        ggplot2::ggtitle(
          # paste0("Distribution of ", weights_method, " weights"))
          sprintf("Distribution of %s weights",
                  weights_method))

      if(verbose){
        print(p)
      }

      if(save.out){
        ggplot2::ggsave(sprintf("%s/weights/histograms/Hist_%s-%s_%s_%s.png",
                                home_dir, exposure, outcome, form_name, weights_method),
                        plot = p, height = 8, width = 14)

        if (verbose){
          cat("Weights have now been saved into the 'weights/values/' folder.")
          cat("\n")
          cat("Weights histograms have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
        }
      }
    }

    if(save.out){
      saveRDS(weights, file = sprintf("%s/weights/%s-%s_%s_%s_fit.rds",
                                      home_dir, exposure, outcome, form_name, weights_method))
      if (verbose){
        cat("\n")
        cat("Weights models have been saved as an .rds object in the 'weights' folder.", "\n")
      }
    }

    weights
  }
}


