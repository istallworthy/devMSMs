#' Creates IPTW balancing weights
#'
#' Creates IPTW balancing weights at each user-specified exposure time point
#' using balancing formulas that relate exposure at each time point to all
#' relevant confounders.
#'
#' @export
#' @seealso {[WeightIt::weightitMSM()],
#'   <https://ngreifer.github.io/WeightIt/reference/weightitMSM.html>}
#' @param home_dir path to home directory (required if 'save.out' = TRUE)
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param formulas list of balancing formulas at each time point output from
#'   createFormulas()
#' @param method (optional) character string of weightitMSM() balancing method
#'   abbreviation (default is Covariate Balancing Propensity Score "cbps")
#' @param SL.library required for superLearner weighting method ("super"); see
#'   SuperLearner::listWrappers() for options
#' @param criterion (optional) criterion used to select best weights (default is
#'   "p.mean" minimizing avg Pearson correlation for continuous exposures and
#'   "smd.mean" for binary exposures) (requird for "gbm" method)
#' @param read_in_from_file (optional) TRUE or FALSE indicator to read in
#'   weights that have been previously run and saved locally (default is FALSE)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param ... for other inputs to weightitMSM()
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
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     save.out = FALSE)
#' @examplesIf requireNamespace("CBPS", quietly = TRUE)
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    save.out = FALSE)
#'                    
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3", "A.1:B.1"),
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     save.out = FALSE)
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
#' 
#' @examplesIf requireNamespace("CBPS", quietly = TRUE)
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    method = "cbps",
#'                    save.out = FALSE)
#' 
#' @examplesIf requireNamespace("gbm", quietly = TRUE)                   
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    method = "gbm",
#'                    save.out = FALSE)                    
#' 
#' @examplesIf requireNamespace("dbarts", quietly = TRUE)
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    method = "bart",
#'                    save.out = FALSE)
#' 
#' @examplesIf requireNamespace("SuperLearner", quietly = TRUE)
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    method = "super",
#'                    save.out = FALSE)   
#'                                                      
createWeights <- function(data, exposure, outcome, formulas, method = "cbps",
                          SL.library = "SL.glm", criterion = NA, home_dir = NULL,
                          read_in_from_file = FALSE, verbose = TRUE, save.out = TRUE, ...) {
  
  # call <- match.call()
  
  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.",
           call. = FALSE)
    }
    if (!is.character(home_dir)) {
      stop("Please provide a valid home directory path as a string if you wish to save output locally.",
           call. = FALSE)
    }
    if (!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.",
           call. = FALSE)
    }
  }
  
  if (missing(data)) {
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  if (!inherits(data, "mids") && !is.data.frame(data) && !is.list(data)) {
    stop("Please provide either a 'mids' object, a data frame, or a list of imputed data frames in the 'data' field.",
         call. = FALSE)
  }
  if (is.list(data) && !is.data.frame(data)  && !inherits(data, "mids") &&
      !all(sapply(data, is.data.frame))) {
    stop("Please supply a list of data frames that have been imputed.",
         call. = FALSE)
  }
  
  if (inherits(data, "mids")) {
    rlang::check_installed("mice")
  }
  
  if (missing(exposure)) {
    stop("Please supply a single exposure.",
         call. = FALSE)
  }
  if (!is.character(exposure) || length(exposure) != 1) {
    stop("Please supply a single exposure as a character.",
         call. = FALSE)
  }
  else if (grepl("\\.", exposure)) {
    stop ("Please supply an exposure without the '.time' suffix or any '.' special characters. Note that the exposure variables in your dataset should be labeled with the '.time' suffix.",
          call. = FALSE)
  }
  
  if (missing(outcome)) {
    stop("Please supply a single outcome.",
         call. = FALSE)
  }
  if (!is.character(outcome) || length(outcome) != 1) {
    stop("Please supply a single outcome as a character.",
         call. = FALSE)
  }
  else if (!grepl("\\.", outcome)) {
    stop("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
         call. = FALSE)
  }
  
  if (missing(formulas)) {
    stop("Please supply a list of balancing formulas.",
         call. = FALSE)
  }
  if (!is.list(formulas) || is.data.frame(formulas)) {
    stop("Please provide a list of formulas for each exposure time point",
         call. = FALSE)
  }
  
  if (!is.character(method) || length(method) != 1) {
    stop("Please provide as a character string a weights method from this list: 'glm', 'gbm', 'bart', 'super', 'cbps'.",
         call. = FALSE)
  }
  if (!method %in% c("glm", "gbm", "bart", "super", "cbps")) {
    stop("Please provide a weights method from this list: 'glm', 'gbm', 'bart', 'super', 'cbps'.",
         call. = FALSE)
  }
  
  if (!is.logical(read_in_from_file)) {
    stop("Please set read_in_from_file to either TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(read_in_from_file) != 1) {
    stop("Please provide a single TRUE or FALSE value to read_in_from_file.",
         call. = FALSE)
  }
  
  if (!is.logical(verbose)) {
    stop("Please set verbose to either TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(verbose) != 1) {
    stop("Please provide a single TRUE or FALSE value to verbose.",
         call. = FALSE)
  }
  
  if (!is.logical(save.out)) {
    stop("Please set save.out to either TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(save.out) != 1) {
    stop("Please provide a single TRUE or FALSE value to save.out.",
         call. = FALSE)
  }
  
  weights_method <- method
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[", 1)
  
  if (save.out) {
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
  
  if (read_in_from_file) {
    
    tryCatch ( {
      weights <- readRDS(file.path(home_dir, "weights", 
                                   sprintf("%s-%s_%s_%s_fit.rds",
                                           exposure, outcome, form_name, 
                                           weights_method)))
      
      if (verbose) {
        message("Reading in balancing weights from the local folder.")
      }
      
      if (!is.list(weights) || !inherits(weights[[1]], "weightitMSM")) {
        stop("The weights saved locally are not the correct class of weightitMSM classed objects saved as a list.",
             call. = FALSE)
      }
      
      #checking weights
      
      if (is.data.frame(data)) {
        l <- 1
        n <- names(data)
      }
      else if (inherits(data, "mids")) {
        l <- data$m
        n <- names(mice::complete(data, 1))
      }
      else if (is.list(data) && !is.data.frame(data)) {
        l <- length(data)
        n <- names(data[[1]])
      }
      
      if (length(weights) != l) {
        stop("The locally saved weights object does not match the data you have supplied.",
             call. = FALSE)
      }
      
      if (weights[[1]]$method != weights_method) {
        stop(sprintf("The weights saved locally are incorrectly labeled and were actually created with %s weighting method.",
                     weights[[1]]$method),
             call. = FALSE)
      }
      
      if (length(weights[[1]]$covs.list) != length(formulas)) {
        stop("The locally saved weights were created using formulas other than what are supplied here.",
             call. = FALSE)
      }
      
      if (as.logical(unlist(lapply(weights[[1]]$covs.list, function(x) {
        !all(names(x) %in% n)
      })))) {
        stop("The locally saved weights were created using data other than what is supplied.",
             call. = FALSE)
      }
      
      
    }, error = function(x) {
      stop("These weights have not previously been saved locally. Please re-run with read_in_from_file = FALSE",
           call. = FALSE)
    })
  }
  else {
    
    # List of formulas for each time point
    
    form <- unname(formulas)
    
    
    #temp workaround while noah is fixing weightitMSM bug...this will only work for continuous exposures tho
    
    if (weights_method == "gbm" && is.na(criterion)) {
      criterion <- "p.mean"
    }
    
    # function to wrap weightitMSM() and calculate weights
    
    calculate_weights <- function(data, form, weights_method, SL.library, 
                                  criterion, verbose, ...) {
      
      
      
      #fitting weights model
      
      if (weights_method == "super") {
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     include.obj = TRUE,
                                     SL.library = SL.library, #required
                                     verbose = verbose,
                                     ...)
      }
      else if (weights_method == "glm") {
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     verbose = verbose,
                                     ...)
        
      }
      else if (weights_method == "gbm") {
        
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     criterion = criterion, #required? even tho doc says there is default
                                     verbose = verbose,
                                     ...)
      }
      else if (weights_method == "cbps") {
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     verbose = verbose, #for cbps
                                     over = FALSE,
                                     ...)
      }
      else {
        fit <- WeightIt::weightitMSM(form,
                                     data = data,
                                     method = weights_method,
                                     stabilize = TRUE,
                                     use.kernel = TRUE,
                                     include.obj = TRUE,
                                     verbose = verbose,
                                     over = FALSE,
                                     ...)
      }
      fit
    }
    
    
    if (inherits(data, "mids")) {
      
      # Cycling through imputed datasets
      
      weights <- lapply(seq_len(data$m), function(i) {
        d <- as.data.frame(mice::complete(data, i))
        
        if (anyNA(d)) {
          stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
               call. = FALSE)
        }
        if (any(duplicated(d[["ID"]]))) {
          stop("Please provide wide imputed datasets with a single row per ID.",
               call. = FALSE)
        }
        
        fit <- calculate_weights(d, form, weights_method, SL.library, 
                                 criterion, verbose, ...)
        
        d$weights <- fit$weights
        
        if (verbose) {
          
          cat(sprintf("For imputation %s and the %s weighting method, the median weight value is %s (SD= %s; range= %s-%s). \n",
                      i,
                      weights_method,
                      round(median(fit$weights), 2),
                      round(sd(fit$weights), 2),
                      round(min(fit$weights), 2),
                      round(max(fit$weights), 2)))
          
          cat('\n')
        }
        
        if (save.out) {
          
          # Save weights merged with ID variable
          
          write.csv(x = d, 
                    file = file.path(home_dir, "weights", "values", 
                                     sprintf("%s-%s_%s_%s_%s.csv",
                                             exposure, outcome, form_name, 
                                             weights_method, i)))
        }
        
        # Writes image of the histogram of weights to assess heavy tails
        
        p <- ggplot2::ggplot(data = as.data.frame(fit$weight), 
                             ggplot2::aes(x = fit$weight)) +
          ggplot2::geom_histogram(color = 'black', bins = 15) +
          ggplot2::xlab("Weights") +
          ggplot2::ggtitle(sprintf("Distribution of %s weights",
                                   weights_method))
        
        if (verbose) {
          print(p)
        }
        
        if (save.out) {
          ggplot2::ggsave(file.path(home_dir, "weights", "histograms", 
                                    sprintf("Hist_%s-%s_%s_%s_%s.png",
                                            exposure, outcome, form_name, 
                                            weights_method, i)),
                          plot = p, height = 8, width = 14)
        }
        
        fit
      })
      
      if (save.out & verbose) {
        cat("Weights for each imputation have now been saved into the 'weights/values/' folder.")
        cat("\n")
        cat("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely have heavy tails.")
      }
    }
    
    else if (is.list(data) && !is.data.frame(data)) {
      
      # Cycling through list of imputed datasets
      
      weights <- lapply(seq_len(length(data)), function(i) {
        d <- data[[i]]
        
        if (anyNA(d)) {
          stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
               call. = FALSE)
        }
        if (any(duplicated(d[["ID"]]))) {
          stop("Please provide wide imputed datasets with a single row per ID.",
               call. = FALSE)
        }
        
        fit <- calculate_weights(d, form, weights_method, SL.library, 
                                 criterion, verbose, ...)
        
        d$weights <- fit$weights
        
        if (verbose) {
          
          cat(sprintf("For imputation %s and the %s weighting method, the median weight value is %s (SD= %s; range= %s-%s ). \n",
                      i,
                      weights_method,
                      round(median(fit$weights), 2),
                      round(sd(fit$weights), 2),
                      round(min(fit$weights), 2),
                      round(max(fit$weights), 2)))
          cat('\n')
        }
        
        if (save.out) {
          
          # Save weights merged with ID variable
          
          write.csv(x = d,
                    file = file.path(home_dir, "weights", "values", 
                                     sprintf("%s-%s_%s_%s_%s.csv",
                                             exposure, outcome, form_name, 
                                             weights_method, i)))
        }
        
        # Writes image of the histogram of weights to assess heavy tails
        
        p <- ggplot2::ggplot(data = as.data.frame(fit$weight), 
                             ggplot2::aes(x = fit$weight)) +
          ggplot2::geom_histogram(color = 'black', bins = 15) +
          ggplot2::xlab("Weights") +
          ggplot2::ggtitle(paste0("Distribution of ", weights_method, 
                                  " weights"))
        
        if (verbose) {
          print(p)
        }
        
        if (save.out) {
          ggplot2::ggsave(file.path(home_dir, "weights", "histograms", 
                                    sprintf("Hist_%s-%s_%s_%s_%s.png",
                                            exposure, outcome, form_name, 
                                            weights_method, i)),
                          plot = p, height = 8, width = 14)
        }
        
        fit
      })
      
      if (save.out & verbose) {
        cat("Weights for each imputation have now been saved into the 'weights/values/' folder.")
        cat("\n")
        cat("Weights histograms for each imputation have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
      }
    }
    
    else if (is.data.frame(data)) {
      if (any(duplicated(data[["ID"]]))) {
        stop("Please provide wide dataset with a single row per ID.",
             call. = FALSE)
      }
      if (anyNA(data)) {
        stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
             call. = FALSE)
      }
      
      # Creating weights
      
      weights <- list(
        calculate_weights(data, form, weights_method, SL.library, criterion, 
                          verbose, ...)
      )
      
      data$weights <- weights[[1]]$weights
      
      if (verbose) {
        cat(sprintf("For the %s weighting method, the median weight value is %s (SD = %s; range = %s-%s). \n",
                    weights_method,
                    round(median(data$weights), 2),
                    round(sd(data$weights), 2),
                    round(min(data$weights), 2),
                    round(max(data$weights))))
        
        cat('\n')
      }
      
      names(weights) <- "0"
      
      if (save.out) {
        
        # Save weights merged with ID variable
        
        write.csv(x = data, 
                  file = file.path(home_dir, "weights", "values", 
                                   sprintf("%s-%s_%s_%s.csv",
                                           exposure, outcome, form_name, 
                                           weights_method)))
      }
      
      # Writes image of the histogram of weights to assess heavy tails
      
      p <- ggplot2::ggplot(data = as.data.frame(weights[[1]]$weights),
                           ggplot2::aes(x = weights[[1]]$weights)) +
        ggplot2::geom_histogram(color = 'black', bins = 15) +
        ggplot2::xlab("Weights") +
        ggplot2::ggtitle( sprintf("Distribution of %s weights",
                                  weights_method))
      
      if (verbose) {
        print(p)
      }
      
      if (save.out) {
        ggplot2::ggsave(file.path(home_dir, "weights", "histograms", 
                                  sprintf("Hist_%s-%s_%s_%s.png",
                                          exposure, outcome, form_name, 
                                          weights_method)),
                        plot = p, height = 8, width = 14)
        
        if (verbose) {
          cat("Weights have now been saved into the 'weights/values/' folder.")
          cat("\n")
          cat("Weights histograms have now been saved in the 'weights/histograms/' folder --likely has heavy tails.")
          cat("\n")
        }
      }
    }
    
    if (save.out) {
      saveRDS(weights, 
              file = file.path(home_dir, "weights", 
                               sprintf("%s-%s_%s_%s_fit.rds",
                                       exposure, outcome, form_name, 
                                       weights_method)))
      if (verbose) {
        cat("\n")
        cat("Weights models have been saved as an .rds object in the 'weights' folder.", "\n")
      }
    }
  }
  
  weights
}


