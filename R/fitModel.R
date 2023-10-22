#' Fit outcome model
#'
#' Fits weighted marginal outcome model as a generalized linear model of the
#' user's choosing, relating exposure main effects to outcome using IPTW
#' weights.
#' @seealso [survey::svyglm()] for more on family/link specifications.
#'   
#' @param home_dir path to home directory (required if 'save.out' = TRUE)
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param weights list of IPTW weights output from createWeights()
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param model character indicating one of the following outcome models:
#'  * "m0" (exposure main effects)
#'  * "m1" (exposure main effects & covariates)
#'  * "m2" (exposure main effects & their interactions)
#'  * "m3" (exposure main effects, their interactions, & covariates)
#' @param family (optional) family function specification for svyglm model
#' @param link (optional) character link function specification for svyglm model
#' @param int_order integer specification of highest order exposure main effects
#'   interaction, required for interaction models ("m2", "m3")
#' @param covariates list of characters reflecting variable names of covariates,
#'   required for covariate models ("m1", "m3")
#' @param epochs (optional) data frame of exposure epoch labels and values
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return list of svyglm model output
#' @export
#' @examples
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
#'
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
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' m <- fitModel(data = test,
#'               weights = w,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               model = "m0",
#'               save.out = FALSE)
#' m <- fitModel(data = test,
#'               weights = w,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               model = "m0",
#'               family = gaussian,
#'               link = "identity",
#'               epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
#'                                   values = I(list(c(1, 2), c(3)))),
#'               save.out = FALSE)
#' m <- fitModel(data = test,
#'               weights = w,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               model = "m1",
#'               covariates = "C",
#'               save.out = FALSE)
#' m <- fitModel(data = test,
#'               weights = w,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               model = "m2",
#'               int_order = 3,
#'               save.out = FALSE)
#' m <- fitModel(data = test,
#'               weights = w,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               model = "m3",
#'               int_order = 3,
#'               covariates = "C",
#'               save.out = FALSE)

fitModel <- function(data, weights, exposure, exposure_time_pts, outcome, model,
                     family = NULL, link = NA, int_order = NA, covariates = NULL, epochs = NULL,
                     home_dir = NULL, verbose = TRUE, save.out = TRUE) {
  
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
  else if (!is.character(outcome) || length(outcome) != 1) {
    stop ("Please supply a single outcome as a character.",
          call. = FALSE)
  }
  else if (!grepl("\\.", outcome)) {
    stop ("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
          call. = FALSE)
  }
  
  if (missing(weights)) {
    stop("Please supply a list of IPTW weights.",
          call. = FALSE)
  }
  if (!is.list(weights) || is.data.frame(weights)) {
    stop("Please supply a list of weights output from the createWeights function.",
          call. = FALSE)
  }
  if (is.list(weights) && !is.data.frame(weights) &&
      !all(sapply(weights, inherits, "weightitMSM"))) {
      stop("Please supply a list of weights output from the createWeights function.",
            call. = FALSE)
  }
  
  if (missing(exposure_time_pts)) {
    stop("Please supply the exposure time points at which you wish to create weights.",
          call. = FALSE)
  }
  if (!is.numeric(exposure_time_pts)) {
    stop("Please supply a list of exposure time points as integers.",
          call. = FALSE)
  }
  else if (!length(exposure_time_pts) > 1) {
    stop ("Please supply at least two exposure time points.",
          call. = FALSE)
  }
  
  if (missing(model)) {
    stop('Please provide an outcome model selection "m" from 0-3 (e.g., "m1")',
          call. = FALSE)
  }
  
  if (!is.character(model)) {
    stop('Please provide as a character string a valid model "m" from 0-3 (e.g., "m1")',
          call. = FALSE)
  }
  if (!is.character(model) || length(model) != 1) {
    stop('Please provide a single outcome model selection "m" from 0-3 (e.g., "m1")',
          call. = FALSE)
  }
  if (!(model %in% c("m0", "m1", "m2", "m3"))) {
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")',
          call. = FALSE)
  }
  if ((model %in% c("m2", "m3")) && (is.na(int_order) || 
                                     !is.numeric(int_order) || 
                                     length(int_order) > 1)) {
    stop("Please provide an integer interaction order if you select a model with interactions.",
         call. = FALSE)
  }
  if ((model %in% c("m1", "m3")) && (is.null(covariates) || 
                                          !is.character(covariates))) {
    stop("Please provide a list of covariates as characters if you select a covariate model.",
          call. = FALSE)
  }
  
  if (is.null(family)) {
    family <- gaussian
  }
  if (!is.function(family)) {
    stop("Please provide a valid family in the form of a function (without quotations).",
          call. = FALSE)
  }
  
  if (anyNA(link)) {
    link <- "identity"
  }
  if (!is.character(link) || length(link) != 1) {
    stop("Please provide as a character a valid link function.",
          call. = FALSE)
  }

  if (!is.null(covariates)) {
    if (!is.character(covariates)) {
      stop("Please provide a list of character strings for covariates.",
            call. = FALSE)
    }
    if (sum(as.numeric(sapply(strsplit(covariates, "\\."), "[", 2)) > 
            exposure_time_pts[1], na.rm = T) > 0) {
      warning("Please only include covariates that are time invariant or measured at the first exposure time point.",
               call. = FALSE)
    }
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
  
  if (verbose || save.out) {
    rlang::check_installed("sjPlot")
  }
  
  if (inherits(data, "mids") || (is.list(data) && !is.data.frame(data))) {
    rlang::check_installed(c("mice", "mitml"))
  }
  
  if (save.out) {
    models_dir <- file.path(home_dir, "models")
    if (!dir.exists(models_dir)) {
      dir.create(models_dir)
    }
  }
  
  # Lists out exposure-epoch combos
  
  if (is.null(epochs)) { #making epochs time pts if not specified by user
    
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
    
  }
  else {
    
    if (!is.data.frame(epochs) || ncol(epochs) != 2 ||
         !all(colnames(epochs) == c("epochs", "values"))) {
      stop("If you supply epochs, please provide a dataframe with two columns of epochs and values.",
            call. = FALSE)
    }
    if (anyNA(epochs$values)) {
      stop("Please provide one or a list of several values for each epoch.",
            call. = FALSE)
    }
  }
  
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, 
                      paste, sep = "", collapse = ".")
  
  #getting null comparisons for LHT
  
  n <- {
    if (model %in% c("m0", "m2")) "int"
    else "covs"
  }
  
  l <- link
  family <- family(link = l)

  if (inherits(data, "mids")) { #imputed dataset
    
    fits <- lapply(seq_len(data$m), function(y) {
      
      d <- mice::complete(data, y)
      d$weights <- weights[[y]]$weights
      
      getModel(d = d, exposure = exposure, exposure_time_pts = exposure_time_pts, 
               outcome = outcome, epochs = epochs, exp_epochs = exp_epochs, 
               int_order = int_order, model = model, fam = family, 
               covariates = covariates, verbose = verbose)
      
    })
    
    fits.null <- lapply(seq_len(data$m), function(y) {
      
      d <- mice::complete(data, y)
      d$weights <- weights[[y]]$weights
      
      getModel(d = d, exposure = exposure, exposure_time_pts = exposure_time_pts, 
               outcome = outcome, epochs = epochs, exp_epochs = exp_epochs, 
               int_order = int_order, model = n, fam = family, 
               covariates = covariates, verbose = verbose)
      
    })
    
    if (verbose) {
      message("Please inspect the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n\n")
      message("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n\n")
    }
    
    print(mitml::testModels(fits, fits.null))
    cat("\n")
  }
  else if (is.list(data) && !is.data.frame(data)) { #imputed dataset
    fits <- lapply(seq_len(length(data)), function(y) {
      
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights
      getModel(d = d, exposure = exposure, exposure_time_pts = exposure_time_pts, 
               outcome = outcome, epochs = epochs, exp_epochs = exp_epochs, 
               int_order = int_order, model = model, fam = family, 
               covariates = covariates, verbose = verbose)
      
    })
    
    fits.null <- lapply(seq_len(length(data)), function(y) {
      
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights
      getModel(d = d, exposure = exposure, exposure_time_pts = exposure_time_pts, 
               outcome = outcome, epochs = epochs, exp_epochs = exp_epochs, 
               int_order = int_order, model = n, fam = family, 
               covariates = covariates, verbose = verbose)
      
    } )
    
    if (verbose) {
      message("Please inspect the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n\n")
      message("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n\n")
    }
    
    print(mitml::testModels(fits, fits.null))
    cat("\n")
  }
  else { #df
    d <- data
    d$weights <- weights[["0"]]$weights
    
    fits <- list(getModel(d = d, exposure = exposure, exposure_time_pts = exposure_time_pts, 
                          outcome = outcome, epochs = epochs, exp_epochs = exp_epochs, 
                          int_order = int_order, model = model, fam = family, 
                          covariates = covariates, verbose = verbose))
    
    fits.null <- list(getModel(d = d, exposure = exposure, exposure_time_pts = exposure_time_pts, 
                               outcome = outcome, epochs = epochs, exp_epochs = exp_epochs, 
                               int_order = int_order, model = n, fam = family, 
                               covariates = covariates, verbose = verbose))

    if (verbose) {
      message("Please inspect the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n\n")
      message("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n\n")
    }
    
    print(anova(fits[[1]], fits.null[[1]]))
    names(fits) <- "0"
    cat("\n")
  }
  
  
  if (inherits(data, "mids") || ((is.list(data)) && !is.data.frame(data))) {
    
    names(fits) <- seq_len(length(fits))
    
    if (verbose) {
      sprintf("The marginal model, %s run for each imputed dataset is summarized below: \n",
              model)
      
      
      cat(sprintf("The marginal model, %s, is summarized below:\n", model))
      
      print(sjPlot::tab_model(fits, auto.label = FALSE, show.se = TRUE))
      
    }
    
    if (save.out) {
      
      print(sjPlot::tab_model(fits, auto.label = TRUE, show.se = TRUE, 
                              dv.labels = paste("Imp", c(1:length(data)), sep = " "),
                              file = file.path(home_dir, "models",
                                               sprintf("%s-%s_%s_table_mod_ev.doc",
                                                       exposure, outcome, model))))
    }
    
    names(fits) <- NULL
  }
  else {
    if (verbose) {
      
      cat(sprintf("The marginal model, %s, is summarized below:\n", model))
      
      print(sjPlot::tab_model(fits, auto.label = FALSE, show.se = TRUE))
      
    }
    
    if (save.out) {
      
      print(sjPlot::tab_model(fits, auto.label = FALSE, show.se = TRUE,
                              file = file.path(home_dir, "models",
                                               sprintf("%s-%s_%s_table_mod_ev.doc",
                                                       exposure, outcome, model))))
    }
    
  }
  
  if (save.out) {
    saveRDS(fits,
            file = file.path(home_dir, "models",
                             sprintf("%s-%s_%s_model.rds",
                                     exposure, outcome, model)))
    cat("\n")
    
    if (verbose) {
      cat("\n")
      print(lapply(fits, summary))
      cat("\n")
      cat("Tables of model evidence have now been saved in the 'models/' folder.\n")
    }
  }
  
  fits
  
}
