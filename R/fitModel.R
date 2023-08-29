#' Fit outcome model
#'
#' Fits weighted marginal outcome model as a generalized linear model of the
#' user's choosing, relating exposure main effects to outcome using IPTW
#' weights.
#' @importFrom survey svydesign svyglm
#' @importFrom jtools export_summs
#' @importFrom dplyr mutate filter select
#' @seealso {[survey::svyglm()] for more on family/link specifications, <url1>}
#' @seealso {[createWeights()], <url1>}
#' @param home_dir path to home directory
#' @param data data in wide format as: a data frame, path to folder of imputed
#'   .csv files, or mids object
#' @param weights list of IPTW weights output from createWeights()
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix
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
#' @examples fitModel(home_dir, data, weights, exposure, exposure_time_pts, outcome, tv_confounders, model = "m0")

#'
#' @return list of svyglm model output

fitModel <- function(home_dir, data, weights, exposure, exposure_time_pts, outcome, tv_confounders, model, family = gaussian, link = "identity", int_order = NA, covariates = NULL, epochs = NULL, verbose = TRUE) {

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
  if (missing(weights)){
    stop("Please supply a list of IPTW weights.", call. = FALSE)
  }
  if (missing(exposure_time_pts)){
    stop("Please supply the exposure time points at which you wish to create weights.", call. = FALSE)
  }
  if (missing(tv_confounders)){
    stop("Please supply a list of time-varying confounders.", call. = FALSE)
  }
  if (missing(model)){
    stop('Please provide an outcome model selection "m" from 0-3 (e.g., "m1")', call. = FALSE)
  }

  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.", call. = FALSE)
  }
  if (!mice::is.mids(data) & !is.data.frame(data) & !is.character(data)) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.", call. = FALSE)
  }

  if (!is.character(model)){
    stop('Please provide as a character string a valid model "m" from 0-3 (e.g., "m1")', call. = FALSE)
  }
  if (!(model %in% c("m0", "m1", "m2", "m3"))) {
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")', call. = FALSE)
  }
  if ((model == "m2" | model == "m3") & (is.na(int_order) | !is.numeric(int_order))){
    stop("Please provide an integer interaction order if you select a model with interactions.", call. = FALSE)
  }
  if ((model == "m1" | model == "m3") & (is.null(covariates) | !is.character(covariates))){
    stop("Please provide a list of covariates as characters if you select a covariate model.", call. = FALSE)
  }

  if(!inherits(family, "function")){
    stop("Please provide a valid family in the form of a function (without quotations).", call. = FALSE)
  }
  if(!inherits(link, "character")){
    stop("Please provide as a character a valid link function.", call. = FALSE)
  }
  if (!is.null(covariates) & sum(as.numeric(sapply(strsplit(covariates, "\\."), "[",2)) > exposure_time_pts[1], na.rm = T) >0){
    warning("Please only include covariates that are time invariant or measured at the first exposure time point.")
  }
  if (!inherits(weights, "list")){
    stop("Please supply a list of weights output from the createWeights function.", call. = FALSE)
  }


  weights_method <- weights[[1]]$method

  #create models directory
  models_dir <- file.path(home_dir, "models")
  if (!dir.exists(models_dir)) {
    dir.create(models_dir)
  }


  if (is.character(data)) {
    if (!dir.exists(data)) {
      stop("Please provide a valid directory path with imputed datasets, a data frame, or a 'mids' object for the 'data' field.", call. = FALSE)
    }
    if (length(dir(data)) < 2) {
      stop("If you specify data as a directory, please supply more than 1 imputed dataset.", call. = FALSE)
    }

    # List imputed files
    files <- list.files(data, full.names = TRUE, pattern = "\\.csv")

    # Read and process imputed datasets
    data <- lapply(files, function(file) {
      imp_data <- read.csv(file)
      imp_data
    })
  }


  # Lists out exposure-epoch combos
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(time_pts),
                         values = time_pts)
  }
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, paste, sep = "", collapse = "_")


  #getting null comparison
  if (model == "m0") {n <- "int"}
  else if (model == "m1") {n <- "covs"}
  else if (model == "m2") {n <- "int"}
  else if (model == "m3") {n <- "covs"}

  # if (family == "gaussian"){
  l <- link
  family <- family(link = l)
  assign("fam", family, envir = .GlobalEnv)
  # family <- family(link = l)
  # fam <- family(link = l)
  # }

  if (mice::is.mids(data)){ #imputed dataset
    fits <- lapply(seq_len(data$m), function(y) {
      d <- complete(data, y)
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model, family, covariates, interactions, verbose)
    })

    fits.null <- lapply(seq_len(data$m), function(y) {
      d <- complete(data, y)
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model = n, family, covariates, interactions, verbose)
    })

    if (verbose){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    print(mitml::testModels(fits, fits.null))
    cat("\n")
  }


  else if (inherits(data, "list")){ #imputed dataset
    fits <- lapply(seq_len(length(data)), function(y) {
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model, family, covariates, interactions, verbose)
    })

    fits.null <- lapply(seq_len(length(data)), function(y) {
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model = n, family, covariates, interactions, verbose)
    })

    if (verbose){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    print(mitml::testModels(fits, fits.null))
    cat("\n")
  }

  # fam = family
  else if (is.data.frame(data)){ #df
    fits <- lapply(1, function(y) {
      d <- data
      d$weights <- NULL
      d$weights <- weights[["0"]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model, family, covariates, interactions, verbose)
    })

    fits.null <- lapply(1, function(y) {
      d <- data
      d$weights <- NULL
      d$weights <- weights[["0"]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model = n, family, covariates, interactions, verbose)
    })

    if (verbose){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    print(anova(fits[[1]], fits.null[[1]]))
    names(fits) <- "0"
    cat("\n")
  }



  if (mice::is.mids(data) | inherits(data, "list")){
    if (verbose){
      cat(paste0("USER ALERT: the marginal model, ", model, ", run for each imputed dataset is summarized below:"), "\n")
    }
    names(fits) <- seq_len(length(fits))

    suppressWarnings(jtools::export_summs(
      fits, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
      model.names = c(paste0("Imp.", seq_len(length(fits)))),
      file.name = file.path(home_dir, "models", paste0(exposure, "-", outcome, "_", model, "_table_mod_ev.docx"))
    ))

    names(fits) <- NULL


  }
  else{
    if (verbose){
      cat(paste0("USER ALERT: the marginal model, ", model, ", is summarized below:"), "\n")
    }

    suppressWarnings(jtools::export_summs(fits, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
                                          file.name = file.path(home_dir, "models", paste0(exposure, "-", outcome, "_", model, "_table_mod_ev.docx"))
    ))

  }


  if (verbose){
    cat("\n")
    print(lapply(fits, function(x) {summary(x)}))
    cat("\n")
    cat("Tables of model evidence have now been saved in the 'models/' folder.\n")
  }


  saveRDS(fits, file = file.path(home_dir, "/models/", paste0(exposure, "-", outcome, "_", model, "_model.rds")))
  cat("\n")

  fits

}
