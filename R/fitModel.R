#' Fit weighted model
#' This code fits a weighted marginal structural model to examine the effects of different exposure histories on outcome
#' @param model user-specified model
#' @return fits
#' @export
#' @importFrom survey svydesign svyglm
#' @importFrom jtools export_summs
#' @importFrom dplyr mutate filter select
#' @return fits
#' @export

fitModel <- function(home_dir, data, weights, exposure, exposure_time_pts, outcome, tv_confounders, model, family = gaussian, link = "identity", int_order = NA, covariates = NULL, epochs = NULL, user.o = TRUE) {

  if (!(model %in% c("m0", "m1", "m2", "m3"))) {
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")')
  }
  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.")
  }
  if (!class(data) %in% c("mids", "data.frame", "character")) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.")
  }
  if (!(model %in% c("m0", "m1", "m2", "m3"))) {
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")')
  }
  if ((model == "m2" | model == "m3") & (is.na(int_order))){
    stop("Please provide an interaction order if you select a model with interactions.")
  }
  if ((model == "m1" | model == "m3") & (is.null(covariates))){
    stop("Please provide a list of covariates if you select a covariate model.")
  }
  if (!is.null(covariates) & sum(as.numeric(sapply(strsplit(covariates, "\\."), "[",2)) > exposure_time_pts[1], na.rm = T) >0){
    warning("Please only include covariates that are time invariant or measured at the first exposure time point.")
  }

  weights_method <- weights[[1]]$method

  #create models directory
  models_dir <- file.path(home_dir, "models")
  if (!dir.exists(models_dir)) {
    dir.create(models_dir)
  }


  if (class(data) == "character") {
    if (!dir.exists(data)) {
      stop("Please provide a valid directory path with imputed datasets, a data frame, or a 'mids' object for the 'data' field.")
    }
    if (length(dir(data)) < 2) {
      stop("If you specify data as a directory, please supply more than 1 imputed dataset.")
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
  if (model == "m1") {n <- "covs"}
  if (model == "m2") {n <- "int"}
  if (model == "m3") {n <- "covs"}

  # if (family == "gaussian"){
    l <- link
    family <- family(link = l)
    assign("fam", family, envir = .GlobalEnv)
    # family <- family(link = l)
    # fam <- family(link = l)
  # }

  if (class(data) == "mids"){ #imputed dataset
    fits <- lapply(1:data$m, function(y) {
      d <- complete(data, y)
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model, family, covariates, interactions, user.o)
    })

    fits.null <- lapply(1:data$m, function(y) {
      d <- complete(data, y)
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model = n, family, covariates, interactions, user.o)
    })

    if (user.o == TRUE){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    print(mitml::testModels(fits, fits.null))
    cat("\n")

  }


  if (class(data) == "list"){ #imputed dataset
    fits <- lapply(1:length(data), function(y) {
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model, family, covariates, interactions, user.o)
    })

    fits.null <- lapply(1:length(data), function(y) {
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model = n, family, covariates, interactions, user.o)
    })

    if (user.o == TRUE){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    print(mitml::testModels(fits, fits.null))
    cat("\n")

  }

# fam = family
  if (class(data) == "data.frame"){ #df
    fits <- lapply(1, function(y) {
      d <- data
      d$weights <- NULL
      d$weights <- weights[["0"]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model, family, covariates, interactions, user.o)
    })

    fits.null <- lapply(1, function(y) {
      d <- data
      d$weights <- NULL
      d$weights <- weights[["0"]]$weights
      getModel(d, outcome, epochs, exp_epochs, int_order, model = n, family, covariates, interactions, user.o)
    })

    if (user.o == TRUE){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    print(anova(fits[[1]], fits.null[[1]]))
    names(fits) = "0"
    cat("\n")

  }



  if (class(data) == "mids" | class(data) == "list"){
    if (user.o == TRUE){
      cat(paste0("USER ALERT: the marginal model, ", model, ", run for each imputed dataset is summarized below:"), "\n")
    }

    suppressWarnings(jtools::export_summs(
      fits, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
      model.names = c(paste0("Imp.", 1:length(fits))),
      file.name = file.path(home_dir, "models", paste0(exposure, "-", outcome, "_", model, "_table_mod_ev.docx"))
    ))

  } else{
    if (user.o == TRUE){
      cat(paste0("USER ALERT: the marginal model, ", model, ", is summarized below:"), "\n")
    }

    suppressWarnings(jtools::export_summs(fits, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
      file.name = file.path(home_dir, "models", paste0(exposure, "-", outcome, "_", model, "_table_mod_ev.docx"))
    ))


  }

  names(fits) <- 1:length(fits)

  if (user.o == TRUE){
    cat("\n")
    print(lapply(fits, function(x) {summary(x)}))
    cat("\n")
    cat("Tables of model evidence have now been saved in the 'models/' folder.\n")
  }

  names(fits) <- NULL

  saveRDS(fits, file = file.path(home_dir, "/models/", paste0(exposure, "-", outcome, "_", model, "_model.rds")))
  cat("\n")

  fits

}
