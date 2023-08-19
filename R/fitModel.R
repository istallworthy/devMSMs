#' Fit weighted model
#' This code fits a weighted marginal structural model to examine the effects of different exposure histories on outcome
#' @param msm_object msm object that contains all relevant user inputs
#' @param data_for_model_with_weights_cutoff dataset with truncated weights see truncateWeights
#' @param balance_stats_final final bal stats conducted w/ full forms
#' @param model user-specified model
#' @return fits
#' @export
#' @importFrom survey svydesign svyglm
#' @importFrom jtools export_summs
#' @importFrom dplyr mutate filter select
#' @return fits
#' @export
#' @seealso [truncateWeights()], [asesssBalance()]
#' @examples fitModel(object, data_for_model_with_weights_cutoff, balance_stats_final, model="m3")


fitModel <- function(home_dir, data, weights, exposure, exposure_time_pts, outcome, tv_confounders, model, family = gaussian, link = "identity", int_order = NA, covariates = NULL, epochs = NULL, user.o = TRUE) {

  # Error checking
  if (!(model %in% c("m0", "m1", "m2", "m3"))) {
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")')
  }

  # exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))

  #error checking
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




  getModel <- function(d, outcome, epochs, exp_epochs, int_order, model, fam, link, covariates, interactions){
    if (sum(duplicated(d$"ID")) > 0){
      stop("Please provide wide data with a single row per ID.")
    }

    # exposure epochs
    if( is.null(epochs)){ #making epochs time pts if not specified by user
      epochs <- data.frame(epochs = as.character(exposure_time_pts),
                           values = exposure_time_pts)
    } else { #add epochs by averaging exposure time points
      #adds exposure epochs

      #calculates the mean value for each exposure for each exposure epoch
      for (e in 1:nrow(epochs)){
        epoch <- epochs[e,1]
        temp <- data.frame(row.names = 1:nrow(d))
        new_var <- paste0(exposure, "_", epoch)
        if (! new_var %in% colnames(d)){
          #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
          for (l in 1:length(as.numeric(unlist(epochs[e,2])))){
            level <- as.numeric(unlist(epochs[e,2]))[l]
            z  <- d[,which(grepl(paste0(exposure, ".", level), names(d)))]
            temp <- cbind(temp, z)
          }
          #adds a new variable of the exposure averaged within epoch
          d  <- d %>% dplyr::mutate(!!new_var := rowMeans(temp, na.rm=T))
          d[,new_var] <- as.numeric(d[,new_var])
          # d <- cbind(d, new)
        }
      }

    }

    # Covariate models checking
    if (model == "m1" | model == "m3") {
      if (sum(covariates %in% colnames(d)) < length(covariates)){
        stop("Please only include covariates that correspond to variables in the wide dataset.")
      }

      covariate_list <- paste(as.character(covariates), sep = "", collapse = " + ")
    } else{
      covariate_list <- NULL

      covariate_list <- paste(as.character(covariates), sep = "", collapse = " + ")
    }


    # interaction model checking
    if (model == "m2" | model == "m3"){
      if (int_order > nrow(epochs)){
        stop("Please provide an interaction order equal to or less than the total number of epochs/time points.")
      }
      interactions <- paste(
        lapply(2:int_order, function(z) {
          paste(apply(combn(exp_epochs, z), 2, paste, sep = "", collapse = ":"), sep = "", collapse = " + ")
        }),
        collapse = " + "
      )
      #create interactions in data
      for (x in 1:length(unlist(strsplit(interactions, "\\+")))) {
        name <- gsub(" ", "", unlist(strsplit(interactions, "\\+"))[x])
        if (! name %in% colnames(d)){
          temp <- d[, c(gsub(" ", "", as.character(unlist(strsplit(unlist(strsplit(interactions, "\\+"))[x], "\\:"))))) ]
          d  <- d %>% dplyr::mutate(!! name := matrixStats::rowProds(as.matrix(temp), na.rm=T))
        }
      }
    } else {
      interactions <- NULL
    }


    if (user.o == TRUE){
      cat("Please insepct the following likelihood ratio test to determine if the exposures collective predict significant variation in the outcome compared to a model without exposure terms.", "\n")
      cat("\n")
      cat("We strongly suggest not conducting history comparisons if the likelihood ratio test is non-significant.", "\n")
      cat("\n")
    }

    s <- survey::svydesign(
      id = ~1, #
      data = d, # Adds list of imputation data?
      weights = ~ weights
    )

    # Fitting baseline model w/ main effects only (m0) for all models
    f0 <- paste(outcome, "~", paste0(exp_epochs, sep = "", collapse = " + "))
    m0 <- survey::svyglm(as.formula(f0), family = fam(link = link),  design = s) # List of model fitted to all imputed datasets

    if (model == "m0") {

      null <- survey::svyglm(as.formula(paste(outcome, "~", 1)), family = fam(link = link),  design = s)
      print(anova(null, m0))

      return(m0) # Save model

    } else {

      # Baseline + sig covar model OR baseline + sig covar + int model
      if (model == "m1" | model == "m3") {
        f1 <- paste(f0, "+", covariate_list) # Baseline + covariate model
        m1 <- survey::svyglm(as.formula(f1), family = fam(link = link), design = s)

        # Baseline + imbalanced covars
        if (model == "m1") {

          null <- survey::svyglm(as.formula(paste(outcome, "~ 1 + ", covariate_list)), family = fam(link = link),  design = s)
          print(anova(null, m1))

          return(m1)
        }
      }

      # Baseline + interaction OR baseline + covars + interactions
      if (model == "m2" | model == "m3") {
        f2 <- paste(f0, "+", paste(interactions, sep = "", collapse = " + "))

        # Baseline + interactions
        if (model == "m2") {
          # Fitting m2
          m2 <- survey::svyglm(as.formula(f2), family = fam(link = link), design = s)
          # Baseline + interactions

          null <- survey::svyglm(as.formula(paste(outcome, "~ 1")), family = fam(link = link),  design = s)
          print(anova(null, m2))

          return(m2)
        }

        # Baseline + covars + interactions
        if (model == "m3") {
          # Fitting m3
          f3 <- paste0(f1, "+", paste(interactions, sep = "", collapse = " + "))
          m3 <- survey::svyglm(as.formula(f3), family = fam(link = link), design = s)
          # Baseline + covars+ interactions

          null <- survey::svyglm(as.formula(paste(outcome, "~ 1 + ", covariate_list)), family = fam(link = link),  design = s)
          print(anova(null, m3))

          # stack= complete(data, "long")
          # # colnames(stack)[colnames(stack) == ".imp"] <- "impvar"
          # psfmi::pool_D4(data = stack, fm1 = f3, nimp = 5, fm0 = as.formula(paste(outcome, "~ 1 + ", covariate_list)), impvar = ".imp")


          return(m3)
        }
      }
    }
  }


  if (class(data) == "mids"){ #imputed dataset
    fits <- lapply(1:data$m, function(y) {
      d <- complete(data, y)
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model, fam, link, covariates, interactions)



    })
  }


  if (class(data) == "list"){ #imputed dataset
    fits <- lapply(1:length(data), function(y) {
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model, fam, link, covariates, interactions)

    })
  }


  if (class(data) == "data.frame"){ #imputed dataset
    fits <- lapply(1, function(y) {
      d <- data
      d$weights <- NULL
      d$weights <- weights[["0"]]$weights

      getModel(d, outcome, epochs, exp_epochs, int_order, model, fam, link, covariates, interactions)
    })
  }
  names(fits) = "0"


  if (user.o == TRUE){
    if (class(data) == "mids" | class(data) == "list"){
      cat(paste0("USER ALERT: the marginal model, ", model, ", run for each imputed dataset is summarized below:"), "\n")

      suppressWarnings(jtools::export_summs(
        fits, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
        model.names = c(paste0("Imp.", 1:length(fits))),
        file.name = file.path(home_dir, "/models/", paste0(exposure, "-", outcome, "_", model, "_table_mod_ev.docx"))
      ))

    } else{
      cat(paste0("USER ALERT: the marginal model, ", model, ", is summarized below:"), "\n")

      suppressWarnings(jtools::export_summs(
        fits, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
        file.name = file.path(home_dir, "/models/", paste0(exposure, "-", outcome, "_", model, "_table_mod_ev.docx"))
      ))

    }

    names(fits) <- 1:length(fits)
    print(lapply(fits, function(x) {summary(x)}))
    names(fits) <- NULL

    cat("Tables of model evidence have now been saved in the 'models/' folder.\n")
  }

  saveRDS(fits, file = file.path(home_dir, "/models/", paste0(exposure, "-", outcome, "_", model, "_model.rds")))
  cat("\n")

  fits

}
