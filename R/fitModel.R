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

fitModel <- function(home_dir, data, weights, exposure, outcome, tv_confounders, model, family, link, int_order = NA< covariates = NULL, epochs = NULL) {

  weights_method <- weights[[1]]$method
  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))

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
    imps <- lapply(files, function(file) {
      imp_data <- read.csv(file)
      imp_data
    })

    # Process imputed datasets
    imps2 <- lapply(1:length(imps), function(x) {
      imp_data <- imps[[x]]
      v <- sapply(strsplit(tv_confounders, "\\."), "[",1)
      v <- v[!duplicated(v)]
      library(splitstackshape)
      imp_data <- merged.stack(
        imp_data,
        id.vars = c(ID, ti_confounders),
        var.stubs = v,
        sep = ".",
        keep.all = TRUE
      )
      colnames(imp_data)[colnames(imp_data) == ".time_1"] <- "WAVE"
      imp_data$.imp <- x - 1
      imp_data <- data.frame(imp_data)
      imp_data
    })
    # Combine imputed datasets into mids object
    imp2 <- do.call(rbind.data.frame, imps2)
    imp2$X <- seq_len(nrow(imp2))
    data <- mice::as.mids(imp2, .imp = ".imp")
  }


  # Listing any baseline imbalanced covariates
  # Renames factors (that were appended w/ level)


  # Covariate models checking
  if (model == "m1" | model == "m3") {
    covariate_list <- paste(as.character(covariates), sep = "", collapse = " + ")
  }

  # Lists out exposure-epoch combos
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(time_pts),
                         values = time_pts)
  }
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, paste, sep = "", collapse = "_")

  # Getting interactions between exposure epoch main effects (e.g., infancy:toddlerhood)
  if (model == "m2" | model == "m3"){
    interactions <- paste(
      lapply(2:length(exp_epochs), function(z) {
        paste(apply(combn(exp_epochs, z), 2, paste, sep = "", collapse = ":"), sep = "", collapse = " + ")
      }),
      collapse = " + "
    )
  }

  if (class == "list"){ #imputed dataset
    fits <- lapply(1:length(data), function(y) {
      d <- data[[y]]
      d$weights <- NULL
      d$weights <- weights[[y]]$weights

      # Getting design info
      s <- survey::svydesign(
        id = ~1, #
        data = d, # Adds list of imputation data?
        weights = ~weights
      )

      # Fitting baseline model w/ main effects only (m0) for all models
      f0 <- paste(paste0(outcome, ".", outcome_time_pt), "~", paste0(exp_epochs, sep = "", collapse = " + "))
      f0 <- paste(paste0(outcome, ".", outcome_time_pt), "~", paste0(exp_epochs, sep = "", collapse = " + "))

      m0 <- survey::svyglm(as.formula(f0), design = s) # List of model fitted to all imputed datasets

      if (model == "m0") {
        return(m0) # Save model
      } else {

        # Baseline + sig covar model OR baseline + sig covar + int model
        if (model == "m1" | model == "m3") {
          # Fitting m1
          f1 <- paste(f0, "+", covariate_list) # Baseline + covariate model
          m1 <- survey::svyglm(as.formula(f1), design = s)
          # Baseline + imbalanced covars
          if (model == "m1") {
            return(m1)
          }
        }

        # Baseline + interaction OR baseline + covars + interactions
        if (model == "m2" | model == "m3") {
          f2 <- paste(f0, "+", paste(interactions, sep = "", collapse = " + "))

          # Baseline + interactions
          if (model == "m2") {
            # Fitting m2
            m2 <- survey::svyglm(as.formula(f2), design = s)
            # Baseline + interactions
            return(m2)
          }

          # Baseline + covars + interactions
          if (model == "m3") {
            # Fitting m3
            f3 <- paste(f1, "+", paste(interactions, sep = "", collapse = " + "))
            m3 <- survey::svyglm(as.formula(f3), design = s)
            # Baseline + covars+ interactions
            return(m3)
          }
        }
      }
    })
  })
names(fits) <- all_cutoffs

cat(paste0("USER ALERT: the marginal model, ", model, ", run for each imputed dataset using the user-specified weights truncation percentile value of ",
           paste(weights_percentile_cutoff), " is summarized below for each imputed dataset:\n"))
print(lapply(fits[paste0(weights_percentile_cutoff)], function(x) lapply(x, summary)))

# Print table of model evidence comparing imputations per cutoff level
lapply(1:length(fits), function(y) {
  i <- fits[[y]]
  suppressWarnings(jtools::export_summs(
    i, to.file = "docx", statistics = c(N = "nobs", AIC = "AIC", R2 = "r.squared"),
    model.names = c(paste0("Imp.", 1:length(i))),
    file.name = file.path(home_dir, "msm", paste0(exposure, "-", outcome, "_", names(fits)[y], "_", model, "_table_mod_ev.docx"))
  ))
})

cat("Tables of model evidence for all truncation percentile values have now been saved in the 'msm' folder.\n")

fits <- NULL

saveRDS(fits, file = file.path(home_dir, "msm", paste0(exposure, "-", outcome, "_", model, "_model.rds")))
cat("\n")

return(fits)

}
