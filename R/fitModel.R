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

fitModel <- function(msm_object, data_for_model_with_weights_cutoff, balance_stats_final, model = "m3") {

  home_dir <- msm_object$home_dir
  exposure <- msm_object$exposure
  exposure_epochs <- msm_object$exposure_epochs
  outcome <- msm_object$outcome
  outcome_time_pt <- msm_object$outcome_time_pt
  factor_covariates <- msm_object$factor_covariates
  weights_percentile_cutoff <- msm_object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity <- as.numeric(unlist(strsplit(msm_object$weights_percentile_cutoffs_sensitivity, " ")))
  m <- msm_object$m
  balance_thresh <- msm_object$balance_thresh
  time_varying_covariates <- msm_object$time_varying_variables
  time_pts <- msm_object$time_pts
  weights_method <- msm_object$weights_method

  # Error checking
  if (!(model %in% c("m0", "m1", "m2", "m3"))) {
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")')
  }

  cat(paste0("Fitting model ", model, " with weights generated using the ", weights_method, " to each of the ", m, " imputed datasets.\n"))
  cat("\n")

  all_cutoffs <- c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  # Averages across all imputed dataset bal stats to determine imbalanced covariates (baseline ones will be used in any covariate models)
  unbalanced_covars <- as.data.frame(rowMeans(do.call(cbind, lapply(balance_stats_final, `[[`, "std_bal_stats"))))
  unbalanced_covars <- data.frame(
    exposure = exposure,
    exp_time = balance_stats_final[[1]]$exp_time,
    covar_time = balance_stats_final[[1]]$covar_time,
    covariate = balance_stats_final[[1]]$covariate,
    avg_bal = unname(unbalanced_covars)
  ) %>%
    mutate(balanced_avg = ifelse(abs(avg_bal) < balance_thresh, 1, 0)) %>%
    filter(balanced_avg == 0) %>%
    select(covariate)
  # Getting only time invariant or t1 imbalanced covariates
  unbalanced_covars$keep <- ifelse(grepl(paste0(".", time_pts[1]), unlist(unbalanced_covars$covariate)), 1, 0)
  unbalanced_covars$keep <- ifelse(
    unbalanced_covars$keep == 0 &
      !gsub(" ", "", sapply(strsplit(unlist(unbalanced_covars$covariate), "\\."), "[", 1)) %in% time_varying_covariates,
    1,
    unbalanced_covars$keep
  )
  baseline_covars <- unlist(unbalanced_covars %>% filter(keep == 1) %>% select(covariate))
  nonb_covars <- unlist(unbalanced_covars %>% filter(keep == 0) %>% select(covariate))

  # Listing any baseline imbalanced covariates
  # Renames factors (that were appended w/ level)
  baseline_covars[sapply(strsplit(sapply(strsplit(baseline_covars, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates] <-
    sapply(strsplit(baseline_covars, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(baseline_covars, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates]
  baseline_covars <- baseline_covars[!grepl(exposure, baseline_covars)] # Excludes exposure
  covariate_list <- paste(as.character(baseline_covars), sep = "", collapse = " + ")

  # Covariate models checking
  if (model == "m1" | model == "m3") {
    # If there are no imbalanced covariates
    if (covariate_list[1] == "") {
      stop("You have selected a covariate model but there are no imbalanced baseline covariates to include. Please choose another model.")
    } else { # There are imbalanced covariates
      cat("USER ALERT: The following imbalanced baseline covariates are included in the model: \n")
      cat(covariate_list, "\n\n")
      cat("USER ALERT: The following imbalanced covariates will NOT be included in the model: \n")
      cat(nonb_covars, "\n\n")
    }
  }

  # Lists out exposure-epoch combos
  exp_epochs <- apply(expand.grid(exposure, as.character(exposure_epochs[, 1])), 1, paste, sep = "", collapse = "_")

  # Getting interactions between exposure epoch main effects (e.g., infancy:toddlerhood)
  interactions <- paste(
    lapply(2:length(exp_epochs), function(z) {
      paste(apply(combn(exp_epochs, z), 2, paste, sep = "", collapse = ":"), sep = "", collapse = " + ")
    }),
    collapse = " + "
  )

  fits <- lapply(1:length(data_for_model_with_weights_cutoff), function(y) { # Cycling through truncation cutoff values
    d <- data_for_model_with_weights_cutoff[[y]]
    cutoff <- all_cutoffs[y]
    cutoff_label <- ifelse(cutoff == weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))

    lapply(d, function(x) { # Cycling through imputed datasets
      data <- x
      data$weights <- NULL
      data$weights <- data[, colnames(data)[grepl("weight", colnames(data))]]

      # Getting design info
      s <- survey::svydesign(
        id = ~1, #
        data = data, # Adds list of imputation data
        weights = ~weights
      )

      # Fitting m0
      f0 <- paste(paste0(outcome, ".", outcome_time_pt), "~", paste0(exp_epochs, sep = "", collapse = " + "))
      m0 <- survey::svyglm(as.formula(f0), design = s) # List of model fitted to all imputed datasets
      # Baseline model (main effects) only
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
