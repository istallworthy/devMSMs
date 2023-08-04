#' Create short forms that contain time invariant covariates and time varying covariates at only t-1 or user-spec lag
#'
#' @param object msm object that contains all relevant user inputs
#' @param full_forms list of full forms from createForms
#' @param keep optional list specifying any time-varying covariates at larger lags to keep in these forms
#' @return short_forms
#' @export
#' @examples createShortForms(object, full_forms, keep=NULL)
#'
createShortForms <- function(object, full_forms, keep = NULL) {
  # Validate inputs
  # if (!inherits(object, "msm")) {
  #   stop("The 'object' parameter must be an msm object.")
  # }

  if (!is.list(full_forms)) {
    stop("The 'full_forms' parameter must be a list.")
  }

  if (!is.null(keep) && !is.list(keep)) {
    stop("The 'keep' parameter must be a list.")
  }

  # Use list destructuring to extract relevant elements from the object
  with(object, {
    home_dir <- home_dir
    exp_time_pts <- exposure_time_pts
    short_form_lag <- short_form_lag
    exposure_var <- exposure
    outcome_var <- outcome
  })

  # Print informative message to the user
  message("Short formulas containing all time invariant covariates and time-varying covariates only at t-", short_form_lag,
          " will be created at each of the following exposure time points: ",
          paste(exp_time_pts[short_form_lag + 2:length(exp_time_pts) - 1], collapse = ", "), ".")

  short_forms <- list()
  forms_csv <- data.frame()

  # Getting relevant forms
  list <- full_forms[names(full_forms)[grepl(exposure_var, names(full_forms))]]

  # Error checking
  if (length(list) != length(exp_time_pts)) {
    stop("Make sure the forms list contains formulas for each exposure time point for each exposure.")
  }

  # Cycling through all full formulas (at each exposure time point)
  for (z in 1:length(list)) {
    time <- as.numeric(sapply(strsplit(names(list)[z], "-"), "[", length(unlist(strsplit(names(list)[z], "-")))))

    # Find any user-specified covariates to keep
    keep_cov <- NULL
    if (!is.null(keep)) {
      matching_keep <- keep$exposure == exposure_var & keep$time == time
      if (any(matching_keep)) {
        keep_cov <- keep$tv_covar[matching_keep]
      }
    }

    if (time == exp_time_pts[1] || time == exp_time_pts[2]) { # Ignore first time point as there are no lagged values and second time point because only t-1 exists
      short_forms[[names(list)[z]]] <- list[[z]] # Populate with full formula
    } else {
      form <- list[[z]]
      dv <- form[[2]]
      covars <- paste(deparse(form[[3]], width.cutoff = 500), collapse = "")
      covar_time <- sapply(strsplit(unlist(strsplit(as.character(covars), "\\+")), "\\."), "[", 2)
      covars <- as.character(unlist(strsplit(covars, "\\+")))

      if (!is.null(keep_cov)) {
        new_covars <- c(covars[!as.numeric(covar_time) < exp_time_pts[z - 1] | is.na(covar_time)], keep_cov) # Add any user-specified covars to keep
      } else {
        new_covars <- covars[!as.numeric(covar_time) < exp_time_pts[z - short_form_lag] | is.na(covar_time)] # Filter out longer lags
      }

      new_form <- as.formula(paste0(dv, "~", paste(new_covars, sep = "", collapse = "+")))

      message(paste0("The short formula for ", names(list)[z], " including time-varying covariates at t-", short_form_lag, " only is:"))
      print(new_form)

      # Write to csv file
      forms_csv_temp <- data.frame()
      forms_csv_temp[1, 1] <- paste0("Short formula for ", exposure_var, "-", outcome_var, " at ", exposure_var, " time point ", as.character(time), ":")
      forms_csv_temp[1, 2] <- paste0(dv, "~", paste(new_covars, sep = "", collapse = "+"))
      forms_csv <- rbind(forms_csv, forms_csv_temp)

      short_forms[[names(list)[z]]] <- new_form
    }
  }

  # Save to files
  forms_dir <- file.path(home_dir, "forms")
  if (!dir.exists(forms_dir)) {
    dir.create(forms_dir)
  }
  forms_csv_file <- file.path(forms_dir, paste0(exposure_var, "-", outcome_var, "_short_balancing_formulas.csv"))
  write.csv(forms_csv, file = forms_csv_file, row.names = FALSE)

  short_forms_file <- file.path(forms_dir, paste0(exposure_var, "-", outcome_var, "_short_forms.rds"))
  saveRDS(short_forms, file = short_forms_file)

  message(paste0("Across all short balancing formulas at all exposure time points, there are a total of ",
                 sum(unlist(lapply(1:nrow(forms_csv), function(x) { length(unlist(strsplit(forms_csv[x, 2], "\\+"))) }))),
                 " covariate confounders."))

  message("Short formulas including time-varying covariates at t-", short_form_lag, " only have now been saved in the 'forms/' folder.")

  return(short_forms)
}
