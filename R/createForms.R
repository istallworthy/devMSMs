# Creates formulas for generating weights
#'
#' Creates 3 types of formulas for each exposure at each time point including all time-varying covariates at all prior lags
#'
#' @param object msm object that contains all relevant user inputs
#' @param wide_long_datasets from formatForWeights
#' @param all_potential_covariates from selectCoariates
#' @return forms
#' @export
#' @seealso [formatForWeights()], [identifyPotentialConfounds] for more on inputs
#' @examples createForms(object, wide_long_datasets, all_potential_covariates)

createForms <- function(object, all_potential_covariates) {
  # Avoid unnecessary assignments, use list destructuring instead
  # Extract the required elements from the object directly
  with(object, {
    # Use meaningful variable names and avoid excessive comments
    exposure_var <- exposure
    outcome_var <- outcome
    exposure_time_pts <- exposure_time_pts
    time_varying_covariates <- time_varying_variables
    keep_concurrent_tv_vars <- keep_concurrent_tv_vars
    home_dir <- home_dir
    ID <- ID

    # Move the user alert message to a separate function or use a message/warning
    message("USER ALERT: Please inspect the full balancing formulas for each exposure time point below. They contain all confounders that will be used to assess balance at each exposure time point.",
            "By default, these formulas exclude any contemporaneous time-varying variables (given that they are difficult to disentangle from mediators).",
            "If there are contemporaneous time-varying variables that you wish to include in the balancing formula, please list them in the 'keep_concurrent_tv_vars' field of the msmObject and re-run this function.")

    cat("\n")

    # Create a new folder for the forms if it doesn't exist
    forms_dir <- file.path(home_dir, "forms")
    if (!dir.exists(forms_dir)) {
      dir.create(forms_dir)
    }

    # Use vector instead of data frame for forms_csv
    forms_csv <- character()

    # Cycles through all exposure time points
    for (x in seq_along(exposure_time_pts)) {
      time <- exposure_time_pts[x]

      # Find variables to include for that time point
      vars_to_include <- all_potential_covariates[
        as.numeric(sapply(strsplit(all_potential_covariates, "\\."), "[", 2)) < time |
          is.na(sapply(strsplit(all_potential_covariates, "\\."), "[", 2))
      ]

      if (x > 1) {
        # Adds lagged exposure
        vars_to_include <- c(
          vars_to_include,
          apply(expand.grid(exposure_var, as.character(exposure_time_pts[exposure_time_pts < time])),
                1, paste, sep = "", collapse = ".")
        )
      }

      vars_to_include <- vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                                 time_varying_covariates,
                                                                 paste(exposure_var, time, sep = "."),
                                                                 apply(expand.grid(time_varying_covariates[!time_varying_covariates
                                                                                                           %in% keep_concurrent_tv_vars],
                                                                                   as.character(time)), 1, paste, sep = "", collapse = "."),
                                                                 paste(outcome_var, time, sep = "."),
                                                                 apply(expand.grid(outcome_var, as.character(time_pts[time_pts > time])),
                                                                       1, paste, sep = "", collapse = "."),
                                                                 time_var_exclude)]

      vars_to_include <- vars_to_include[!duplicated(vars_to_include)]

      # Creates form for the given exposure time point
      f <- as.formula(paste(paste0(exposure_var, ".", time, " ~ "),
                            paste0(vars_to_include[order(vars_to_include)], sep = "", collapse = " + ")))

      # Prints form for user inspection
      cat(paste0("The full formula for ", exposure_var, "-", outcome_var, " at ", exposure_var, " time point ",
                 as.character(time), " is:"), "\n")
      print(f)
      cat("\n")

      # Appends the form string to forms_csv
      forms_csv <- c(forms_csv, paste0("Full formula for ", exposure_var, "-", outcome_var, " at ", exposure_var,
                                       " time point ", as.character(time), ":"))
      forms_csv <- c(forms_csv, paste(exposure_var, "~", paste0(vars_to_include[order(vars_to_include)], sep = "", collapse = " + ")))

      # Assigns the form to forms list
      assign(paste("form_", exposure_var, "-", outcome_var, "-", time, sep = ""), f, envir = parent.frame())
    }

    # Writes forms_csv to a CSV file
    forms_csv_file <- file.path(forms_dir, paste0(exposure_var, "-", outcome_var, "_full_balancing_formulas.csv"))
    writeLines(forms_csv, con = forms_csv_file)

    # Return the list of forms
    return(ls(pattern = "form_", envir = parent.frame()))
  })
}
