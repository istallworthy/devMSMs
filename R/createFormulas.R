

<<<<<<< HEAD
<<<<<<< Updated upstream
createFormulas <- function(exposure, outcome, tv_confounders, ti_confounders, type, bal_stats = NULL, concur_conf = NA, keep_conf=NA, ug=F ){
=======
createFormulas <- function(home_dir, exposure, outcome, tv_confounders, ti_confounders, type, bal_stats = NULL, concur_conf = NA, keep_conf= NA, custom = NULL, ug=F ){
>>>>>>> Stashed changes
=======
createFormulas <- function(home_dir, exposure, outcome, tv_confounders, ti_confounders, type, bal_stats = NULL, concur_conf = NA, keep_conf=NA, ug=F ){
>>>>>>> main

  #error checking
  if(! type %in% c("short", "full", "update")){
    stop("Please provide a type from the following list: 'short', 'full', or 'update'")
  }
<<<<<<< Updated upstream
  if (type != "update" & !is.null(bal_stats)){
    stop ("Please only provide balance statistics for the type 'update'.")
  }
  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.")
  }
  #create parent directory
  forms_dir <- file.path(home_dir, "formulas")
  if (!dir.exists(forms_dir)) {
    dir.create(forms_dir)
  }
  # Create type directory
  forms_dir <- file.path(home_dir, "formulas", type)
  if (!dir.exists(forms_dir)) {
    dir.create(forms_dir)
  }

=======

  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.")
  }
>>>>>>> Stashed changes

  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))
  time_varying_covariates <- tv_confounders
  all_covars <- c(tv_confounders, ti_confounders)

  if (!is.null(custom)){
    if (length(cutom) != length(exposure_time_pts) | class(custom) != "list"){
      stop("If you wish to supply custom formulas, please provide a list of formulas for each exposure time point.")
    }

    forms <- custom

  } else{

    if (type != "update" & !is.null(bal_stats)){
      stop ("Please only provide balance statistics for the type 'update'.")
    }

    #create parent directory
    forms_dir <- file.path(home_dir, "formulas")
    if (!dir.exists(forms_dir)) {
      dir.create(forms_dir)
    }
    # Create type directory
    forms_dir <- file.path(home_dir, "formulas", type)
    if (!dir.exists(forms_dir)) {
      dir.create(forms_dir)
    }



    factor_covariates <- colnames(data)[which(sapply(data, class) == "factor")]
    factor_covariates <- factor_covariates[!factor_covariates %in% "ID"]

    forms_csv <- character()
    forms <- list()

    for (x in seq_along(exposure_time_pts)) {
      time <- exposure_time_pts[x]

      #identifying lagged tv confounders relative to time point based on formula type
      if (type == "full"){
        message("Please manually inspect the full balancing formula below:")
        time_var_include <- time_varying_covariates[as.numeric(sapply(strsplit(time_varying_covariates, "\\."), "[", 2)) < time]
      }

      if (type == "short"){
        message("Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:")
        time_var_include <- time_varying_covariates[as.numeric(sapply(strsplit(time_varying_covariates, "\\."), "[", 2)) == exposure_time_pts[x-1]]
      }

      if (type == "update"){
        if(is.null(bal_stats)){
          stop("Please provide balance statistics if you wish to run the update verison of this function")
        }

        message("Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:")

        time_var_include <- time_varying_covariates[as.numeric(sapply(strsplit(time_varying_covariates, "\\."), "[", 2)) == exposure_time_pts[x-1]]

        new <- bal_stats %>%
          dplyr::mutate(balanced_avg = ifelse(abs(avg_bal) < balance_thresh, 1, 0)) %>%
          dplyr::filter(balanced_avg == 0)%>%
          dplyr::filter(exp_time == time, as.numeric(covar_time) < exposure_time_pts[x - 1],
                        as.numeric(covar_time) > 0) %>% # Finds any lagged imbalanced covars
          dplyr::select(covariate)

        # Renames factors (that were appended w/ level)
        if (nrow(new) > 0) {
          new$covariate[sapply(strsplit(new$covariate, "_"), `[`, 1) %in% factor_covariates] <-
            sapply(strsplit(new$covariate, "_"), `[`, 1)[sapply(strsplit(new$covariate, "_"), `[`, 1) %in% factor_covariates]

          new <- as.character(unlist(new))

          time_var_include=c(time_var_include, new)

          message(paste0("For ", exposure, " at exposure time point ", exposure_time_pt,
                         ", the following covariate(s) will be added to the short balancing formula: "), temp, "\n")
        }else{
          message(paste0("For ", exposure, " at exposure time point ", exposure_time_pt,
                         "no time-varying confounders at additional lags were added."))
        }
      }

      vars_to_include <- c(ti_confounders, time_var_include)

      #adding in any user-specified concurrent confounders (default is only lagged)
      if(!is.na(concur_conf)){
        if(sapply(strsplit(concur_conf, "\\."), "[", 1) %in% exposure){
          stop("Do not include the exposure concurrently. Please revise the concur_conf field.")
        }
        if(sum(concur_conf %in% tv_confounders) != length(concur_conf)){
          stop("The variables in the concur_conf field must be included in the tv_confounders.")
        }
        if(as.numeric(sapply(strsplit(concur_conf, "\\."), "[", 2)) %in% time){
          vars_to_include <- c(vars_to_include,
                               concur_conf[as.numeric(sapply(strsplit(concur_conf, "\\."), "[", 2)) %in% time] )
        }
      }

      #adding in any user-specified confounders to retain in all formulas
      if(!is.na(keep_conf)){
        if(sum(keep_conf %in% tv_confounders) != length(keep_conf)){
          stop("The variables in the keep_conf field must be included in tv_confounders.")
        }
        keep_conf <- keep_conf[!paste0(exposure, ".", time) %in% keep_conf]
        vars_to_include <- c(vars_to_include, keep_conf[as.numeric(sapply(strsplit(keep_conf, "\\."), "[", 2)) < time])
      }

      # Creates form for the given exposure time point
      f <- as.formula(paste(paste0(exposure, ".", time, " ~ "),
                            paste0(vars_to_include[order(vars_to_include)], sep = "", collapse = " + ")))

      # Prints form for user inspection
      cat(paste0("The ", type, " formula for ", exposure, "-", outcome, " at ", exposure, " time point ",
                 as.character(time), " is:"), "\n")
      print(f)
      cat("\n")

      # Appends the form string to forms_csv
      forms_csv <- c(forms_csv, paste0(type, " formula for ", exposure, "-", outcome, " at ", exposure,
                                       " time point ", as.character(time), ":"))
      forms_csv <- c(forms_csv, paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep = "", collapse = " + ")))

      # Assigns the form to forms list
      assign(paste(type, "_form_", exposure, "-", outcome, "-", time, sep = ""), f, envir = parent.frame())
      forms[[paste(type, "_form_", exposure, "-", outcome, "-", time, sep = "")]] <- f
    }

    # Writes forms_csv to a CSV file
    forms_csv_file <- file.path(forms_dir, paste0(exposure, "-", outcome, "_", type, "_balancing_formulas.csv"))
    writeLines(forms_csv, con = forms_csv_file)

    # writes to rds
    forms_rds_file <- file.path(forms_dir, paste0(exposure, "-", outcome, "_", type, "_balancing_formulas.rds"))
    saveRDS(ls(pattern = type, "_form_", envir = parent.frame()), file = forms_rds_file)
  }

  forms
}
