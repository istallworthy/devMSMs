
#' Create balancing formulas
#'
#' Creates balancing formulas relating exposure to all relevant time-varying and
#' time invariant confounders at each exposure time point to be used to create
#' IPTW weights.
#'
#' @param home_dir path to home directory
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix
#' @param ti_confounders list of time invariant confounders
#' @param type type of formula to create from 'full' (includes all lagged
#'   time-varying confounders), 'short' (includes time-varying confounders at
#'   t-1 lag only), or 'update' (adds to 'short' formulas any imbalanced
#'   time-varying confounders at lags great than t-1)
#' @param bal_stats list of balance statistics from assessBalance(), required
#'   for 'update' type
#' @param concur_conf (optional) list of variable names reflecting time-varying
#'   confounders to retain in formulas contemporaneously (default is none)
#' @param keep_conf (optional) list of variable names reflecting confounders to
#'   always retain in formulas (default depends on type)
#' @param custom (optional) custom list of formulas at each exposure time point
#'   (default is to create automatically according to type)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @return list of balancing formulas at each exposure time point
#' @export
#'
#' @examples
#' test <- data.frame(A.1 = 1:10,
#' A.2 = 21:30,
#' A.3 = 1:10,
#' B.1 = 2:11,
#' B.2 = 1:10,
#' B.3 = 4:13,
#' C = 3:12,
#' D.3 = 4:13)
#' createFormulas(getwd(), "A", c(1, 2, 3), "D.3", c("B.1", "B.2", "B.3"), "C", "full")


createFormulas <- function(home_dir, exposure, exposure_time_pts, outcome, tv_confounders, ti_confounders, type, bal_stats = NULL, concur_conf = NULL, keep_conf = NULL, custom = NULL, verbose = TRUE ){

  if (missing(home_dir)){
    stop("Please supply a home directory.", call. = FALSE)
  }
  if (missing(exposure)){
    stop("Please supply a single exposure.", call. = FALSE)
  }
  if (missing(outcome)){
    stop("Please supply a single outcome.", call. = FALSE)
  }
  if (missing(exposure_time_pts)){
    stop("Please supply the exposure time points at which you wish to create weights.", call. = FALSE)
  }
  if (missing(tv_confounders)){
    stop("Please supply a list of time-varying confounders.", call. = FALSE)
  }
  if (missing(ti_confounders)){
    stop("Please supply a list of time invariant confounders.", call. = FALSE)
  }
  if (missing(type)){
    stop("Please supply a 'full', 'short', or 'update' type", call. = FALSE)
  }

  if(!is.character(type)){
    stop("Please provide a character string type from the following list: 'short', 'full', or 'update'", call. = FALSE)
  }
  if(! type %in% c("short", "full", "update")){
    stop("Please provide a type from the following list: 'short', 'full', or 'update'", call. = FALSE)
  }
  if (type != "update" & !is.null(bal_stats)){
    stop ("Please only provide balance statistics for the type 'update'.", call. = FALSE)
  }

  if(!is.null(bal_stats) & !is.data.frame(bal_stats)){
    stop("Please provide a data frame of balance statistics from the assessBalance function.", call. = FALSE)
  }
  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.", call. = FALSE)
  }


  time_varying_covariates <- tv_confounders
  all_covars <- c(tv_confounders, ti_confounders)

  if (!is.null(custom)){
    if (length(custom) != length(exposure_time_pts) | !inherits(custom, "list")){
      stop("If you wish to supply custom formulas, please provide a list of formulas for each exposure time point.", call. = FALSE)
    }

    forms <- custom
  }
  else{
    if (type != "update" & !is.null(bal_stats)){
      stop ("Please only provide balance statistics for the type 'update'.", call. = FALSE)
    }

    #create parent directory
    forms_dir1 <- file.path(home_dir, "formulas")
    if (!dir.exists(forms_dir1)) {
      dir.create(forms_dir1)
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

      else if (type == "short"){
        message("Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:")
        time_var_include <- time_varying_covariates[as.numeric(sapply(strsplit(time_varying_covariates, "\\."), "[", 2)) ==
                                                      exposure_time_pts[x-1]]
      }

      else if (type == "update"){
        if(is.null(bal_stats)){
          stop("Please provide balance statistics if you wish to run the update verison of this function", call. = FALSE)
        }

        message("Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:")

        time_var_include <- time_varying_covariates[as.numeric(sapply(strsplit(time_varying_covariates, "\\."), "[", 2)) ==
                                                      exposure_time_pts[x-1]]

        if (x > 1) {
          new <- bal_stats %>%
            dplyr::filter(balanced == 0) %>%
            dplyr::filter(exp_time == time, as.numeric(covar_time) < exposure_time_pts[x - 1],
                          as.numeric(covar_time) > 0) %>% # Finds any lagged imbalanced covars
            dplyr::select(covariate)


          # Renames factors (that were appended w/ level)
          if (nrow(new) > 0) {
            new$covariate[sapply(strsplit(new$covariate, "_"), `[`, 1) %in% factor_covariates] <-
              sapply(strsplit(new$covariate, "_"), `[`, 1)[sapply(strsplit(new$covariate, "_"), `[`, 1) %in% factor_covariates]

            new <- as.character(unlist(new))

            time_var_include <- c(time_var_include, new)

            if (verbose){
              cat(paste0("For ", exposure, " at exposure time point ", time ,
                         ", the following covariate(s) will be added to the short balancing formula: "), paste(new, collapse = ", "), "\n")
              cat("\n")
            }
          }
          else{
            if (verbose) {
              cat(paste0("For ", exposure, " at exposure time point ", time ,
                         " no time-varying confounders at additional lags were added."), "\n")
              cat("\n")
            }
          }
        }
      }

      vars_to_include <- c(ti_confounders, time_var_include)

      #adding in any user-specified concurrent confounders (default is only lagged)
      if(!is.null(concur_conf)){
        if(!inherits(concur_conf, "character")){
          stop("Please provide as a character string a list of concurrent confounders to include.", call. = FALSE)
        }
        if(sapply(strsplit(concur_conf, "\\."), "[", 1) %in% exposure){
          stop("Do not include the exposure concurrently. Please revise the concur_conf field.", call. = FALSE)
        }

        if(sum(concur_conf %in% tv_confounders) != length(concur_conf)){
          stop("The variables in the concur_conf field must be included in the tv_confounders.", call. = FALSE)
        }
        if(as.numeric(sapply(strsplit(concur_conf, "\\."), "[", 2)) %in% time){
          vars_to_include <- c(vars_to_include,
                               concur_conf[as.numeric(sapply(strsplit(concur_conf, "\\."), "[", 2)) %in% time] )
        }
      }

      #adding in any user-specified confounders to retain in all formulas
      if(!is.null(keep_conf)){
        if(!inherits(keep_conf, "character")){
          stop("Please provide as a character string a list of confounders to include in all formulas.", call. = FALSE)
        }
        if(sum(keep_conf %in% tv_confounders) != length(keep_conf)){
          stop("The variables in the keep_conf field must be included in tv_confounders.", call. = FALSE)
        }
        keep_conf <- keep_conf[!paste0(exposure, ".", time) %in% keep_conf]

        if (length(keep_conf[as.numeric(sapply(strsplit(keep_conf, "\\."), "[", 2)) < time]) > 0){
          if(! keep_conf[as.numeric(sapply(strsplit(keep_conf, "\\."), "[", 2)) < time] %in% vars_to_include){
            vars_to_include <- c(vars_to_include, keep_conf[as.numeric(sapply(strsplit(keep_conf, "\\."), "[", 2)) < time])
          }
        }
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
    forms_csv_file <- paste0(forms_dir, "/", type, "_", exposure, "-", outcome, "_", type, "_balancing_formulas.csv")
    writeLines(forms_csv, con = forms_csv_file)

    # writes to rds
    forms_rds_file <- paste0(forms_dir, "/", type, "_", exposure, "-", outcome, "_", type, "_balancing_formulas.rds")
    saveRDS(ls(pattern = type, "_form_", envir = parent.frame()), file = forms_rds_file)
  }

  forms
}
