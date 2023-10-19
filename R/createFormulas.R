
#' Create balancing formulas
#'
#' Creates balancing formulas relating exposure to all relevant time-varying and
#' time invariant confounders at each exposure time point to be used to create
#' IPTW weights.
#'
#' @param home_dir path to home directory (required if 'save.out' = TRUE)
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix, should include exposure and outcome variables (at least
#'   time-varying exposure variables required here)
#' @param ti_confounders list of time invariant confounders (at least one required)
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
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return list of balancing formulas at each exposure time point
#' @export
#'
#' @examples
#' #Full Formulas
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
#'
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
#'
#' #Short Formulas
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3"),
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     save.out = FALSE)
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     save.out = FALSE)
#'
#' c <- list("short_form-1" = as.formula(A.1 ~ C),
#'           "short_form-2" = as.formula(A.2 ~ A.1 + B.1 + C),
#'           "short_form-3" = as.formula(A.3 ~ A.2 + B.2 + C))
#'
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     custom = c,
#'                     save.out = FALSE)
#'
#' #Update Formulas
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
#' b <- assessBalance(data = test,
#'                    exposure = "A",
#'                    exposure_time_pts = c(1, 2, 3),
#'                    outcome = "D.3",
#'                    type = "weighted",
#'                    weights = w,
#'                    formulas = f,
#'                    save.out = FALSE)
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3"),
#'                     ti_confounders = "C",
#'                     type = "update",
#'                     bal_stats = b,
#'                     save.out = FALSE)
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "update",
#'                     bal_stats = b,
#'                     save.out = FALSE)
#' 

createFormulas <- function(home_dir, exposure, exposure_time_pts, outcome, type, ti_confounders, 
                           tv_confounders, bal_stats = NULL, concur_conf = NULL,
                           keep_conf = NULL, custom = NULL, verbose = TRUE, save.out = TRUE ) {
  
  if (save.out) {
    if (missing(home_dir)) {
      stop ("Please supply a home directory.",
            call. = FALSE)
    }
    else if (!is.character(home_dir)) {
      stop ("Please provide a valid home directory path as a string if you wish to save output locally.",
            call. = FALSE)
    }
    else if (!dir.exists(home_dir)) {
      stop ("Please provide a valid home directory path if you wish to save output locally.",
            call. = FALSE)
    }
  }
  
  if (missing(exposure)) {
    stop ("Please supply a single exposure.",
          call. = FALSE)
  }
  else if (!is.character(exposure) || length(exposure) != 1) {
    stop ("Please supply a single exposure as a character.",
          call. = FALSE)
  }
  else if (grepl("\\.", exposure)) {
    stop ("Please supply an exposure without the '.time' suffix or any '.' special characters. Note that the exposure variables in your dataset should be labeled with the '.time' suffix.",
          call. = FALSE)
  }
  
  if (missing(outcome)) {
    stop ("Please supply a single outcome.",
          call. = FALSE)
  }
  else if (!is.character(outcome) || length(outcome) != 1) {
    stop ("Please supply a single outcome as a character with a '.time' suffix denoting the outcome time point.",
          call. = FALSE)
  }
  else if (!grepl("\\.", outcome)) {
    stop ("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
          call. = FALSE)
  }
  else if (as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) != 
           exposure_time_pts[length(exposure_time_pts)] && 
           !as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) > 
           exposure_time_pts[length(exposure_time_pts)] ) {
    stop ("Please supply an outcome variable with a time point that is equal to or greater than the last exposure time point.",
          call. = FALSE)
  }
  
  if (missing(exposure_time_pts)) {
    stop ("Please supply the exposure time points at which you wish to create weights.",
          call. = FALSE)
  }
  else if (!is.numeric(exposure_time_pts)) {
    stop ("Please supply a list of exposure time points as integers.",
          call. = FALSE)
  }
  else if (!length(exposure_time_pts) > 1) {
    stop ("Please supply at least two exposure time points.",
          call. = FALSE)
  }
  
  if (missing(tv_confounders)) {
    stop( "Please list at least all wide exposure variables as tv_confounders.",
             call. = FALSE)
  }
  else if (!is.character(tv_confounders)) {
    stop ("Please provide a list of time-varying confounders as character strings.",
          call. = FALSE)
  }
  else if (!length(tv_confounders) > 1) {
    warning ("Please make sure to include all exposure variables as time-varying confounders.",
             call. = FALSE)
  }
  else if (any(!grepl("\\.", tv_confounders))) {
    stop ("Please list all time-varying confounders with suffix '.time' that should match variables in dataset.",
          call. = FALSE)
  }
  else if (any(!paste(exposure, exposure_time_pts, sep = ".") %in% tv_confounders)) {
    stop ("Please include all measured exposure variables in wide format in tv_confounders.",
          call. = FALSE)
  }

  else if (is.null(custom)) {
    if (any(!exposure_time_pts %in% 
            as.numeric(unlist(sapply(strsplit(tv_confounders, "\\."), "[", 2))))) {
      stop ("Exposure time points and the time points at which time-varying confounders are measured must fully overlap.",
            call. = FALSE) 
    }
  }
  
  if (missing(ti_confounders)) {
    stop ("You have not specified time invariant confounders.",
          call. = FALSE)
  }
  else if (any(grepl("\\.", ti_confounders))) {
    stop ("Time invariant confounders should not include the suffix '.time' or any '.' special characters.",
          call. = FALSE)
  }
  
  if (missing(type)) {
    stop ("Please supply a 'full', 'short', or 'update' type.",
          call. = FALSE)
  }
  else  if (!is.character(type)) {
    stop ("Please provide a character string type from the following list: 'short', 'full', or 'update'.",
          call. = FALSE)
  }
  else if (! type %in% c("short", "full", "update") || length(type) != 1) {
    stop ("Please provide a single type from the following list: 'short', 'full', or 'update'.",
          call. = FALSE)
  }
  
  if (type != "update" && !is.null(bal_stats)) {
    stop ("Please only provide balance statistics for the type 'update'.",
          call. = FALSE)
  }
  
  if (!is.null(bal_stats) && !is.data.frame(bal_stats)) {
    stop ("Please provide a data frame of balance statistics from the assessBalance function.",
          call. = FALSE)
  }
  
  if (!is.logical(verbose)) {
    stop ("Please set verbose to either TRUE or FALSE.",
          call. = FALSE)
  }
  else if (length(verbose) != 1) {
    stop ("Please provide a single TRUE or FALSE value to verbose.",
          call. = FALSE)
  }
  
  if (!is.logical(save.out)) {
    stop ("Please set save.out to either TRUE or FALSE.",
          call. = FALSE)
  }
  else if (length(save.out) != 1) {
    stop ("Please provide a single TRUE or FALSE value to save.out.",
          call. = FALSE)
  }
  
  
  all_covars <- c(tv_confounders, ti_confounders)
  
  #for custom formulas
  
  if (!is.null(custom)) {
    if (length(custom) != length(exposure_time_pts) || !is.list(custom) || 
        is.data.frame(custom)) {
      stop ("If you wish to supply custom formulas, please provide a list of formulas for each exposure time point.",
            call. = FALSE)
    }
    
    forms <- custom
    
    if (any(as.logical(unlist(lapply(forms, function(x) {
      !inherits(x, "formula") } ))))) {
      stop ("Please make sure each entry of your custom formulas list is a formula.",
            call. = FALSE)
    }
    
    if (any(names(forms) != paste(paste0(type, "_form"), exposure_time_pts, 
                                  sep = "-"))) {
      if (verbose) {
        message("Renaming custom formulas.")
      }
      
      names(forms) = paste(paste0(type, "_form"), exposure_time_pts, sep = "-")
    }
    
    
    #checking custom formulas
    
    covars <- lapply(forms, function(f) {
      deparse1(f, collapse = "") # gets covariates
    })
    
    if (any(as.logical(unlist(lapply(covars, function(x) {
      sapply(strsplit(sapply(strsplit(x, "\\~"), "[", 1), "\\."), "[", 1) != 
        exposure }))))) {
      stop ("Please make each formula in your custom formula list a function of exposure at each time point.",
            call. = FALSE)
    }
    
    if (any(as.numeric(unlist(lapply(covars, function(x) {
      sapply(strsplit(sapply(strsplit(x, "\\~"), "[", 1), "\\."), "[", 2)}))) != 
      exposure_time_pts)) {
      stop ("Please make each formula in your custom formula list a function of exposure at each time point.",
            call. = FALSE)
    }
    
    covars <- lapply(covars, function(x) {
      sapply(strsplit(x, "\\~"), "[", 2)
    })
    covars <- lapply(covars, function(x){
      as.character(unlist(strsplit(x, "\\+")))
    })
    covars <- lapply(covars, function(x) {
      gsub(" ", "", x)
    })
    
    if (any(as.logical(unlist(lapply(covars, function(x) {
      any(!x %in% tv_confounders & !x  %in% 
          ti_confounders)
    }))))) {
      miss <- lapply(covars, function(x){
        x[!x %in% tv_confounders & !x %in% 
            ti_confounders]
      })
      stop (sprintf("Please make sure all variables in your custom formulas are included as either time-varying or time invariant confounders.
           The following variables are not: %s",
                    paste(miss, sep = ", ")),
            call. = FALSE)
    }
    
  }
  
  #create formulas 
  
  else {
    if (type != "update" && !is.null(bal_stats)) {
      stop ("Please only provide balance statistics for the type 'update'.",
            call. = FALSE)
    }
    
    if (save.out) {
      forms_dir1 <- file.path(home_dir, "formulas")
      if (!dir.exists(forms_dir1)) {
        dir.create(forms_dir1)
      }
      forms_dir <- file.path(home_dir, "formulas", type)
      if (!dir.exists(forms_dir)) {
        dir.create(forms_dir)
      }
    }
    
    factor_covariates <- names(data)[sapply(data, is.factor)]
    factor_covariates <- setdiff(factor_covariates, "ID")
    
    forms_csv <- character()
    forms <- list()
    
    for (x in seq_along(exposure_time_pts)) {
      time <- exposure_time_pts[x]
      
      #identifying lagged tv confounders relative to time point based on formula type
      
      if (type == "full") {
        
        if (verbose) {
          message("USER ALERT: Please manually inspect the full balancing formula below:")
        }
        
        time_var_include <- tv_confounders[as.numeric(sapply(strsplit(tv_confounders, 
                                                                      "\\."), "[", 2)) < time]
        
      }
      
      else if (type == "short") {
        
        if (verbose) {
          message("USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:")
        }
        
        time_var_include <- tv_confounders[as.numeric(sapply(strsplit(tv_confounders, "\\."), "[", 2)) ==
                                             exposure_time_pts[x - 1]]
        
      }
      
      else if (type == "update") {
        
        if (is.null(bal_stats)) {
          stop ("Please provide balance statistics from the assessBalance() function if you wish to run the update version of this function",
                call. = FALSE)
        }
        else if (!is.data.frame(bal_stats)) {
          stop ("Please provide balance statistics from the assessBalance() function if you wish to run the update version of this function",
                call. = FALSE)
        }
        
        if (verbose) {
          message ("USER ALERT: Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and
                  those greater at further lags that remained imbalanced:")
        }
        
        time_var_include <- tv_confounders[as.numeric(sapply(strsplit(tv_confounders, "\\."), "[", 2)) ==
                                             exposure_time_pts[x - 1]]
        
        if (x > 1) {
          
          new <- bal_stats[bal_stats$balanced == 0 &
                             bal_stats$exp_time == time &
                             as.numeric(bal_stats$covar_time) < exposure_time_pts[x - 1] &
                             as.numeric(bal_stats$covar_time) > 0, ]
          
          
          if (nrow(new) > 0) {
            
            new <- new[, "covariate"]
            
            new <- as.character(unlist(new))
            
            time_var_include <- c(time_var_include, new)
            
            if (verbose) {
              
              cat(sprintf("For %s at exposure time point %s the following covariate(s) will be added to the short balancing formula:
                          %s \n",
                          exposure, time, paste(new, collapse = ", ")))
              
              cat("\n")
            }
          }
          else {
            if (verbose) {
              
              cat(sprintf("For %s at exposure time point %s no time-varying confounders at additional lags were added. \n",
                          exposure, time ))
              cat("\n")
            }
          }
        }
      }
      
      vars_to_include <- c(ti_confounders, time_var_include)
      
      #adding in any user-specified concurrent confounders (default is only lagged)
      
      if (!is.null(concur_conf)) {
        
        if (!inherits(concur_conf, "character")) {
          stop ("Please provide as a character string a list of concurrent confounders to include.",
                call. = FALSE)
        }
        
        else if (sapply(strsplit(concur_conf, "\\."), "[", 1) %in% exposure) {
          stop ("Do not include the exposure concurrently as a confounder. Please revise the concur_conf field.",
                call. = FALSE)
        }
        
        if (sum(concur_conf %in% tv_confounders) != length(concur_conf)) {
          stop ("The variables in the concur_conf field must be included in the tv_confounders.",
                all. = FALSE)
        }
        
        if (as.numeric(sapply(strsplit(concur_conf, "\\."), "[", 2)) %in% time) {
          vars_to_include <- c(vars_to_include,
                               concur_conf[as.numeric(sapply(strsplit(concur_conf, 
                                                                      "\\."), "[", 2)) %in% time] )
        }
      }
      
      #adding in any user-specified confounders to retain in all formulas
      
      if (!is.null(keep_conf)) {
        
        if (!is.character(keep_conf)) {
          stop ("Please provide as a character string a list of confounders to include in all formulas.",
                call. = FALSE)
          
        }
        
        if (sum(keep_conf %in% tv_confounders) != length(keep_conf)) {
          stop ("The variables in the keep_conf field must be included in tv_confounders.",
                call. = FALSE)
        }
        
        keep_conf <- keep_conf[!paste0(exposure, ".", time) %in% keep_conf]
        
        if (length(keep_conf[as.numeric(sapply(strsplit(keep_conf, "\\."), 
                                               "[", 2)) < time]) > 0) {
          
          if (! keep_conf[as.numeric(sapply(strsplit(keep_conf, "\\."), "[", 2)) 
                          < time] %in% vars_to_include) {
            vars_to_include <- c(vars_to_include, 
                                 keep_conf[as.numeric(sapply(strsplit(keep_conf, 
                                                                      "\\."), "[", 2)) < time])
          }
        }
      }
      
      if (any(duplicated(vars_to_include))) {
        stop (sprintf("The following variable(s) are duplicated: %s.",
                      paste(vars_to_include[duplicated(vars_to_include)], 
                            collapse = ", ")))
      }
      
      # Creates form for the given exposure time point
      
      f <- as.formula(paste(paste0(exposure, ".", time, " ~ "),
                            paste0(vars_to_include[order(vars_to_include)], 
                                   sep = "", collapse = " + ")))
      
      # Prints form for user inspection
      if (verbose) {
        cat(sprintf("The %s formula for %s - %s at %s time point %s is: \n",
                    type, exposure, outcome, exposure, as.character(time)))
        print(f)
        cat("\n")
      }
      
      # Appends the form string to forms_csv
      
      if (save.out) {
        forms_csv <- c(forms_csv,
                       sprintf("%s formula for %s-%s at %s time point %s:",
                               type, exposure, outcome, exposure, 
                               as.character(time)))
        
        forms_csv <- c(forms_csv, paste(exposure, "~",
                                        paste0(vars_to_include[order(vars_to_include)], 
                                               sep = "", collapse = " + ")))
        
      }
      
      # Assigns the form to forms list
      
      forms[[paste(type, "_form", "-", time, sep = "")]] <- f
      
    }
    
    if (save.out) {
      
      # Writes forms_csv to a CSV file
      
      forms_csv_file <- file.path(forms_dir, 
                                  sprintf("%s_%s-%s_%s_balancing_formulas.csv",
                                          type, exposure, outcome, type))
      
      writeLines(forms_csv, con = forms_csv_file)
      
      # writes to rds
      
      forms_rds_file <- file.path(forms_dir, 
                                  sprintf("%s_%s-%s_%s_balancing_formulas.rds",
                                          type, exposure, outcome, type))
      
      saveRDS(forms, file = forms_rds_file)
    }
  }
  
  forms
}
