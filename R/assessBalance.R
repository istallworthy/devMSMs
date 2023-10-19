#' Assesses confounder balancing
#'
#' Draws on functions from the cobalt package to quantify the relations between
#' exposure and confounders at each exposure time point according to the
#' guidelines from Jackson, 2016 on how to assess balance for time-varying
#' exposures.
#'
#' @seealso {[cobalt] package,
#'   <https://cran.r-project.org/web/packages/cobalt/index.html>}
#' @seealso {Jackson, 2016 for more on assessing balance for time-varying
#'   exposures, <https://pubmed.ncbi.nlm.nih.gov/27479649/>}
#' @param home_dir (optional) path to home directory (required if save.out =
#'   TRUE)
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param formulas list of balancing formulas at each time point output from
#'   createFormulas()
#' @param weights list of IPTW weights output from createWeights, required for
#'   type 'weighted'
#' @param type type of balance assessment; 'prebalance' or 'weighted'
#' @param balance_thresh (optional) one or two numbers between 0 and 1
#'   indicating a single balancing threshold or thresholds for more and less
#'   important confounders, respectively (default = 0.1)
#' @param imp_conf (optional) list of variable names reflecting important
#'   confounders, required if two balance thresholds are supplied
#' @param verbose (optiona) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @returns a data frame of balance statistics
#' @export
#' @examples
#'
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
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "short",
#'                     save.out = FALSE)
#'
#' #Prebalance
#' b <- assessBalance(data = test,
#'                    exposure = "A",
#'                    exposure_time_pts = c(1, 2, 3),
#'                    outcome = "D.3",
#'                    type = "prebalance",
#'                    formulas = f,
#'                    save.out = FALSE)
#' b <- assessBalance(data = test,
#'                    exposure = "A",
#'                    exposure_time_pts = c(1, 2, 3),
#'                    outcome = "D.3",
#'                    type = "prebalance",
#'                    formulas = f,
#'                    balance_thresh = 0.2,
#'                    save.out = FALSE)
#' b <- assessBalance(data = test,
#'                    exposure = "A",
#'                    exposure_time_pts = c(1, 2, 3),
#'                    outcome = "D.3",
#'                    type = "prebalance",
#'                    formulas = f,
#'                    balance_thresh = c(0.1, 0.2),
#'                    imp_conf = "B.1",
#'                    save.out = FALSE)
#'
#' # Weighted
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
#' b <- assessBalance(data = test,
#'                    exposure = "A",
#'                    exposure_time_pts = c(1, 2, 3),
#'                    outcome = "D.3",
#'                    type = "weighted",
#'                    weights = w,
#'                    formulas = f,
#'                    balance_thresh = 0.2,
#'                    save.out = FALSE)
#' b <- assessBalance(data = test,
#'                    exposure = "A",
#'                    exposure_time_pts = c(1, 2, 3),
#'                    outcome = "D.3",
#'                    type = "weighted",
#'                    weights = w,
#'                    formulas = f,
#'                    balance_thresh = c(0.1, 0.2),
#'                    imp_conf = "B.1",
#'                    save.out = FALSE)


assessBalance <- function(home_dir, data, exposure, exposure_time_pts, outcome, type, formulas, weights = NULL, balance_thresh = NULL,
                          imp_conf = NULL, verbose = TRUE, save.out = TRUE) {
  
  if (!is.logical(save.out)) {
    stop ("`save.out` must be a flag (TRUE or FALSE)",
          call. = FALSE)
  }
  if (save.out) {
    if (missing(home_dir)) {
      stop ("Please supply a home directory.",
            call. = FALSE)
    }
    
    else if(!is.character(home_dir)){
      stop ("Please provide a valid home directory path as a string if you wish to save output locally.",
            call. = FALSE)
    }
    else if(!dir.exists(home_dir)) {
      stop ("Please provide a valid home directory path if you wish to save output locally.",
            call. = FALSE)
    }
  }
  
  if (missing(data)) {
    stop ("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
          call. = FALSE)
  }
  
  else if (!inherits(data, "mids") && !is.data.frame(data) &&
           !(is.list(data) && all(vapply(data, is.data.frame, logical(1L))))) {
    stop ("Please provide either a 'mids' object, a data frame, or a list of imputed csv files in the 'data' field.",
          call. = FALSE)
  }
  
  else if (is.list(data) && !is.data.frame(data)  && !mice::is.mids(data)) {
    if (sum(sapply(data, is.data.frame)) != length(data)) {
      stop ("Please supply a list of data frames that have been imputed.",
            call. = FALSE)
    }
  }
  
  if (missing(exposure)) {
    stop ("Please supply a single exposure.",
          call. = FALSE)
  }
  else if (!is.character(exposure) | length(exposure) != 1) {
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
  else if (!is.character(outcome) | length(outcome) != 1) {
    stop ("Please supply a single outcome as a character.",
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
  else if(!is.numeric(exposure_time_pts)) {
    stop ("Please supply a list of exposure time points as integers.",
          call. = FALSE)
  }
  else if (!length(exposure_time_pts) > 1) {
    stop ("Please supply at least two exposure time points.",
          call. = FALSE)
  }
  
  if (missing(formulas)) {
    stop ("Please supply a list of balancing formulas.",
          call. = FALSE)
  }
  else if(!is.list(formulas) | is.data.frame(formulas)) {
    stop ("Please provide a list of formulas for each exposure time point",
          call. = FALSE)
  }
  else if(length(formulas) != length(exposure_time_pts)) {
    stop ("Please provide a list of formulas for each exposure time point",
          call. = FALSE)
  }
  else if (is.list(formulas) && !is.data.frame(formulas)) {
    if (sum(sapply(formulas, function(x) {
      
      inherits(x, "formula")})) != length(formulas)) {
      stop ("Please supply a list of formulas for each exposure time point.",
            call. = FALSE)
    }
  }
  
  if (missing(type)) {
    stop ("Please supply a 'weighted', 'prebalance' type",
          call. = FALSE)
  }
  if (!inherits(type, "character") || length(type) != 1 ) {
    stop ("Please provide a single type as a character string from the following list: 'prebalance', 'weighted'",
          call. = FALSE)
  }
  else if (!type %in% c("prebalance", "weighted")) {
    stop ("Please provide a type from the following list: 'prebalance', 'weighted'",
          call. = FALSE)
  }
  else if (type == "prebalance" && !is.null(weights)) {
    stop ("The 'prebalance' mode of this function assesses balance prior to weighting and thus does not take weights.",
          call. = FALSE)
  }
  else if (type == "weighted" && (is.null(weights) || missing(weights))) {
    stop ("The 'weighted' mode of this function requires weights be supplied in the form of output from createWeights.",
          call. = FALSE)
  }
  
  if (!inherits(data, "mids") && !is.data.frame(data) &&
      !(is.list(data) && all(vapply(data, is.data.frame, logical(1L))))) {
    stop ("Please provide either a 'mids' object, a data frame, or a list of imputed csv files in the 'data' field.",
          call. = FALSE)
  }
  
  
  if (!is.null(weights) && (!is.list(weights) || is.data.frame(weights))) {
    stop ("Please supply a list of weights output from the createWeights function.",
          call. = FALSE)
  }
  else if (is.list(weights) && !is.data.frame(weights)) {
    if (sum(sapply(weights, function(x) {
      
      inherits(x, "weightitMSM")})) != length(weights)) {
      stop ("Please supply a list of weights output from the createWeights function.",
            call. = FALSE)
    }
  }
  
  if (!is.numeric(balance_thresh) && !is.null(balance_thresh)) {
    stop ("Please provide one or two balance thresholds as numbers from 0-1.",
          call. = FALSE)
  }
  else if (!is.null(balance_thresh) && !length(balance_thresh < 3)) {
    stop ("Please provide one or two balance thresholds as numbers from 0-1.",
          call. = FALSE)
  }
  else if (length(balance_thresh) == 2 && is.null(imp_conf)) {
    stop ("If you wish to provide different balance threshold for important and less important confounders, please provide a list of important confounders in the 'imp_conf' field.",
          call. = FALSE)
  }
  if (is.null(balance_thresh)) {
    balance_thresh <- 0.1
  }
  
  if (!is.null(imp_conf) && length(balance_thresh) == 1 && !is.null(balance_thresh) ) {
    stop ("If you provide a list of important confounders, please provide a list of two balance thresholds for important and less important confounders, respectively",
          call. = FALSE)
  }
  else if (!is.null(imp_conf) && !is.character(imp_conf)) {
    stop ("Please provide a list variable names as characters that are important confounders.",
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
    stop("Please set save.out to either TRUE or FALSE.",
         all. = FALSE)
  }
  else if (length(save.out) != 1) {
    stop ("Please provide a single TRUE or FALSE value to save.out.",
          call. = FALSE)
  }
  
  mi <- !is.data.frame(data)
  
  folder <- switch (type, "prebalance" = "prebalance/", "weighted/")
  
  
  if (save.out) {
    balance_dir <- file.path(home_dir, "balance")
    if (!dir.exists(balance_dir)) {
      dir.create(balance_dir)
    }
    folder_dir <- file.path(home_dir, "balance", folder)
    if (!dir.exists(folder_dir)) {
      dir.create(folder_dir)
    }
    plots_dir <- file.path(home_dir, "balance", folder, "plots")
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir)
    }
  }
  
  weights_method <- switch(type, "prebalance" = "no weights", weights[[1]]$method)
  
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[", 1)
  
  
  ### Getting balance statistics
  
  if (type == "prebalance") {
    
    if (verbose) {
      message(sprintf("USER ALERT: The following statistics display covariate imbalance at each exposure time point prior to weighting,
            using %s formulas.\n",
                      form_name))
      cat("\n")
    }
    
    # Running balance stats function, unweighted, on each imputed dataset
    
    if (mi) {
      
      if (inherits(data, "mids")) {
        m <- data$m
        
        bal_stats <- lapply(seq_len(m), function(k) {
          
          d <- as.data.frame(mice::complete(data, k))
          
          if (anyNA(d)) {
            stop ("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
                  call. = FALSE)
          }
          
          exposure_type <- if (is.numeric(d[, paste0(exposure, '.', exposure_time_pts[1])])) "continuous" else "binary"
          
          if (sum(duplicated(d$"ID")) > 0) {
            stop ("Please provide wide imputed datasets with a single row per ID.",
                  call. = FALSE)
          }
          
          calcBalStats(home_dir, d, formulas, exposure, exposure_time_pts, outcome,
                       balance_thresh, k = k, weights = NULL, imp_conf, verbose, save.out)
        })
      }
      else {
        
        m <- length(data)
        
        bal_stats <- lapply(seq_len(m), function(k) {
          d <- data[[k]]
          
          if (anyNA(d)) {
            stop ("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
                  call. = FALSE)
          }
          
          exposure_type <- if (is.numeric(d[[paste0(exposure, '.', exposure_time_pts[1])]])) "continuous" else "binary"
          
          if (sum(duplicated(d[["ID"]])) > 0) {
            stop ("Please provide wide imputed datasets with a single row per ID.", 
                  call. = FALSE)
          }
          
          calcBalStats(home_dir, d, formulas, exposure, exposure_time_pts, outcome,
                       balance_thresh, k = k, weights = NULL, imp_conf, verbose, save.out)
        })
      }
      
      # Save out balance stats for each imputed dataset
      
      bal_stats_all_imp <- do.call(dplyr::bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      
      
      if (save.out) {
        write.csv(bal_stats_all_imp, 
                  file.path(home_dir, "balance", type, sprintf("%s_all_imps_balance_stat_summary.csv",
                                                               exposure)))
        if (verbose ) {
          cat("Pre balance statistics for each imputed dataset have now been saved in the 'balance/prebalance/' folder", "\n")
        }
        
      }
      
      
      # Gathering imbalanced covariate statistics to average across imputed datasets for the final list/assessment of imbalanced covariates
      # Averaging across imputed datasets
      
      all_bal_stats <- data.frame(
        exposure = exposure,
        exp_time = bal_stats[[1]]$exp_time,
        covar_time = bal_stats[[1]]$covar_time,
        covariate = bal_stats[[1]]$covariate,
        avg_bal = rowMeans(do.call(cbind, lapply(bal_stats, `[[`, "std_bal_stats"))))
      
      #adds custom bal thresh info
      
      if (!is.null(imp_conf)) {
        
        all_bal_stats$bal_thresh <- ifelse(all_bal_stats$covariate %in% imp_conf,
                                           balance_thresh[1], balance_thresh[2])
        
        all_bal_stats$balanced <- ifelse(abs(all_bal_stats$avg_bal) <
                                           all_bal_stats$bal_thresh, 1, 0)
      }
      else {
        all_bal_stats$bal_thresh <- balance_thresh
        
        all_bal_stats$balanced <- ifelse(abs(all_bal_stats$avg_bal) <
                                           all_bal_stats$bal_thresh, 1, 0)
      }
      
      
      if (verbose) {
        cat("\n")
        cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      }
      
      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }
    else {
      
      if (sum(duplicated(data[["ID"]])) > 0) {
        stop ("Please provide wide dataset with a single row per ID.",
              call. = FALSE)
      }
      
      if (anyNA(data)) {
        stop ("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
              call. = FALSE)
      }
      
      exposure_type <- if (is.numeric(data[[paste0(exposure, '.', exposure_time_pts[1])]])) "continuous" else "binary"
      
      bal_stats <- calcBalStats(home_dir, data, formulas, exposure, 
                                exposure_time_pts, outcome, balance_thresh, 
                                k = 0, weights = NULL, imp_conf, verbose, save.out)
      
      # Gathering imbalanced covariate statistics for the final list/assessment of imbalanced covariates
      
      all_bal_stats <- data.frame(
        exposure = exposure,
        exp_time = bal_stats$exp_time,
        covar_time = bal_stats$covar_time,
        covariate = bal_stats$covariate,
        avg_bal = bal_stats$std_bal_stats,
        bal_thresh = bal_stats$bal_thresh,
        balanced = bal_stats$balanced)
      
      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }
  } #ends prebalance
  
  
  # Weighted
  
  else if (type == "weighted") {
    
    if (verbose) {
      message(sprintf("USER ALERT: The following statistics display covariate imbalance at each exposure time point following IPTW weighting,
            using %s formulas.\n", form_name))
      cat("\n")
    }
    
    if (mi) {
      
      # Running balance stats function, unweighted, on each imputed dataset w/ no weights
      
      if (inherits(data, "mids")) {
        m <- data$m
        
        bal_stats <- lapply(seq_len(m), function(k) {
          d <- as.data.frame(mice::complete(data, k))
          
          if (anyNA(d)) {
            stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
                 call. = FALSE)
          }
          if (sum(duplicated(d[["ID"]])) > 0) {
            stop ("Please provide wide imputed datasets with a single row per ID.", 
                  call. = FALSE)
          }
          
          w <- weights[[k]]
          
          calcBalStats(home_dir, d, formulas, exposure, exposure_time_pts, outcome,
                       balance_thresh, k = k, weights = w, imp_conf, verbose, save.out)
        })
      }
      else {
        m <- length(data)
        
        bal_stats <- lapply(seq_len(m), function(k) {
          d <- data[[k]]
          
          if (anyNA(d)) {
            stop ("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).",
                  call. = FALSE)
          }
          if (sum(duplicated(d[["ID"]])) > 0) {
            stop ("Please provide wide imputed datasets with a single row per ID.", 
                  call. = FALSE)
          }
          
          w <- weights[[k]]
          
          calcBalStats(home_dir, d, formulas, exposure, exposure_time_pts, outcome,
                       balance_thresh, k = k, weights = w, imp_conf, verbose, save.out)
        })
      }
      
      
      # Save out balance stats for each imputed dataset
      
      bal_stats_all_imp <- do.call(dplyr::bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      
      if (save.out) {
        write.csv(bal_stats_all_imp, 
                  file.path(home_dir, "balance", type, sprintf("%s_all_imps_balance_stat_summary.csv",
                                                               exposure)))
        if (verbose) {
          cat("Weighted balance statistics for each imputed dataset have now been saved in the 'balance/weighted/' folder", "\n")
        }
      }
      
      # Gathering imbalanced covariate statistics to average across imputed datasets for the final list/assessment of imbalanced covariates
      # Averaging across imputed datasets
      
      all_bal_stats <- data.frame(
        exposure = exposure,
        exp_time = bal_stats[[1]]$exp_time,
        covar_time = bal_stats[[1]]$covar_time,
        covariate = bal_stats[[1]]$covariate,
        avg_bal = rowMeans(do.call(cbind, lapply(bal_stats, `[[`, "std_bal_stats"))))
      
      # #adds custom bal thresh info
      
      if (!is.null(imp_conf)) {
        
        all_bal_stats$bal_thresh <-ifelse(all_bal_stats$covariate %in% imp_conf,
                                          balance_thresh[1], balance_thresh[2])
        
        all_bal_stats$balanced <- ifelse(abs(all_bal_stats$avg_bal) <
                                           all_bal_stats$bal_thresh, 1, 0)
      }
      else {
        all_bal_stats$bal_thresh <- balance_thresh
        
        all_bal_stats$balanced <- ifelse(abs(all_bal_stats$avg_bal) <
                                           all_bal_stats$bal_thresh, 1, 0)
      }
      
      if (verbose) {
        cat("\n")
        cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      }
      
      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }
    else {
      
      if (sum(duplicated(data$"ID")) > 0) {
        stop ("Please provide wide dataset with a single row per ID.",
              call. = FALSE)
      }
      
      w <- weights[[1]]
      bal_stats <- calcBalStats(home_dir, data, formulas, exposure, 
                                exposure_time_pts, outcome, balance_thresh, k = 0, 
                                weights = w, imp_conf, verbose, save.out)
      
      # Gathering imbalanced covariate statistics for the final list/assessment of imbalanced covariates
      
      all_bal_stats  <- data.frame(
        exposure = exposure,
        exp_time = bal_stats$exp_time,
        covar_time = bal_stats$covar_time,
        covariate = bal_stats$covariate,
        avg_bal = bal_stats$std_bal_stats,
        bal_thresh = bal_stats$bal_thresh,
        balanced = bal_stats$balanced)
      
      # Getting totals
      
      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }
  } #ends weighted
  
  ### Plotting and summarizing
  
  tot_cons <- unique(tot_covars) # Total domains/constructs
  
  # Make love plot to summarize imbalance at each exposure time point
  
  data_type <- if (mi) "imputed" else "single"
  
  weights_method <- if (type == "weighted") weights[[1]]$method else "no weights"
  
  sapply(seq_along(exposure_time_pts), function(i) {
    exposure_time_pt <- exposure_time_pts[i]
    
    if (inherits(data, "mids")) {
      exposure_type <- if (is.numeric(mice::complete(data, 1)[, paste0(exposure, '.', exposure_time_pts[1])])) "continuous" else "binary"
      
    }
    else if (!is.data.frame(data)) {
      exposure_type <- if (is.numeric(data[[1]][, paste0(exposure, '.', exposure_time_pts[1])])) "continuous" else "binary"
    }
    else {
      exposure_type <- if (is.numeric(data[[paste0(exposure, '.', exposure_time_pts[1])]])) "continuous" else "binary"
    }
    
    
    temp <- all_bal_stats[all_bal_stats$exp_time == exposure_time_pt, , drop = FALSE]
    
    make_love_plot(home_dir, folder, exposure, exposure_time_pt, exposure_type, 
                   k = 0, form_name, temp, data_type, balance_thresh, 
                   weights_method, imp_conf, verbose, save.out)
  })
  
  if (save.out) {
    outfile <- file.path(home_dir, "balance", type, 
                         sprintf("%s_%s-%s_all_%s_%s_associations.html",
                                 form_name, exposure, outcome, type, 
                                 weights_method))
    
    # Save out all correlations/std mean differences
    
    sink(outfile)
    stargazer::stargazer(all_bal_stats,
                         type = "html",
                         digits = 2,
                         column.labels = colnames(all_bal_stats),
                         summary = FALSE,
                         rownames = FALSE,
                         header = FALSE,
                         out = outfile)
    sink()
    
    
    if (verbose) {
      if (mi) {
        cat(sprintf("Summary plots for %s %s averaged across all imputations have been saved out for each time point in the 'balance/%s/plots/' folder.\n",
                    form_name, exposure, type))
      }
      else {
        cat(sprintf("Summary plots for %s %s have now been saved out for each time point in the 'balance/%s/plots/' folder.\n",
                    form_name, exposure, type))
      }
      cat("\n")
    }
  }
  
  if (save.out) {
    
    # Saving out all pre-balance associations
    
    write.csv(all_bal_stats, 
              file.path(home_dir, "balance", type, 
                        sprintf("%s_%s_%s_stat_summary.csv",
                                exposure, type, weights_method)),
              row.names = FALSE)
    
    
    if (verbose) {
      if (mi) {
        cat(sprintf("Check 'balance/%s/' folder for a table of all %s correlations or standardized mean differences averaged across imputed datasets.\n",
                    type, type))
      }
      else {
        cat(sprintf("Check 'balance/%s/' folder for a table of all %s correlations or standardized mean differences.\n",
                    type, type))
      }
      cat("\n")
    }
  }
  
  
  # Finding all imbalanced variables
  
  unbalanced_covars <- all_bal_stats[all_bal_stats$balanced == 0, , drop = FALSE]
  
  unbalanced_constructs  <- sapply(strsplit(unbalanced_covars$covariate, "\\."),
                                   "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, 
                                                                       "\\."), "[", 1))]
  
  if (save.out) {
    summ <- rbind(exposure = exposure,
                  form_type = form_name,
                  weights_method = weights_method,
                  med_bal_stat = round(median(abs(all_bal_stats$avg_bal)), 4),
                  min_bal_stat = round(median(abs(all_bal_stats$avg_bal)), 4),
                  max_bal_stat = round(max(all_bal_stats$avg_ba), 4))
  }
  
  if (mi) {
    if (verbose) {
      
      cat(sprintf("USER ALERT: Averaging across all imputed datasets for exposure %s using the %s formulas and %s :\n",
                  exposure, form_name, weights_method))
      
      cat(sprintf("The median absolute value relation between exposure and confounder is %s (range = %s-%s).\n",
                  round(median(abs(all_bal_stats$avg_bal)), 2),
                  round(min(all_bal_stats$avg_ba), 2),
                  round(max(all_bal_stats$avg_ba), 2)))
    }
    
    if (nrow(unbalanced_covars) > 0) {
      
      if (save.out) {
        summ <- rbind(summ,
                      n_imbal_confounders = nrow(unbalanced_covars),
                      n_confounders = length(tot_covars),
                      perc_imbal_confounders = round(nrow(unbalanced_covars) / 
                                                       length(tot_covars) * 100, 4),
                      n_imbal_domains = length(unbalanced_constructs),
                      n_domains =  length(tot_cons),
                      perc_imbal_domains = round(length(unbalanced_constructs) / 
                                                   length(tot_cons) * 100, 4),
                      imbal_med_bal_stat = round(median(abs(as.numeric(unbalanced_covars$avg_bal[unbalanced_covars$balanced == 0]))), 4),
                      imbal_min_bal_stat = round(min(unbalanced_covars$avg_bal), 4),
                      imbal_max_bal_stat = round(max(unbalanced_covars$avg_bal), 4)
        )
      }
      
      if (verbose) {
        cat("\n")
        cat(sprintf("As shown below, the following %s covariates across time points out of %s total (%s%%) spanning %s domains out of %s (%s%%) are imbalanced with a remaining median absolute value correlation/std mean difference in relation to %s of %s (range=%s-%s) :",
                    nrow(unbalanced_covars),
                    length(tot_covars),
                    round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2),
                    length(unbalanced_constructs),
                    length(tot_cons),
                    round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
                    exposure,
                    round(median(abs(as.numeric(unbalanced_covars$avg_bal[unbalanced_covars$balanced == 0]))), 2),
                    round(min(unbalanced_covars$avg_bal), 2),
                    round(max(unbalanced_covars$avg_bal), 2)))
      }
      
    }
    else {
      if (save.out) {
        summ <- rbind(summ,
                      n_imbal_confounders = 0)
      }
      if (verbose) {
        cat("There are no imbalanced covariates.", "\n")
      }
      
    }
  }
  else { #non imputed
    if (verbose) {
      cat(sprintf("USER ALERT: For exposure %s using the %s formulas and %s :\n",
                  exposure, form_name, weights_method))
      
      cat(sprintf("The median absolute value relation between exposure and confounder is %s (range = %s -%s). \n",
                  round(median(abs(all_bal_stats$avg_bal)), 2),
                  round(min(all_bal_stats$avg_ba), 2),
                  round(max(all_bal_stats$avg_ba), 2)))
    }
    
    if (nrow(unbalanced_covars) > 0) {
      
      if (save.out) {
        summ <- rbind(summ,
                      n_imbal_confounders = nrow(unbalanced_covars),
                      n_confounders = length(tot_covars),
                      perc_imbal_confounders = round(nrow(unbalanced_covars) / 
                                                       length(tot_covars) * 100, 4),
                      n_imbal_domains = length(unbalanced_constructs),
                      n_domains =  length(tot_cons),
                      perc_imbal_domains = round(length(unbalanced_constructs) / 
                                                   length(tot_cons) * 100, 4),
                      imbal_med_bal_stat =  round(median(abs(as.numeric(unlist(unbalanced_covars[unbalanced_covars$balanced == 0, ][, "avg_bal"])))), 4),
                      imbal_min_bal_stat = round(min(unbalanced_covars$avg_bal), 4),
                      imbal_max_bal_stat = round(max(unbalanced_covars$avg_bal), 4)
        )
      }
      
      if (verbose) {
        cat(sprintf("As shown below, the following %s covariates across time points out of %s total (%s%%) spanning %s domains out of %s (%s%%) are imbalanced with a remaining median absolute value correlation/std mean difference in relation to %s of %s (range=%s-%s) :\n",
                    nrow(unbalanced_covars),
                    length(tot_covars),
                    round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2),
                    length(unbalanced_constructs),
                    length(tot_cons),
                    round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
                    exposure,
                    round(median(abs(as.numeric(unlist(unbalanced_covars[unbalanced_covars$balanced == 0, ][, "avg_bal"])))), 2),
                    round(min(unbalanced_covars$avg_bal), 2),
                    round(max(unbalanced_covars$avg_bal), 2)))
      }
    }
    else {
      if (save.out) {
        summ <- rbind(summ,
                      n_imbal_confounders = 0)
      }
      if (verbose) {
        cat("There are no imbalanced covariates.", "\n")
      }
    }
  }
  
  
  if (nrow(unbalanced_covars) > 0) {
    if (verbose){
      cat("\n")
      cat("\n")
      cat(knitr::kable(unbalanced_covars,
                       caption = "Imbalanced Covariates",
                       format = 'pipe'), sep = "\n")
      
      cat("\n")
    }
    
    if (save.out) {
      
      # save out summary balance assessment
      
      k <- knitr::kable(summ,
                        caption = sprintf("Summary of %s Exposure for %s using %s",
                                          exposure, form_name, weights_method),
                        format = 'html')
      kableExtra::save_kable(k,
                             file = file.path(home_dir, "balance", type, 
                                              sprintf("%s_%s-%s_%s_%s_balance_summary.html",
                                                      form_name, exposure, outcome, 
                                                      type, weights_method)))
      
      
      # Save out only imbalanced covariates
      
      outfile <- file.path(home_dir, "balance", type, 
                           sprintf("%s_%s-%s_%s_%s_all_covariates_imbalanced.html",
                                   form_name, exposure, outcome, type, 
                                   weights_method))
      
      sink(outfile)
      stargazer::stargazer(unbalanced_covars,
                           type = "html",
                           digits = 2,
                           column.labels = colnames(unbalanced_covars),
                           summary = FALSE,
                           rownames = FALSE,
                           header = FALSE,
                           out = outfile)
      sink()
    }
  }
  else {
    if (verbose) {
      cat("There are no imbalanced covariates.")
    }
  }
  
  all_bal_stats
}



