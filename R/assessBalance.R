#' Assesses confounder balancing
#'
#' @importFrom knitr kable
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 guide_axis
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr bind_rows
#' @importFrom dplyr %>%
#' @seealso {[cobalt] package, <url1>}
#' @seealso {Jackson, 2016 for more on assessing balance for time-varying exposures, <url1>}
#' @param home_dir path to home directory
#' @param data data in wide format as: a data frame, path to folder of imputed .csv files, or mids object
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint" suffix
#' @param formulas list of balancing formulas at each time point output from createFormulas()
#' @param weights list of IPTW weights output from createWeights, required for type 'weighted'
#' @param type type of balance assessment 'prebalance' or 'weighted'
#' @param balance_thresh (optional) one or two numbers between 0 and 1 indicating a single balancingn threshold or thresholds for more and less important confounders, respectively (default = 0.1)
#' @param imp_conf (optional) list of variable names reflecting important confounders, required if two balance thresholds are supplied
#' @param verbose TRUE or FALSE indicator for user output
#' @returns a list data frame of balance statistics

assessBalance <- function(home_dir, data, exposure, exposure_time_pts, outcome, tv_confounders, type, formulas, weights = NULL, balance_thresh = 0.1, imp_conf = NULL, verbose = TRUE){

  if (missing(home_dir)){
    stop("Please supply a home directory.", call. = FALSE)
  }
  if (missing(data)){
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
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
  if (missing(type)){
    stop("Please supply a 'weighted', 'prebalance' type", call. = FALSE)
  }
  if (missing(formulas)){
    stop("Please supply a list of balancing formulas.", call. = FALSE)
  }

  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.", call. = FALSE)
  }

  if (!inherits(type, "character") | length(type) != 1 ){
    stop("Please provide a single type as a character string from the following list: 'prebalance', 'weighted'", call. = FALSE)
  }
  else if(!type %in% c("prebalance", "weighted")){
    stop("Please provide a type from the following list: 'prebalance', 'weighted'", call. = FALSE)
  }
  else  if (type == "prebalance" & !is.null(weights)){
    stop("The 'prebalance' mode of this function assesses balance prior to weighting and thus does not take weights.", call. = FALSE)
  }

  if (!mice::is.mids(data) & !is.data.frame(data) & !is.character(data)) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.", call. = FALSE)
  }

  if (!is.null(weights) & ! inherits(weights, "list")){
    stop("Please supply a list of weights output from the createWeights function.", call. = FALSE)
  }

  if(!is.numeric(balance_thresh)){
    stop("Please provide one or two balance thresholds as numbers from 0-1.")
  }
  if (length(balance_thresh) == 2 & is.null(imp_conf)){
    stop("If you wish to provide different balance threshold for important and less important confounders, please provide a list of important confounders in the 'imp_conf' field.", call. = FALSE)
  }

  if (!is.null(imp_conf) & length(balance_thresh) == 1){
    stop("If you provide a list of important confounders, please provide a list of two balance tresholds for important and less important confounders, respectively", call. = FALSE)
  }
  else if(!is.null(imp_conf) & !is.character(imp_conf)){
    stop("Please provide a list variable names as characters that are important confounders.", call. = FALSE)
  }

  if(!inherits(formulas, "list")){
    stop("Please provide a list of formulas for each exposure time point", call. = FALSE)
  }


  folder <- ifelse(type == "prebalance", "prebalance/", "weighted/")

  #creating required directories
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

  weights_method <- ifelse(type == "prebalance", "no weights", weights[[1]]$method)

  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[",1)


  if (is.character(data)) {
    if (!dir.exists(data)) {
      stop("Please provide a valid directory path with imputed .csv's, a data frame, or a 'mids' object for the 'data' field.", call. = FALSE)
    }
    else if (length(dir(data)) < 2) {
      stop("If you specify data as a directory, please supply more than 1 imputed dataset in addition to the original data.", call. = FALSE)
    }

    # List imputed files
    files <- list.files(data, full.names = TRUE, pattern = "\\.csv")

    # Read and process imputed datasets as list
    data <- lapply(files, function(file) {
      imp_data <- read.csv(file)

      if (sum(duplicated(imp_data$"ID")) > 0){
        stop("Please provide a directory to data that are imputed in wide format with a single row per ID.", call. = FALSE)
      }

      imp_data
    })
  }


  ### Getting balance statistics
  if (type == "prebalance"){
    if (verbose){
      message(paste0("USER ALERT: The following statistics display covariate imbalance at each exposure time point prior to weighting,
            using ", form_name, " formulas."), "\n")
      cat("\n")
    }

    # Running balance stats function, unweighted, on each imputed dataset
    if (mice::is.mids(data) | inherits(data, "list")){

      if (mice::is.mids(data)){
        m = data$m

        bal_stats <- lapply(seq_len(m), function(k) {
          d <- as.data.frame(mice::complete(data,k))
          exposure_type <- ifelse(inherits(d[, paste0(exposure, '.',
                                                   exposure_time_pts[1])], "numeric"), "continuous", "binary")

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputed datasets with a single row per ID.", call. = FALSE)
          }

          calcBalStats(d, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = k, weights = NULL, imp_conf, verbose)
        })
      }

      else if (inherits(data, "list")){
        m = length(data)

        bal_stats <- lapply(seq_len(m), function(k) {
          d <- data[[k]]

          exposure_type <- ifelse(inherits(d[, paste0(exposure, '.',
                                                   exposure_time_pts[1])], "numeric"), "continuous", "binary")

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputd datasets with a single row per ID.", call. = FALSE)
          }

          calcBalStats(d, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = k, weights = NULL, imp_conf, verbose)
        })
      }


      # Save out balance stats for each imputed dataset
      bal_stats_all_imp <- do.call(bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      write.csv(bal_stats_all_imp, paste0(home_dir, "/balance/", type, "/",
                                          exposure, "_all_imps_balance_stat_summary.csv"), row.names = FALSE)

      if (verbose){
        message("Pre balance statistics for each imputed dataset have now been saved in the 'balance/prebalance/' folder", "\n")
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
      if (!is.null(imp_conf)){
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(bal_thresh = ifelse(all_bal_stats$covariate %in% imp_conf,
                                                                             balance_thresh[1], balance_thresh[2]))
        all_bal_stats <- all_bal_stats %>%   dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) <
                                                                                all_bal_stats$bal_thresh, 1, 0) )
      }
      else{
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(bal_thresh = balance_thresh)
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) <
                                                                              all_bal_stats$bal_thresh, 1, 0) )
      }

      cat("\n")

      if (verbose){
        cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      }

      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }

   else if (is.data.frame(data)){
      if (sum(duplicated(data$"ID")) > 0){
        stop("Please provide wide dataset with a single row per ID.", call. = FALSE)
      }

      exposure_type <- ifelse(inherits(data[, paste0(exposure, '.', exposure_time_pts[1])], "numeric"), "continuous", "binary")
      bal_stats <- calcBalStats(data, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = 0, weights = NULL, imp_conf, verbose)

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
  else if (type == "weighted"){
    message(paste0("USER ALERT: The following statistics display covariate imbalance at each exposure time point following IPTW weighting,
            using ", form_name, " formulas."), "\n")
    cat("\n")

    if (mice::is.mids(data) | inherit(data, list)){
      # Running balance stats function, unweighted, on each imputed dataset w/ no weights
      if (mice::is.mids(data)){
        m <- data$m

        bal_stats <- lapply(seq_len(m), function(k) {
          d <- as.data.frame(mice::complete(data, k))

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputed datasets with a single row per ID.", call. = FALSE)
          }

          w <- weights[[k]]

          calcBalStats(d, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = k, weights = w, imp_conf, verbose) #are weights just in data at this point?
        })
      }


      else if (inherits(data, "list")){
        m = length(data)
        bal_stats <- lapply(seq_len(m), function(k) {
          d <- data[[k]]

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputed datasets with a single row per ID.", call. = FALSE)
          }

          w <- weights[[k]]

          calcBalStats(d, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = k, weights = w, imp_conf, verbose)
        })
      }


      # Save out balance stats for each imputed dataset
      bal_stats_all_imp <- do.call(bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      write.csv(bal_stats_all_imp, paste0(home_dir, "/balance/", type, "/", exposure,
                                          "_all_imps_balance_stat_summary.csv"), row.names = FALSE)

      if (verbose){
        message("Weighted balance statistics for each imputed dataset have now been saved in the 'balance/weighted/' folder", "\n")
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
      if (!is.null(imp_conf)){
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(bal_thresh = ifelse(all_bal_stats$covariate %in% imp_conf,
                                                                             balance_thresh[1], balance_thresh[2]))
        all_bal_stats <- all_bal_stats %>%dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) <
                                                                             all_bal_stats$bal_thresh, 1, 0) )
      }
      else{
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(bal_thresh = balance_thresh)
        all_bal_stats <- all_bal_stats %>% dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) <
                                                                              all_bal_stats$bal_thresh, 1, 0) )
      }

      cat("\n")
      if (verbose){
        cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      }

      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }


    if (is.data.frame(data)){
      #error checking
      if (sum(duplicated(data$"ID")) > 0){
        stop("Please provide wide dataset with a single row per ID.", call. = FALSE)
      }

      w <- weights[[1]]
      bal_stats <- calcBalStats(data, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = 0, weights = w, imp_conf, verbose)

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
  tot_cons <- tot_covars[!duplicated(tot_covars)] # Total domains/constructs

  # Make love plot to summarize imbalance at each exposure time point
  data_type <- ifelse(mice::is.mids(data) | inherits(data, "list"), "imputed", "single")

  weights_method = ifelse(type == "weighted", weights[[1]]$method, "no weights")

  sapply(seq_along(exposure_time_pts), function(i) {
    exposure_time_pt <- exposure_time_pts[i]

    if (mice::is.mids(data)){
      exposure_type <- ifelse(inherits(mice::complete(data,1)[, paste0(exposure, '.', exposure_time_pts[1])], "numeric"), "continuous", "binary")
    }
    else if (inherits(data, "list")){
      exposure_type <- ifelse(inherits(data[[1]][, paste0(exposure, '.', exposure_time_pts[1])], "numeric"), "continuous", "binary")
    }
    else if (is.data.frame(data)){
      exposure_type <- ifelse(inherits(data[, paste0(exposure, '.', exposure_time_pts[1])], "numeric"), "continuous", "binary")
    }

    temp <- all_bal_stats %>% filter(exp_time == exposure_time_pt)

    make_love_plot(home_dir, folder, exposure, exposure_time_pt, exposure_type, k = 0, form_name, temp,
                   data_type, balance_thresh, weights_method, imp_conf)
  })


  if (verbose){
    if (mice::is.mids(data) | inherits(data, "list")){
      cat(paste0("USER ALERT: Summary plots for ", form_name, " ", exposure,
                 " averaged across all imputations have been saved out for each time point in the 'balance/",
                 type, "/plots/' folder."), "\n")
    }
    else{
      message(paste0("USER ALERT: Summary plots for ", form_name, " ", exposure,
                     " have now been saved out for each time point in the 'balance/", type, "/plots/' folder."), "\n")
    }
  }


  # Save out all correlations/std mean differences
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_all_", type, "_", weights_method, "_associations.html"))
  stargazer::stargazer(all_bal_stats, type = "html", digits = 2, column.labels = colnames(all_bal_stats), summary = FALSE,
                       rownames = FALSE, header = FALSE, out = paste0(home_dir, "/balance/", type, "/", exposure, "-",
                                                                      outcome, "_all_", type,"_", weights_method, "_associations.html"))
  sink()

  if (verbose){
    if (mice::is.mids(data) | inherits(data, "list")){
      message(paste0("USER ALERT: Check 'balance/", type, "/' folder for a table of all ", type, " correlations or
                   standardized mean differences averaged across imputed datasets."), "\n")
    }
    else {
      message(paste0("USER ALERT: Check 'balance/", type, "/' folder for a table of all ", type, " correlations or
                   standardized mean differences."), "\n")
    }
  }


  # Saving out all pre-balance associations
  write.csv(all_bal_stats, paste0(home_dir, "/balance/", type, "/", exposure, "_", type, "_", weights_method, "_stat_summary.csv"), row.names = FALSE)

  if (verbose){
    message(paste0("All associations between exposure and covariates for ", form_name, " ", exposure, "-",
                   outcome, " have been saved in the '", type, "/' folder"), "\n")
    cat("\n")
  }

  # Finding all imbalanced variables
  unbalanced_covars <- all_bal_stats %>%
    filter(balanced == 0)
  unbalanced_constructs  <- sapply(strsplit(unbalanced_covars$covariate, "\\."),
                                   "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1))]


  if (mice::is.mids(data) | inherits(data, "list")){
    cat(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ",
               form_name, " formulas and ", weights_method, " :"), "\n")

    cat(paste0("The median absolute value relation between exposure and confounder is ", round(median(abs(all_bal_stats$avg_bal)), 2), " (range = ",
               round(min(all_bal_stats$avg_ba), 2), "-", round(max(all_bal_stats$avg_ba), 2), ")."), "\n")

    cat(paste0("As shown below, the following ", nrow(unbalanced_covars), " covariates across time points out of ",
               length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), "%) spanning ",
               length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
               "%) are imbalanced with a remaining median absolute value correlation/std mean difference in relation to ",
               exposure, " of ", round(median(abs(as.numeric(unlist(unbalanced_covars %>% dplyr:: filter(unbalanced_covars$balanced == 0) %>%
                                                                      dplyr:: select(avg_bal))))), 2), " (range=",
               round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "), "\n")

  }
  else {
    cat(paste0("USER ALERT: For exposure ", exposure, " using the ",form_name," formulas and ", weights_method, " :"), "\n")

    cat(paste0("The median absolute value relation between exposure and confounder is ", round(median(abs(all_bal_stats$avg_bal)), 2), " (range = ",
               round(min(all_bal_stats$avg_ba), 2), "-", round(max(all_bal_stats$avg_ba), 2), ")."), "\n")

    cat(paste0("As shown below, the following ", nrow(unbalanced_covars), " covariates across time points out of ",
               length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), "%) spanning ",
               length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
               "%) are imbalanced with a remaining median absolute value correlation/std mean difference in relation to ",
               exposure, " of ", round(median(abs(as.numeric(unlist(unbalanced_covars %>% dplyr:: filter(unbalanced_covars$balanced == 0) %>%
                                                                      dplyr:: select(avg_bal))))), 2), " (range=",
               round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "), "\n")
  }

  cat("\n")
  cat("\n")
  cat(knitr::kable(unbalanced_covars, caption = "Imbalanced Covariates", format = 'pipe'), sep = "\n")
  cat("\n")

  # Save out only imbalanced covariates
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_",type,"_", weights_method, "_all_covariates_imbalanced.html"))
  stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE, rownames = FALSE, header = FALSE,
            out = paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_", type, "_", weights_method, "_all_covariates_imbalanced.html"))
  sink()

  all_bal_stats
}



