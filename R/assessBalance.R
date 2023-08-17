#' Assess success of covariate balancing
#' Assesses how well balance was achieved for each of the covariates/potential confounds in relation to each of the exposure,
#' and returns a list of unbalanced covariates for each exposure to add to future models
#'
#' @param object msm object that contains all relevant user inputs
#' @param forms formula for assessing balance
#' @param data_for_model_with_weights imputed data with weights
#' @param histories optional binary indicator of whether to print histories for bal stats
#' @return list of unbalanced_covariates_for_models for each exposure
#' @seealso [msmObject()] for more on the weights_models param
#' @seealso [createWeights()] for more on the weights_models param
#' @export
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
#' @examples assessBalance(object, forms, data_for_model_with_weights, histories=1)

assessBalance <- function(home_dir, data, exposure,  outcome, tv_confounders, type, formulas, weights = NULL, balance_thresh = 0.1, imp_conf = NULL, user.o = T){

  # Error checking
  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.")
  }
  if(! type %in% c("prebalance", "weighted")){
    stop("Please provide a type from the following list: 'prebalance', 'weighted'")
  }
  if (!class(data) %in% c("mids", "data.frame", "character")) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.")
  }
  if (type == "prebalance" & !is.null(weights)){
    stop("The 'prebalance' mode of this function assesses balance prior to weighting and thus does not take weights.")
  }
  if (length(balance_thresh) == 2 & is.null(imp_conf)){
    stop("If you wish to provide different balance threshold for important and less important confounders, please provide a list of important confounders in the 'imp_conf' field.")
  }
  if (!is.null(imp_conf) & length(balance_thresh) == 1){
    stop("If you provide a list of important confounders, please provide a list of two balance tresholds for important and less important confounders, respectively")
  }

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

  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))

  exposure_type <- ifelse(class(data[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

  weights_method <- ifelse(type == "prebalance", "no weights", weights[[1]]$method)

  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[",1)


   if (class(data) == "character") {
    if (!dir.exists(data)) {
      stop("Please provide a valid directory path with imputed .csv's, a data frame, or a 'mids' object for the 'data' field.")
    }
    if (length(dir(data)) < 2) {
      stop("If you specify data as a directory, please supply more than 1 imputed dataset.")
    }

    # List imputed files
    files <- list.files(data, full.names = TRUE, pattern = "\\.csv")

    # Read and process imputed datasets as list
    data <- lapply(files, function(file) {
      imp_data <- read.csv(file)

      if (sum(duplicated(imp_data$"ID")) > 0){
        stop("Please provide a directory to data that are imputed in wide format with a single row per ID.")
      }
      imp_data
    })
  }


  ### Getting balance statistics
  if (type == "prebalance"){
    if (user.o == TRUE){
    message(paste0("USER ALERT: The following statistics display covariate imbalance at each exposure time point prior to weighting,
            using ", form_name, " formulas."), "\n")
    cat("\n")
    }

    # Running balance stats function, unweighted, on each imputed dataset
    if (class(data) == "mids" | class(data) =="list"){

      if (class(data) == "mids"){
        m = data$m

        bal_stats <- lapply(1:m, function(k) {
          d <-as.data.frame(mice::complete(data,k))

          exposure_type <- ifelse(class(d[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputed datasets with a single row per ID.")
          }

          calcBalStats(d, formulas, exposure, outcome, balance_thresh, k = k, weights = NULL, imp_conf, user.o)
        })
      }

      if (class(data) == "list"){
        m = length(data)

        bal_stats <- lapply(1:m, function(k) {
          d <- data[[k]]

          exposure_type <- ifelse(class(d[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputd datasets with a single row per ID.")
          }

          calcBalStats(d, formulas, exposure, outcome, balance_thresh, k = k, weights = NULL, imp_conf, user.o)
        })
      }


      # Save out balance stats for each imputed dataset
      bal_stats_all_imp <- do.call(bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      write.csv(bal_stats_all_imp, paste0(home_dir, "/balance/", type, "/", exposure, "_all_imps_balance_stat_summary.csv"), row.names = FALSE)

      if (user.o == TRUE){
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
        all_bal_stats <- all_bal_stats %>%   dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) < all_bal_stats$bal_thresh, 1, 0) )
      } else{
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(bal_thresh = balance_thresh)
        all_bal_stats <- all_bal_stats %>% dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) < all_bal_stats$bal_thresh, 1, 0) )
      }

      cat("\n")

      if (user.o == TRUE){
        cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      }

      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }


    if (class(data) == "data.frame"){
      if (sum(duplicated(data$"ID")) > 0){
        stop("Please provide wide dataset with a single row per ID.")
      }

      exposure_type <- ifelse(class(data[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

      bal_stats <- calcBalStats(data, formulas, exposure, outcome, balance_thresh, k = 0, weights = NULL, imp_conf, user.o)

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
  if (type == "weighted"){
    message(paste0("USER ALERT: The following statistics display covariate imbalance at each exposure time point following IPTW weighting,
            using ", form_name, " formulas."), "\n")
    cat("\n")

    if (class(data) == "mids" | class(data) =="list"){
      # Running balance stats function, unweighted, on each imputed dataset w/ no weights
      if (class(data) == "mids"){
        m <- data$m

        bal_stats <- lapply(1:m, function(k) {
          d <- as.data.frame(mice::complete(data, k))

          exposure_type <- ifelse(class(d[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputed datasets with a single row per ID.")
          }

          w <- weights[[k]]

          calcBalStats(d, formulas, exposure, outcome, balance_thresh, k = k, weights = w, imp_conf, user.o) #are weights just in data at this point?
        })
      }

      if (class(data) == "list"){
        m = length(data)
        bal_stats <- lapply(1:m, function(k) {
          d <- data[[k]]

          exposure_type <- ifelse(class(d[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

          if (sum(duplicated(d$"ID")) > 0){
            stop("Please provide wide imputed datasets with a single row per ID.")
          }

          w <- weights[[k]]

          calcBalStats(d, formulas, exposure, outcome, balance_thresh, k = k, weights = w, imp_conf, user.o)
        })
      }

      # Save out balance stats for each imputed dataset
      bal_stats_all_imp <- do.call(bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      write.csv(bal_stats_all_imp, paste0(home_dir, "/balance/", type, "/", exposure, "_all_imps_balance_stat_summary.csv"), row.names = FALSE)

      if (user.o == TRUE){
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
        all_bal_stats <- all_bal_stats %>%dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) < all_bal_stats$bal_thresh, 1, 0) )
      } else{
        all_bal_stats <- all_bal_stats %>% dplyr::mutate(bal_thresh = balance_thresh)
        all_bal_stats <- all_bal_stats %>% dplyr:: mutate(balanced = ifelse(abs(all_bal_stats$avg_bal) < all_bal_stats$bal_thresh, 1, 0) )
      }

      cat("\n")
      if (user.o == TRUE){
        cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      }

      tot_covars <- sapply(strsplit(all_bal_stats$covariate, "\\."), `[`, 1)
    }


    if (class(data) == "data.frame"){
      #error checking
      if (sum(duplicated(data$"ID")) > 0){
        stop("Please provide wide dataset with a single row per ID.")
      }

      w <- weights[[1]]

      bal_stats <- calcBalStats(data, formulas, exposure, outcome, balance_thresh, k = 0, weights = w, user.o)

      exposure_type <- ifelse(class(data[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

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
  # x_lab <- ifelse(exposure_type == "continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")


  # Make love plot to summarize imbalance at each exposure time point
  data_type <- ifelse(class(data) == "mids" | class(data) == "list", "imputed", "single")
  folder <- ifelse(type == "prebalance", "prebalance/", "weighted/")

  weights_method = ifelse(type == "weighted", w$method, "no weights")

  sapply(seq_along(exposure_time_pts), function(i) {
    exposure_time_pt <- exposure_time_pts[i]
    temp <- all_bal_stats %>% filter(exp_time == exposure_time_pt)

    # make_love_plot(home_dir, folder, exposure, exposure_type, temp, data_type, balance_thresh, imp_conf)
    make_love_plot(home_dir, folder, exposure, exposure_time_pt, exposure_type, k, form_name, temp, data_type, balance_thresh, weights_method, imp_conf)

  })

    # labels <- ifelse(temp$balanced == 0, temp$covariate, "")
    # min_val <- ifelse(min(temp$avg_bal) < 0, min(temp$avg_bal) - 0.02, min(balance_thresh) - 0.02)
    # max_val <- ifelse(max(temp$avg_bal) > 0, max(temp$avg_bal) + 0.02, max(balance_thresh) + 0.02)


    # lp <- ggplot2::ggplot(temp, aes(x = avg_bal, y = covariate)) +
    #   ggplot2::geom_point(aes(fill = "white", alpha = 1), na.rm = TRUE) +
    #   ggplot2::geom_text(aes(label = labels, hjust = -0.2, vjust = 0.2), size = 1.5, color = "red") +
    #   ggplot2::xlab(x_lab) +
    #   ggplot2::ylab("Covariate") +
    #   ggplot2::xlim(min_val, max_val) +
    #   ggplot2::ggtitle(paste0(exposure, "(t=", exposure_time_pt, ") Balance")) +
    #   ggplot2::theme(
    #     panel.background = ggplot2::element_rect(fill = "white"),
    #     axis.text.x = ggplot2::element_text(color = "black"),
    #     axis.text.y = ggplot2::element_text(color = "black"),
    #     axis.text = ggplot2::element_text(size = 8),
    #     panel.border = ggplot2::element_rect(fill = NA, color = "black"),
    #     plot.background = ggplot2::element_blank(),
    #     plot.title = ggplot2::element_text(size = 10),
    #     legend.background = ggplot2::element_blank(),
    #     legend.key = ggplot2::element_blank(),
    #     legend.position = "none"
    #   ) +
    #   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    #
    # if (nrow(temp) > 40) { # Stagger covariate labels if there are many to ease viewing
    #   lp <- lp + ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
    # }
    #
    # if (!is.null(imp_conf)){ #adding threshold lines
    #   lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh[1], linetype = "dashed", color = "red")
    #   lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh[1], linetype = "dashed", color = "red")
    #   lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh[2], linetype = "dashed", color = "red")
    #   lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh[2], linetype = "dashed", color = "red")
    # } else{
    #   lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh, linetype = "dashed", color = "red")
    #   lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh, linetype = "dashed", color = "red")
    # }
#
#     if (class(data) == "mids" | class(data) == "list"){
#       suppressMessages(ggplot2::ggsave(lp, filename = paste0(home_dir, "/balance/", type, "/plots/", form_name, "_", weights_method, "_", exposure, "_all_imps_",
#                                                              exposure_time_pt, "_summary_", type, "_plot.jpeg")))
#     } else{
#       suppressMessages(ggplot2::ggsave(lp, filename = paste0(home_dir, "/balance/", type, "/plots/", form_name, "_", weights_method, "_", exposure, "_",
#                                                              exposure_time_pt, "_summary_", type, "_plot.jpeg")))
#     }
#   })

  if (user.o == TRUE){
    if (class(data) == "mids" | class(data) == "list"){
      cat(paste0("USER ALERT: Summary plots for ", form_name, " ", exposure,
                 " averaged across all imputations have been saved out for each time point in the 'balance/", type, "/plots/' folder."), "\n")
    } else{

      message(paste0("USER ALERT: Summary plots for ", form_name, " ", exposure, " have now been saved out for each time point in the 'balance/", type, "/plots/' folder."), "\n")
    }
  }


  # Save out all correlations/std mean differences
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_all_", type, "_", weights_method, "_associations.html"))
  stargazer::stargazer(all_bal_stats, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE,
                       rownames = FALSE, header = FALSE, out = paste0(home_dir, "/balance/", type, "/", exposure, "-",
                                                                      outcome, "_all_", type,"_", weights_method, "_assocations.html"))
  sink()

  if (user.o == TRUE){
    if (class(data) == "mids" | class(data) == "list"){
      message(paste0("USER ALERT: Check 'balance/", type, "/' folder for a table of all ", type, " correlations or
                   standardized mean differences averaged across imputed datasets."), "\n")
    } else {
      message(paste0("USER ALERT: Check 'balance/", type, "/' folder for a table of all ", type, " correlations or
                   standardized mean differences."), "\n")
    }
  }

  # <<<<<<< Updated upstream
  #   # Renames factor covariates
  #   unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates] <-
  #     sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates]
  #   unbalanced_constructs <- sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1))]
  # =======

  # Save out all correlations/std mean differences
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_all_", type, "_", weights_method, "_associations.html"))
  stargazer(all_bal_stats, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE,
            rownames = FALSE, header = FALSE, out = paste0(home_dir, "/balance/", type, "/", exposure, "-",
                                                           outcome, "_all_", type,"_", weights_method, "_assocations.html"))
  sink()

  # Saving out all pre-balance associations
  write.csv(all_bal_stats, paste0(home_dir, "/balance/", type, "/", exposure, "_", type, "_", weights_method, "_stat_summary.csv"), row.names = FALSE)

  if (user.o == TRUE){
    message(paste0("All associations between exposure and covariates for ", form_name, " ", exposure, "-",
                   outcome, " have been saved in the '", type, "/' folder"), "\n")
    cat("\n")
  }

  # Finding all imbalanced variables\
  unbalanced_covars <- all_bal_stats %>%
    filter(balanced == 0)

  #   if (class(data) == "mids" | class(data) == "list"){
  #     cat(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ", form_name, " formulas and ", weights_method, " :"), "\n")
  #     cat(paste0("The median absolute value relation between exposure and confounder is ", round(median(abs(all_bal_stats$avg_bal)), 2), " (range = ",
  #                round(min(all_bal_stats$avg_ba), 2), "-", round(max(all_bal_stats$avg_ba), 2), ")."), "\n")
  #     cat("\n")
  #
  #     message(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ",
  #                    form_name, ", the following ", nrow(unbalanced_covars), " covariates across time points out of ",
  #                    length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), ",%) spanning ",
  #                    length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
  #                    "%) are imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
  #                    exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)), 2), " (range=",
  #                    round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "), "\n")
  #   }else {
  #     cat(paste0("USER ALERT: For exposure ", exposure, " using the ", form_name, " formulas and ", weights_method, " :"), "\n")
  #     cat(paste0("The median absolute value relation between exposure and confounder is ", round(median(abs(all_bal_stats$avg_bal)), 2), " (range = ",
  #                round(min(all_bal_stats$avg_ba), 2), "-", round(max(all_bal_stats$avg_ba), 2), ")."), "\n")
  #     cat("\n")
  #
  #     message(paste0("USER ALERT: For exposure ", exposure, " using the ",
  #                    form_name, ", the following ", nrow(unbalanced_covars), " covariates across time points out of ",
  #                    length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), ",%) spanning ",
  #                    length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2), "%) are imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
  #                    exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)), 2), " (range=",
  #                    round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "))
  #   }
  # # =======

  # =======
  #   unbalanced_covars <- all_bal_stats %>%
  #     filter(balanced_avg == 0)
  #   cat("\n")
  #   unbalanced_constructs <- sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1))]


  if (class(data) == "mids" | class(data) == "list"){
    cat(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ", form_name, " formulas and ", weights_method, " :"), "\n")
    cat(paste0("The median absolute value relation between exposure and confounder is ", round(median(abs(all_bal_stats$avg_bal)), 2), " (range = ",
               round(min(all_bal_stats$avg_ba), 2), "-", round(max(all_bal_statss$avg_ba), 2), ")."), "\n")
    cat(paste0("As shown below, the following ", nrow(unbalanced_covars), " covariates across time points out of ",
               length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), "%) spanning ",
               length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
               "%) are imbalanced with a remaining median absolute value correlation/std mean difference in relation to ",
               exposure, " of ", round(median(abs(as.numeric(unlist(unbalanced_covars %>% dplyr:: filter(unbalanced_covars$balanced == 0) %>%
                                                                      dplyr:: select(avg_bal))))), 2), " (range=",
               round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "), "\n")
  }else {
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
  # >>>>>>> Stashed changes
  # =======
  #                exposure, " of ", round(median(abs(as.numeric(unlist(unbalanced_covars %>% dplyr:: filter(unbalanced_covars$balanced_avg == 0) %>%
  #                                                                   dplyr:: select(avg_bal))))), 2), " (range=",
  #                round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "), "\n")

  cat("\n")
  cat(knitr::kable(unbalanced_covars, caption = "Imbalanced Covariates", format = 'pipe'), sep = "\n")

  # Save out only imbalanced covariates
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_",type,"_", weights_method, "_all_covariates_imbalanced.html"))
  stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE, rownames = FALSE, header = FALSE,
            out = paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_", type, "_", weights_method, "_all_covariates_imbalanced.html"))
  sink()

  all_bal_stats

}



