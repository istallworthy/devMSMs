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
assessBalance <- function(data, exposure, outcome, tv_confounders, type, formulas, weights = NULL, balance_thresh = 0.1){

  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))
  ID <- "S_ID"

  if(! type %in% c("prebalance", "weighted")){
    stop("Please provide a type from the following list: 'prebalance', 'weighted'")
  }

  weights_method <- ifelse(type == "prebalance", "no weights", weights$method)

  # Getting form type from input
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[",1)

  exposure_type <- ifelse(class(data[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")


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

  ### Assessing data type
  # Check data type
  if (!class(data) %in% c("mids", "data.frame", "character")) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.")
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
      # imp_data <- imp_data %>%
      #   dplyr::select(-c(X, contains("ESETA1_")))
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
    # Combine imputed datasets
    imp2 <- do.call(rbind.data.frame, imps2)
    imp2$X <- seq_len(nrow(imp2))
    data <- mice::as.mids(imp2, .imp = ".imp")
  }


  ### Getting balance statistics

  if (type == "prebalance"){

    message(paste0("USER ALERT: The following statistics display covariate imbalance at each exposure time point prior to weighting,
            using ", form_name, "formulas."), "\n")
    cat("\n")

    # Running balance stats function, unweighted, on each imputed dataset w/ no weights
    if (class(data) == "mids"){
      m=data$m
      bal_stats <- lapply(1:m, function(k) {
        d <- data[[k]]
        calcBalStats(d, formulas, exposure, outcome, balance_thresh, k = k, weights = NULL)
      })

      # Save out balance stats for each imputed dataset
      bal_stats_all_imp <- do.call(bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      write.csv(bal_stats_all_imp, paste0(home_dir, "/balance/", type, "/", exposure, "_all_imps_balance_stat_summary.csv"), row.names = FALSE)
      message("Pre balance statistics for each imputed dataset have now been saved in the 'balance/prebalance/' folder", "\n")

      # Gathering imbalanced covariate statistics to average across imputed datasets for the final list/assessment of imbalanced covariates
      # Averaging across imputed datasets
      unbalanced_covars <- data.frame(
        exposure = exposure,
        exp_time = bal_stats[[1]]$exp_time,
        covar_time = bal_stats[[1]]$covar_time,
        covariate = bal_stats[[1]]$covariate,
        avg_bal = rowMeans(do.call(cbind, lapply(bal_stats, `[[`, "std_bal_stats")))
      ) %>%
        mutate(balanced_avg = ifelse(abs(avg_bal) < balance_thresh, 1, 0)) # Assessing new averaged bal stat in relation to the balance threshold


      cat("\n")
      cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      # Getting totals
      tot_covars <- sapply(strsplit(bal_stats[[1]]$covariate, "\\."), `[`, 1)
    }

    if (class(data) == "data.frame"){
      bal_stats <- calcBalStats(data, formulas, exposure, outcome, balance_thresh, k = 0, weights = NULL)

      # Gathering imbalanced covariate statistics for the final list/assessment of imbalanced covariates
      unbalanced_covars <- data.frame(
        exposure = exposure,
        exp_time = bal_stats$exp_time,
        covar_time = bal_stats$covar_time,
        covariate = bal_stats$covariate,
        avg_bal = bal_stats$std_bal_stats,
        balanced_avg = bal_stats$balanced)

      # Getting totals
      tot_covars <- sapply(strsplit(bal_stats$covariate, "\\."), `[`, 1)
    }
  } #ends prebalance


  if (type == "weighted"){

    message(paste0("USER ALERT: The following statistics display covariate imbalance at each exposure time point following IPTW weighting,
            using ", form_name, "formulas."), "\n")
    cat("\n")

    # Running balance stats function, unweighted, on each imputed dataset w/ no weights
    if (class(data) == "mids"){
      m <- data$m
      bal_stats <- lapply(1:m, function(k) {
        d <- data[[k]]
        w <- weights[[k]]
        calcBalStats(d, formulas, exposure, outcome, balance_thresh, k = k, weights = w) #are weights just in data at this point?
      })

      # Save out balance stats for each imputed dataset
      bal_stats_all_imp <- do.call(bind_rows, bal_stats)
      bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
      write.csv(bal_stats_all_imp, paste0(home_dir, "/balance/", type, "/", exposure, "_all_imps_balance_stat_summary.csv"), row.names = FALSE)
      message("Weighted balance statistics for each imputed dataset have now been saved in the 'balance/weighted/' folder", "\n")

      # Gathering imbalanced covariate statistics to average across imputed datasets for the final list/assessment of imbalanced covariates
      # Averaging across imputed datasets
      unbalanced_covars <- data.frame(
        exposure = exposure,
        exp_time = bal_stats[[1]]$exp_time,
        covar_time = bal_stats[[1]]$covar_time,
        covariate = bal_stats[[1]]$covariate,
        avg_bal = rowMeans(do.call(cbind, lapply(bal_stats, `[[`, "std_bal_stats")))
      ) %>%
        mutate(balanced_avg = ifelse(abs(avg_bal) < balance_thresh, 1, 0)) # Assessing new averaged bal stat in relation to the balance threshold


      cat("\n")
      cat(paste0("*** Averaging Across All Imputations ***"), "\n")
      # Getting totals
      tot_covars <- sapply(strsplit(bal_stats[[1]]$covariate, "\\."), `[`, 1)
    }

    if (class(data) == "data.frame"){

      bal_stats <- calcBalStats(data, formulas, exposure, outcome, balance_thresh, k = 0, weights = weights)

      # Gathering imbalanced covariate statistics for the final list/assessment of imbalanced covariates
      unbalanced_covars <- data.frame(
        exposure = exposure,
        exp_time = bal_stats$exp_time,
        covar_time = bal_stats$covar_time,
        covariate = bal_stats$covariate,
        avg_bal = bal_stats$std_bal_stats,
        balanced_avg = bal_stats$balanced)

      # Getting totals
      tot_covars <- sapply(strsplit(bal_stats$covariate, "\\."), `[`, 1)
    }
  } #ends weighted




  ### Plotting and summarizing

  tot_cons <- tot_covars[!duplicated(tot_covars)] # Total domains/constructs
  x_lab <- ifelse(exposure_type == "continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")

  # Make love plot to summarize imbalance at each exposure time point
  sapply(seq_along(exposure_time_pts), function(i) {
    exposure_time_pt <- exposure_time_pts[i]
    temp <- unbalanced_covars %>% filter(exp_time == exposure_time_pt)
    labels <- ifelse(temp$balanced_avg == 0, temp$covariate, "")
    min_val <- ifelse(min(temp$avg_bal) < 0, min(temp$avg_bal) - 0.1, balance_thresh - 0.05)
    max_val <- ifelse(max(temp$avg_bals) > 0, max(temp$avg_bal) + 0.1, balance_thresh + 0.05)

    lp <- ggplot(temp, aes(x = avg_bal, y = covariate)) +
      geom_point(aes(fill = "white", alpha = 1), na.rm = TRUE) +
      geom_text(aes(label = labels, hjust = -0.2, vjust = 0.2), size = 1.5, color = "red") +
      geom_vline(xintercept = balance_thresh, linetype = "dashed", color = "red") +
      geom_vline(xintercept = -balance_thresh, linetype = "dashed", color = "red") +
      xlab(x_lab) +
      ylab("Covariate") +
      xlim(min_val, max_val) +
      ggtitle(paste0(exposure, "(t=", exposure_time_pt, ") Balance")) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.background = element_blank(),
        plot.title = element_text(size = 10),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "none"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
    if (nrow(temp) > 40) { # Stagger covariate labels if there are many to ease viewing
      lp <- lp + scale_y_discrete(guide = guide_axis(n.dodge = 2))
    }
    suppressMessages(ggsave(lp, filename = paste0(home_dir, "/balance/", type, "/plots/", form_name, "_", weights_method, "_", exposure, "_all_imps_",
                                                  exposure_time_pt, "_summary_", type, "_plot.jpeg")))
    if (class(data) == "mids"){
      message(paste0("USER ALERT: A balance summary plot for ", form_name, " ", exposure, " averaged across all imputations at time ",
                     exposure_time_pt, " has now been saved in the 'balance/", type, "/plots/' folder."), "\n")
    } else{
      message(paste0("USER ALERT: A balance summary plot for ", form_name, " ", exposure, "  at time ",
                     exposure_time_pt, " has now been saved in the 'balance/", type, "/plots/' folder."), "\n")
    }
  })

  # Save out all correlations/std mean differences
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_all_", type, "_", weights_method, "_associations.html"))
  stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE,
            rownames = FALSE, header = FALSE, out = paste0(home_dir, "/balance/", type, "/", exposure, "-",
                                                           outcome, "_all_", type,"_", weights_method, "_assocations.html"))
  sink()

  if (class(data) == "mids"){
    message(paste0("USER ALERT: Check 'balance/", type, "/' folder for a table of all ", type, " correlations or
                   standardized mean differences averaged across imputed datasets."), "\n")
  } else {
    message(paste0("USER ALERT: Check 'balance/", type, "/' folder for a table of all ", type, " correlations or
                   standardized mean differences."), "\n")
  }

  # Renames factor covariates
  unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates] <-
    sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates]
  unbalanced_constructs <- sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1))]

  # Saving out all pre-balance associations
  write.csv(unbalanced_covars, paste0(home_dir, "/balance/", type, "/", exposure, "_", type, "_", weights_method, "_stat_summary.csv"), row.names = FALSE)
  message(paste0("All associations between exposure and covariates for ", form_name, " ", exposure, "-",
                 outcome, " have been saved in the '", type, "/' folder"), "\n")
  cat("\n")

  # Finding all imbalanced variables
  unbalanced_covars <- unbalanced_covars %>%
    filter(balanced_avg == 0)
  cat("\n")

  if (class(data) == "mids"){
    message(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ",
                   form_name, ", the following ", nrow(unbalanced_covars), " covariates across time points out of ",
                   length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), ",%) spanning ",
                   length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2),
                   "%) are imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
                   exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)), 2), " (range=",
                   round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "), "\n")
  }else {
    message(paste0("USER ALERT: For exposure ", exposure, " using the ",
                   form_name, ", the following ", nrow(unbalanced_covars), " covariates across time points out of ",
                   length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), ",%) spanning ",
                   length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2), "%) are imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
                   exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)), 2), " (range=",
                   round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "))
  }
  cat(knitr::kable(unbalanced_covars, caption = "Imbalanced Covariates", format = 'pipe'), sep = "\n")

  # Save out only imbalanced covariates
  sink(paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_",type,"_", weights_method, "_all_imbalanced_covariates.html"))
  stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE, rownames = FALSE, header = FALSE,
            out = paste0(home_dir, "/balance/", type, "/", exposure, "-", outcome, "_", type, "_", weights_method, "_all_imbalanced_covariates.html"))
  sink()

  unbalanced_covars

}



#
# # Comparing balance to prebalance stats if using full balance forms (which were used for prebalance assessment)
# if (form_type == "full_forms") {
#   prebal <- read.csv(paste0(home_dir, "pre balance/", exposure, "_prebalance_stat_summary.csv"), row.names = NULL)
#   colnames(prebal)[colnames(prebal) == "avg_bal"] <- "avg_pre_bal"
#   colnames(prebal)[colnames(prebal) == "balanced_avg"] <- "pre_balanced_avg"
#   comp <- merge(unbalanced_covars, prebal, by = c("exposure", "exp_time", "covar_time", "covariate"))
#   cat("\n")
#   cat(paste0("Prior to weighting using the full balancing formula, ", sum(comp$pre_balanced_avg), " out of ", nrow(comp), " (", round((sum(comp$pre_balanced_avg) / nrow(comp)) * 100, 2), "%) covariates were balanced compared to ",
#              sum(comp$balanced_avg), " out of ", nrow(comp), " (", round((sum(comp$balanced_avg) / nrow(comp)) * 100, 2), "%) after weighting. ",
#              sum(comp$balanced_avg - comp$pre_balanced_avg == -1), " covariates were balanced prior to weighting and are now imbalanced after weighting:"
#   ), "\n")
#   cat(knitr::kable(comp[comp$balanced_avg - comp$pre_balanced_avg == -1, ], format = 'pipe'), sep = "\n")
# }
#
# # Getting just imbalanced covars
# unbalanced_covars <- unbalanced_covars %>% filter(balanced_avg == 0)
#
# # Renames factor covariates
# unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates] <- sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates]
# unbalanced_constructs <- sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1))]
#
# cat("\n")
# cat(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ", weights_method, " weighting method and the ", form_name,
#            ", the following ", nrow(unbalanced_covars), " covariates across time points out of ",
#            length(tot_covars), " total (", round((nrow(unbalanced_covars) / length(tot_covars)) * 100, 2), "%) spanning ",
#            length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round((length(unbalanced_constructs) / length(tot_cons)) * 100, 2),
#            " %) remain imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
#            exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)), 2), " (range=",
#            round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), "): "
# ), "\n")
# cat("\n")
# cat(knitr::kable(unbalanced_covars, caption = paste0("Imbalanced Covariates using ", weights_method, " and ", form_name), format = 'pipe'), sep = "\n")
# cat("\n")
#
# # Save out only imbalanced covariates
# sink(paste0(home_dir, "balance/comparison values/", form_name, "_", weights_method, "_all_imbalanced_covariates.html"))
# stargazer::stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE, rownames = FALSE, header = FALSE,
#                      out = paste0(home_dir, "balance/comparison values/", "_", weights_method, "_", form_name, "_all_imbalanced_covariates.html")
# )
# sink()
#
# write.csv(unbalanced_covars, paste0(home_dir, "balance/comparison values/", form_name, "_", exposure, "_unbalanced_stat_summary.csv"))
# cat(paste0("Imbalanced covariates for ", exposure, " using the ", form_name, " and ", weights_method, " weights method have been saved in the 'balance/comparison values/' folder"), "\n")
#
# return(bal_stats)


