#' Examining initial imbalance for each exposure, each time point, each covariate for each history using the Jackson method
#'
#' @param object msm object that contains all relevant user inputs
#' @param wide_long_datasets from formatForWeights()
#' @param all_forms from createForms()
#' @param histories binary indicator of whether to print history sample distributions
#' @return list of unbalanced_covariates_for_models for each exposure
#' @export
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
#' @importFrom stargazer stargazer
#' @importFrom dplyr bind_rows
#' @examples preBalanceChecking(object, wide_long_datasets, all_forms, histories = 1)
preBalanceChecking <- function(object, wide_long_datasets, all_forms, histories = 1) {
  exposure <- object$exposure
  outcome <- object$outcome
  m <- object$m
  exposure_time_pts <- object$exposure_time_pts
  balance_thresh <- object$balance_thresh
  factor_covariates <- object$factor_covariates
  home_dir <- object$home_dir

  # Getting form type from input
  form_type <- match.call()[[4]]
  form_name <- as.character(form_type)

  exposure_type <- ifelse(length(unique(wide_long_datasets[[1]][, paste0(exposure, ".", exposure_time_pts[1])])) < 3, "binary", "continuous")

  cat("USER ALERT: The following statistics display covariate imbalance at each exposure time point prior to weighting, using full formulas that reflect covariates at all lagged time points.", "\n")
  cat("\n")

  # Running balance stats function, unweighted, on each imputed dataset
  bal_stats <- lapply(1:m, function(k) {
    calcBalStats(object, wide_long_datasets, all_forms, form_name, exposure, outcome, k = k, weighted = 0, histories)
  })

  # Save out balance stats for each imputed dataset
  bal_stats_all_imp <- do.call(bind_rows, bal_stats)
  bal_stats_all_imp <- bal_stats_all_imp[order(bal_stats_all_imp$covariate), ]
  write.csv(bal_stats_all_imp, paste0(home_dir, "pre balance/", exposure, "_all_imps_balance_stat_summary.csv"), row.names = FALSE)
  cat("Pre balance statistics for each imputed dataset have now been saved in the 'pre balance/' folder", "\n")

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

  # Getting totals
  tot_covars <- sapply(strsplit(bal_stats[[1]]$covariate, "\\."), `[`, 1)
  tot_cons <- tot_covars[!duplicated(tot_covars)] # Total domains/constructs
  x_lab <- ifelse(exposure_type == "continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")

  cat("\n")
  cat(paste0("*** Averaging Across All Imputations ***"), "\n")

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
    suppressMessages(ggsave(lp, filename = paste0(home_dir, "pre balance/plots/", form_name, "_", exposure, "_all_imps_", exposure_time_pt, "_summary_pre_balance_plot.jpeg")))
    cat(paste0("USER ALERT: A balance summary plot for ", form_name, " ", exposure, " averaged across all imputations at time ", exposure_time_pt, " has now been saved in the '", "pre balance/plots/' folder."), "\n")
  })

  # Save out all pre-balance correlations/std mean differences
  sink(paste0(home_dir, "pre balance/", exposure, "-", outcome, "_all_pre-balance_associations.html"))
  stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE, rownames = FALSE, header = FALSE,
            out = paste0(home_dir, "pre balance/", exposure, "-", outcome, "_all_pre-balance_assocations.html"))
  sink()
  cat(paste0("USER ALERT: Check 'pre balance/comparison values/' folder for a table of all pre-balance correlations or standardized mean differences and averaged across imputed datasets."), "\n")

  # Renames factor covariates
  unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates] <-
    sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[", 1) %in% factor_covariates]
  unbalanced_constructs <- sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[", 1))]

  # Saving out all pre-balance associations
  write.csv(unbalanced_covars, paste0(home_dir, "pre balance/", exposure, "_prebalance_stat_summary.csv"), row.names = FALSE)
  cat(paste0("All associations between exposure and covariates for ", form_name, " ", exposure, "-", outcome, " prior to balancing have been saved in the 'pre balance/' folder"), "\n")
  cat("\n")

  # Finding all imbalanced variables
  unbalanced_covars <- unbalanced_covars %>%
    filter(balanced_avg == 0)
  cat("\n")

  cat(paste0("USER ALERT: Before weighting, averaging across all imputed datasets for exposure ", exposure, " using the ", form_name, ", the following ", nrow(unbalanced_covars), " covariates across time points out of ",
             length(tot_covars), " total (", round(nrow(unbalanced_covars) / length(tot_covars) * 100, 2), ",%) spanning ",
             length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs) / length(tot_cons) * 100, 2), "%) are imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
             exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)), 2), " (range=",
             round(min(unbalanced_covars$avg_bal), 2), "-", round(max(unbalanced_covars$avg_bal), 2), ") : "
  ), "\n")
  cat(knitr::kable(unbalanced_covars, caption = "Imbalanced Covariates Before Weighting", format = 'pipe'), sep = "\n")

  # Save out only imbalanced covariates
  sink(paste0(home_dir, "pre balance/", exposure, "-", outcome, "_all_imbalanced_covariates.html"))
  stargazer(unbalanced_covars, type = "html", digits = 2, column.labels = colnames(unbalanced_covars), summary = FALSE, rownames = FALSE, header = FALSE,
            out = paste0(home_dir, "pre balance/", exposure, "-", outcome, "_all_imbalanced_covariates.html")
  )
  sink()
}
