#function to create love plots
make_love_plot <- function(home_dir, folder, exposure, exposure_time_pt, exposure_type, k = 0, form_name, balance_stats, data_type, balance_thresh, weights_method, imp_conf){

  stat_var <- colnames(balance_stats)[grepl("_bal", colnames(balance_stats))]
  colnames(balance_stats)[colnames(balance_stats) == stat_var] <- "avg_bal"

  x_lab <- ifelse(exposure_type == "continuous", "Correlation with Exposure", "Standardized Mean Difference Between Exposures")

  labels <- ifelse(balance_stats$balanced == 0, balance_stats$covariate, "")

  min_val <- ifelse(min(balance_stats[, "avg_bal"]) < 0, min(balance_stats[, "avg_bal"]) - 0.05, min(balance_thresh) - 0.05)
  max_val <- ifelse(max(balance_stats[, "avg_bal"]) > 0, max(balance_stats[, "avg_bal"]) + 0.05, max(balance_thresh) + 0.05)

  # Make love plot per exposure time point
  lp <- ggplot2::ggplot(aes(x = avg_bal, y = covariate), data = balance_stats) +
    ggplot2::geom_point(aes(y = as.factor(covariate), x = avg_bal, fill = "white", alpha = 1)) +
    ggplot2::geom_text(aes(label = labels, hjust = -0.2, vjust = 0.2), size = 1.5, color = "red") +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab("Covariate") +
    ggplot2::xlim(min_val, max_val) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                   axis.text.x = ggplot2::element_text(color = "black"),
                   axis.text.y = ggplot2::element_text(color = "black"),
                   axis.text = ggplot2::element_text(size = 8),
                   panel.border = ggplot2::element_rect(fill = NA, color = "black"),
                   plot.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10),
                   legend.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (nrow(balance_stats) > 40) { # stagger covariate labels if there are many
    lp <- lp + ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
  }



  if (!is.null(imp_conf)){ #adding threshold lines
    lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh[1], linetype = "dashed", color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh[1], linetype = "dashed", color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh[2], linetype = "dashed", color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh[2], linetype = "dashed", color = "red")
  } else{
    lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh, linetype = "dashed", color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh, linetype = "dashed", color = "red")

  }

  if (data_type == "imputed"){
    lp <- lp + ggplot2::ggtitle(paste0(exposure, " (t = ", exposure_time_pt, ") Balance for Imputation ", k))

    suppressMessages(ggplot2::ggsave(lp, filename = paste0(home_dir, "/balance/", folder, "plots/", form_name, "_imp_", k, "_", exposure, "_",
                                                           exposure_time_pt, "_", weights_method, "_summary_balance_plot.jpeg"), width = 6, height = 8))
  } else {
    lp <- lp + ggplot2::ggtitle(paste0(exposure, " (t = ", exposure_time_pt, ") Balance"))

    suppressMessages(ggplot2::ggsave(lp, filename = paste0(home_dir, "/balance/", folder, "plots/", form_name, "_", exposure, "_",
                                                           exposure_time_pt, "_", weights_method, "_summary_balance_plot.jpeg"), width = 6, height = 8))
  }
}
