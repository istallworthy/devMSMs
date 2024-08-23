#' Create love plots showing balancing statistics
#'
#' @param balance_stats data frame of balance statistics
#' @param exposure_type String either "continuous" or "binary"
#' @param k imputation number
#' @param weight_method name of weighting method (NULL for unweighted)
#' @return ggplot2 object
#' @examples 
#' print("tbd")
#'
#' @noRd
make_love_plot <- function(obj, balance_stats, k = 0, t, weight_method = NULL) {
  
  exposure <- obj[["exposure"]]
  exposure_type <- obj[["exposure_type"]]
  exposure_time <- obj[["exposure_time_pts"]]
  # exposure <- balance_stats$exposure[1]
  # exposure_time_pt <- balance_stats$exposure_time[1]
  # exposure_type <- balance_stats$exposure_type[1]
  
  # Symmetric
  balance_thresh <- unique(balance_stats$bal_thresh)
  xrange <- c(
    balance_stats$std_bal_stats,
    balance_thresh,
    -1 * balance_thresh
  )
  xrange <- c(-max(abs(xrange)), max(abs(xrange)))
  
  # timing title
  if (length(t) > 1){
    t_t <- " For all Exposure Timepoints"
  } else {
    t_t <- sprintf(" For Exposure Timepoint %s", exposure_time[t])
  }
  
  # data.frame, list, mice
  if (is.na(k)) {
    title_imputation_note <- sprintf(" Averaging Across Imputed Datasets")
  } else if (k == 0) {
    title_imputation_note <- ""
  } else {
    title_imputation_note <- sprintf(" for Imputation %s", k)
  }
  
  if (is.null(weight_method)) {
    title <- sprintf("Covariate Balance%s%s", title_imputation_note, t_t)
  } else {
    title <- sprintf("Covariate Balance%s%s using `%s` weights ", title_imputation_note, t_t, weight_method)
  }
  
  x_lab <- if (exposure_type == "continuous") {
    "Correlation with Exposure"
  } else {
    "Standardized Mean Difference Between Exposures"
  }
  
  y_scale <- ggplot2::scale_y_discrete()
  if (nrow(balance_stats) > 40) { # stagger covariate labels if there are many
    y_scale <- ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
  }
  
  # For correct y-axis ordering
  balance_stats <- balance_stats[order(balance_stats$std_bal_stats), , drop = FALSE]
  balance_stats$covariate <- factor(
    balance_stats$covariate, 
    levels = unique(balance_stats$covariate)
  )
  names(balance_stats)[names(balance_stats) == "exposure"] <- "Exposure"

  balance_stats$Exposure <- factor(
    balance_stats$Exposure, 
    levels = exposure
  )
  
  # Make love plot per exposure time point
  lp <- ggplot2::ggplot(data = balance_stats) +
    ggplot2::geom_vline(
      xintercept = c(balance_thresh, -balance_thresh),
      linetype = "dashed", color = "red"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "red"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        y = .data$covariate,
        x = .data$std_bal_stats,
        color = factor(.data$balanced == 0), #added 
        fill = factor(.data$balanced == 0)  # added
      ),
      # fill = "white", alpha = 1, size = 1.2
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +  
    ggplot2::facet_wrap(
      ~ .data$Exposure, 
      scales = "free_y",
      labeller = "label_both"
    ) +
    ggplot2::labs(x = x_lab, y = "Covariate", title = title) + 
    ggplot2::scale_x_continuous(limits = xrange) +
    y_scale +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(size = 8),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "black"
      ),
      plot.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.position = "none"
    ) + # adding imbalanced labels
    ggplot2::geom_text(
      data = balance_stats[balance_stats$balanced == 0,],
      ggplot2::aes(
        y = .data$covariate,
        x = .data$std_bal_stats,
        label = .data$covariate
      ),
      color = "red",
      vjust = 1,  
      hjust = -0.2,
      size = 1.5    
    )
  
  return(lp)
}
