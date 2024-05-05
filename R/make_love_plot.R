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
#' @export
make_love_plot <- function(balance_stats, exposure_type = c("continuous", "binary"), k = 0, weight_method = NULL) {
  exposure <- balance_stats$exposure[1]
  exposure_time_pt <- balance_stats$exposure_time[1]
  # exposure_type <- balance_stats$exposure_type[1]

  balance_thresh <- unique(balance_stats$bal_thresh)

  # Symmetric
  xrange <- c(
    balance_stats$std_bal_stats,
    balance_thresh,
    -1 * balance_thresh
  )
  xrange <- c(-max(abs(xrange)), max(abs(xrange)))

  # TODO: Check about `data_type` = "single"
  # data.frame, list, mice
  if (k == 0) {
    if (is.null(weight_method)) {
      title <- sprintf("%s Balance", exposure)
    } else {
      title <- sprintf("%s Balance using `%s` weights", exposure, weight_method)
    }
  } else {
    if (is.null(weight_method)) {
      title <- sprintf("%s Balance for Imputation %s", exposure, k)
    } else {
      title <- sprintf("%s Balance for Imputation %s using `%s` weights", exposure, k, weight_method)
    }
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
  balance_stats$covariate <- factor(balance_stats$covariate, levels = balance_stats$covariate, ordered = TRUE)

  # Make love plot per exposure time point
  lp <- ggplot2::ggplot(data = balance_stats) +
    ggplot2::geom_point(
      ggplot2::aes(
        y = .data$covariate,
        x = .data$std_bal_stats
      ),
      fill = "white", alpha = 1
    ) +
    ggplot2::geom_vline(
      xintercept = c(balance_thresh, -balance_thresh),
      linetype = "dashed", color = "red"
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
    )

  return(lp)
}
