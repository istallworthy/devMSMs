#' Create love plots showing balancing statistics
#'
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param folder folder path for saving
#' @param exposure name of exposure variable
#' @param exposure_time_pt exposure time point integer
#' @param exposure_type character indicating binary or continuous exposure type
#' @param k imputation number
#' @param form_name formula name
#' @param balance_stats data frame of balance statistics
#' @param data_type single or imputed data type
#' @param balance_thresh one or two numbers between 0 and 1 indicating a single
#'   balancing threshold or thresholds for more and less important confounders,
#'   respectively
#' @param weights_method method character string of WeightItMSM() balancing
#'   method abbreviation
#' @param imp_conf list of variable names reflecting important confounders
#' @param verbose TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out TRUE or FALSE indicator to save output and intermediary
#'   output locally
#' @return none
#' @export
#' @examples
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
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
#'
#' p <- make_love_plot(folder = "prebalance/",
#'                     exposure = "A",
#'                     exposure_time_pt = 1,
#'                     exposure_type = "continuous",
#'                     form_name = "form_name",
#'                     balance_stats = b,
#'                     data_type = "single",
#'                     balance_thresh = 0.1,
#'                     imp_conf = NULL,
#'                     weights_method = w[[1]]$method,
#'                     save.out = FALSE,
#'                     verbose = TRUE)
#' p <- make_love_plot(folder = "weighted/",
#'                     exposure = "A",
#'                     exposure_time_pt = 2,
#'                     exposure_type = "continuous",
#'                     form_name = "form_name",
#'                     balance_stats = b,
#'                     data_type = "single",
#'                     balance_thresh = c(0.05, 0.1),
#'                     imp_conf = "A.2",
#'                     weights_method = w[[1]]$method,
#'                     save.out = FALSE,
#'                     verbose = TRUE)


make_love_plot <- function(home_dir, folder, exposure, exposure_time_pt, exposure_type, k = 0, form_name, balance_stats, data_type,
                           balance_thresh, weights_method, imp_conf, verbose, save.out) {

  stat_var <- colnames(balance_stats)[grepl("_bal", colnames(balance_stats))]
  colnames(balance_stats)[colnames(balance_stats) == stat_var] <- "avg_bal"
  # balance_stats <- balance_stats %>% dplyr::arrange(avg_bal)
  balance_stats <- balance_stats[order(balance_stats$avg_bal), , drop = FALSE]

  x_lab <- if (exposure_type == "continuous") "Correlation with Exposure" else "Standardized Mean Difference Between Exposures"

  labels <- ifelse(balance_stats$balanced == 0, balance_stats$covariate, "")

  min_val <- if (min(balance_stats[, "avg_bal"]) < 0) min(balance_stats[, "avg_bal"]) - 0.05 else min(balance_thresh) - 0.05

  if (min_val > -(max(balance_thresh))){
    min_val <- -(max(balance_thresh)) - 0.05 #to make sure user-supplied balance thresh is on the figure
  }

  max_val <- if (max(balance_stats[, "avg_bal"]) > 0) max(balance_stats[, "avg_bal"]) + 0.05 else max(balance_thresh) + 0.05

  if (max_val < max(balance_thresh)) {
    max_val <- max(balance_thresh) + 0.05 #to make sure user-supplied balance thresh is on the figure
  }

  # Make love plot per exposure time point
  lp <- ggplot2::ggplot(ggplot2::aes(x = avg_bal,
                                     y = reorder(as.factor(covariate), avg_bal)),
                        data = balance_stats) +
    ggplot2::geom_point(ggplot2::aes(y = reorder(as.factor(covariate), avg_bal),
                                     x = avg_bal,
                                     fill = "white",
                                     alpha = 1)) +
    ggplot2::geom_text(ggplot2::aes(label = labels,
                                    hjust = -0.2,
                                    vjust = 0.2),
                       size = 1.5,
                       color = "red") +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab("Covariate") +
    ggplot2::xlim(min_val, max_val) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                   axis.text.x = ggplot2::element_text(color = "black"),
                   axis.text.y = ggplot2::element_text(color = "black"),
                   axis.text = ggplot2::element_text(size = 8),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        color = "black"),
                   plot.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10),
                   legend.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (nrow(balance_stats) > 40) { # stagger covariate labels if there are many
    lp <- lp + ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
  }

  if (!is.null(imp_conf)) { #adding threshold lines
    lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh[1],
                                   linetype = "dashed",
                                   color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh[1],
                                   linetype = "dashed",
                                   color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh[2],
                                   linetype = "dashed",
                                   color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh[2],
                                   linetype = "dashed",
                                   color = "red")
  }
  else {
    lp <- lp + ggplot2::geom_vline(xintercept = balance_thresh,
                                   linetype = "dashed",
                                   color = "red")
    lp <- lp + ggplot2::geom_vline(xintercept = -balance_thresh,
                                   linetype = "dashed",
                                   color = "red")

  }

  if (data_type == "imputed") {
    lp <- lp + ggplot2::ggtitle(paste0(exposure, " (t = ", exposure_time_pt, ") Balance for Imputation ", k))

    if (save.out) {
      suppressMessages(ggplot2::ggsave(lp,
                                       filename = sprintf("%s/balance/%splots/%s_imp_%s_%s_%s_%s_summary_balance_plot.jpeg",
                                                              home_dir, folder, form_name, k, exposure, exposure_time_pt, weights_method),
                                       width = 6,
                                       height = 8))
    }
  }
  else {
    lp <- lp + ggplot2::ggtitle(paste0(exposure, " (t = ", exposure_time_pt, ") Balance"))

    if (save.out) {
      suppressMessages(ggplot2::ggsave(lp,
                                       filename =  sprintf("%s/balance/%splots/%s_%s_%s_%s_summary_balance_plot.jpeg",
                                                               home_dir, folder, form_name, exposure, exposure_time_pt, weights_method),
                                       width = 6,
                                       height = 8))
    }
  }

  if (verbose) {
    print(lp)
  }

}
