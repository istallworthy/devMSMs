#' Plot results from history comparisons
#' Code to plot predicted values of the different exposure histories by history and colored by dosage of exposure
#' @param msm_object msm object that contains all relevant user inputs
#' @param history_estimates pooled predicted values for the marginal model
#' @return plots
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_errorbarh
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @examples plotResults(msm_object, history_estimates)
plotResults <- function(msm_object, history_estimates) {
  home_dir <- msm_object$home_dir
  exposure <- msm_object$exposure
  outcome <- msm_object$outcome
  exposure_labels <- msm_object$exposure_labels
  outcome_labels <- msm_object$outcome_labels
  colors <- msm_object$colors
  weights_percentile_cutoff <- msm_object$weights_percentile_cutoff
  dose_level <- msm_object$dose_level

  # Error checking
  if (length(colors) > 1 & length(colors) != nrow(msm_object$exposure_epochs) + 1) {
    stop(paste0('Please provide either: ', nrow(msm_object$exposure_epochs) + 1,
                ' different colors, a color palette, or leave this entry blank in the msm object'))
  }

  for (x in seq_along(history_estimates)) {
    cutoff <- names(history_estimates)[x]
    comparisons <- data.frame(history_estimates[[x]])
    comparisons$term <- gsub(paste0(exposure, "_"), "", comparisons$term)
    comparisons$low_ci <- comparisons$estimate - (1.96 * comparisons$std.error)
    comparisons$high_ci <- comparisons$estimate + (1.96 * comparisons$std.error)
    comparisons$history <- as.factor(comparisons$history)
    comparisons$dose <- as.factor(comparisons$dose)

    comparisons <- comparisons %>%
      arrange(dose) # Order by dose

    if (length(colors) > 1) { # If user input a list of colors
      p <- ggplot(data = comparisons, aes(x = estimate, y = history, color = dose)) +
        geom_point(size = 5) +
        scale_color_manual(values = colors) +
        scale_y_discrete(limits = c(as.character(comparisons$history)), expand = c(0, 0.2)) +
        geom_errorbarh(aes(xmin = low_ci, xmax = high_ci), height = 0.6) +
        xlab(paste0("Predicted ", outcome_labels, " Value")) +
        ylab(paste0(exposure_labels, " Exposure History")) +
        xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci), max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
        theme(text = element_text(size = 14)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

      ggsave(paste0(home_dir, "results figures/", exposure, "-", outcome, " cutoff= ", cutoff, ".jpeg"), plot = p)
    } else { # If user lists a palette
      p <- ggplot(data = comparisons, aes(x = estimate, y = history, color = dose)) +
        geom_point(size = 5) +
        scale_color_palette(name = "Dosage") +
        scale_y_discrete(limits = c(as.character(comparisons$history)), expand = c(0, 0.2)) +
        geom_errorbarh(aes(xmin = low_ci, xmax = high_ci), height = 0.6) +
        xlab(paste0("Predicted ", outcome_labels, " Value")) +
        ylab(paste0(exposure_labels, " Exposure History")) +
        xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci), max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
        theme(text = element_text(size = 14)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

      ggsave(paste0(home_dir, "results figures/", exposure, "-", outcome, " cutoff= ", cutoff, ".jpeg"), plot = p)
    }

    message("\n")
    message("See the 'results figures/' folder for graphical representations of results.")
  }
}
