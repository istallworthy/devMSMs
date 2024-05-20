#' Visualize distribution of sample across exposure histories
#'
#' Create customized, user-specified exposure histories and tables displaying
#' sample distribution across them for user inspection.
#'
#' @param epoch_history character vector of length `nrow(data)` containing epochi histories for 
#'  each unit 
#' @param epoch character vector of epoch values
#' @inheritParams compareHistories 
#' 
#' @return invisibly NULL
#' 
#' @keywords internal
print_eval_hist <- function(epoch_history, epoch, hi_lo_cut, reference = NULL, comparison = NULL) {

  if (length(hi_lo_cut) == 1) hi_lo_cut = c(hi_lo_cut, hi_lo_cut)

  epoch_history_tab <- data.frame(table(epoch_history, useNA = "ifany"))
  colnames(epoch_history_tab) <- c("epoch_history", "n")
  n_total <- sum(epoch_history_tab$n)
  n_total_hist <- nrow(epoch_history_tab)
  
  if (!is.null(reference) && !is.null(comparison)) {
    keep_idx <- epoch_history_tab$epoch_history %in% c(reference, comparison)
    epoch_history_tab <- epoch_history_tab[keep_idx, , drop = FALSE]
  }
  n_included <- sum(epoch_history_tab$n)
  n_included_hist <- nrow(epoch_history_tab)

  message(sprintf(
    "USER ALERT: Out of the total of %s individuals in the sample, below is the distribution of the %s (%.0f%%) individuals that fall into %s user-selected exposure histories (out of the %s total) created from %sth and %sth percentile values for low and high levels of exposure-epoch %s. \n",
    n_total,
    n_included,
    n_included / n_total * 100,
    n_included_hist,
    n_total_hist,
    hi_lo_cut[2] * 100,
    hi_lo_cut[1] * 100,
    paste(unique(epoch), collapse = ", ")
  ))

  message("USER ALERT: Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision:\n")

  tab_caption <- sprintf(
    "Summary of user-selected exposure histories based on exposure main effects %s:",
    paste(unique(epoch), collapse = ", ")
  )
  tab <- tinytable::tt(epoch_history_tab, caption = tab_caption)  
  print(tab, "markdown")
  cat("\n")
  
  return(invisible(NULL))
}
