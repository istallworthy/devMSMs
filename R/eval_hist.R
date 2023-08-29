#' Visualize distribution of sample across exposure histories
#'
#' Create customized, user-specified exposure histories and tables displaying
#' sample distribution across them for user inspection.
#'
#' @param data
#' @param exposure name of exposure variable
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix
#' @param epochs
#' @param time_pts
#' @param hi_lo_cut list of two numbers indicating quantile values that reflect
#'   high and low values, respectively, for continuous exposure
#' @param ref
#' @param comps
#'
#' @return none
#' @export
#'
#' @examples
eval_hist <- function(data, exposure, tv_confounders, epochs, time_pts, hi_lo_cut, ref, comps){

  exposure_type <- ifelse(inherits(data[, paste0(exposure, '.', time_pts[1])], "numeric"), "continuous", "binary")

  epochs$epochs <- as.character(epochs$epochs)
  time_varying_wide <- tv_confounders
  time_varying_wide <- sort(time_varying_wide)
  time_varying_wide <- c("ID", time_varying_wide)
  data_wide <- data
  new <- data.frame(ID = data_wide[, "ID"])
  colnames(new) <- "ID"

  # Averages exposure across time points that constitute the exposure epochs (e.g., infancy = 6 & 15)
  for (e in seq_len(nrow(epochs))) {
    epoch <- epochs[e, 1]
    temp <- data.frame(row.names = seq_len(nrow(data_wide)))
    new_var <- paste0(exposure, "_", epoch)
    # Finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
    for (l in seq_len(length(as.numeric(unlist(epochs[e, 2]))))) {
      level <- as.numeric(unlist(epochs[e, 2]))[l]
      z <- data_wide[, which(grepl(paste0(exposure, ".", level), names(data_wide)))]
      temp <- cbind(temp, z)
    }
    new <- new %>%
      dplyr::mutate(!!new_var := rowMeans(temp, na.rm = TRUE))
  }

  # Assigning history (e.g., h-h-h) based on user-specified hi/lo cutoffs
  tot_hist <- apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                    paste, sep = "", collapse = "-")

  if( !is.na(ref) & !is.null(comps)){
    tot_hist <- tot_hist[tot_hist %in% c(ref, comps)]
  }


  if (exposure_type == "continuous"){
    if(is.null(hi_lo_cut)){
      epoch_info$low[t] <- as.numeric(median(data[, var_name], na.rm = T))
      epoch_info$high[t] <- as.numeric(median(data[, var_name], na.rm = T))

      new$history <- lapply(seq_len(nrow(new)), function(x) {
        paste(lapply(seq_len(nrow(epochs)), function(y) {
          if (is.na(new[x, y + 1])) {
            return(NA)
          }
          if (new[x, y + 1] >= as.numeric(median(new[, y + 1], na.rm = TRUE))) {
            return("h")
          }
          if (new[x, y + 1] <= as.numeric(median(new[, y + 1], na.rm = TRUE))) {
            return("l")
          }
        }), collapse = "-")
      })

    }
    else{
      if (length(hi_lo_cut) == 2){
        hi_cutoff <- hi_lo_cut[1]
        lo_cutoff <- hi_lo_cut[2]

        if (hi_cutoff > 1 || hi_cutoff < 0) {
          stop('Please select a high cutoff value between 0 and 1', call. = FALSE)
        }
        if (lo_cutoff > 1 || lo_cutoff < 0) {
          stop('Please select low cutoff value between 0 and 1', call. = FALSE)
        }
      }
      else{
        if (hi_lo_cut > 1 || hi_lo_cut < 0) {
          stop('Please select a hi_lo cutoff value between 0 and 1', call. = FALSE)
        }
        hi_cutoff <- hi_lo_cut
        lo_cutoff <- hi_lo_cut
      }
    }

    new$history <- lapply(seq_len(nrow(new)), function(x) {
      paste(lapply(seq_len(nrow(epochs)), function(y) {
        if (is.na(new[x, y + 1])) {
          return(NA)
        }
        if (new[x, y + 1] >= as.numeric(quantile(new[, y + 1], probs = hi_lo_cut[1], na.rm = TRUE))) {
          return("h")
        }
        if (new[x, y + 1] <= as.numeric(quantile(new[, y + 1], probs =  hi_lo_cut[2], na.rm = TRUE))) {
          return("l")
        }
      }), collapse = "-")
    })

  }

  else if (exposure_type == "binary"){
    new$history <- lapply(seq_len(nrow(new)), function(x) {
      paste(lapply(seq_len(nrow(epochs)), function(y) {
        if (is.na(new[x, y + 1])) {
          return(NA)
        }
        if (new[x, y + 1] == 1) {
          return("h")
        }
        if (new[x, y + 1] == 0) {
          return("l")
        }
      }), collapse = "-")
    })}

  # Summarizing n's by history
  his_summ <- new %>%
    dplyr::group_by(history) %>%
    dplyr::summarize(n = dplyr::n())

  if( !is.na(ref) & !is.null(comps)){
    his_summ <- his_summ %>%
      dplyr::filter(history %in% c(ref, comps))
  }

  his_summ <- his_summ[! grepl("NA", his_summ$history),]
  his_summ <- his_summ[! grepl("NULL", his_summ$history),]

  cat(paste0("USER ALERT: Out of the total of ", nrow(data_wide), " individuals in the sample, below is the distribution of the ", sum(his_summ$n), " (",
             round((sum(his_summ$n) / nrow(data_wide)) * 100, 2), "%) that fall into ", nrow(his_summ), " out of the ", length(tot_hist),
             " the total user-defined exposure histories created from ",
             hi_lo_cut[2] * 100, "th and ", hi_lo_cut[1] * 100, "th percentile values for low and high levels of exposure ", exposure,
             ", respectively, across ", paste(epochs$epochs, collapse = ", ")), "\n")

  cat("Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision:", "\n")

  if (nrow(his_summ) != length(tot_hist)) {
    cat(paste0("USER ALERT: There are no individuals in your sample that fall into ",
               paste(tot_hist[!tot_hist %in% his_summ$history], collapse = " & "),
               " exposure history/histories. You may wish to consider different high/low cutoffs (for continuous exposures), alternative epochs, or choose a different measure to avoid extrapolation."), "\n")
    cat("\n")
  }

  cat("\n")

  cat(knitr::kable(his_summ, caption = paste0("Summary of User-Specified Exposure ", exposure,
                                              " Histories Based on Exposure Epochs ", paste(epochs$epochs, collapse = ", "),
                                              " containing time points ", paste(epochs$values, collapse = ", "),
                                              " and the following high/low cutoffs: ", paste(hi_lo_cut, collapse = " & ")),
                   format = 'pipe', row.names = F), sep = "\n")

  cat("\n")
}
