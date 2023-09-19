#' Visualize distribution of sample across exposure histories
#'
#' Create customized, user-specified exposure histories and tables displaying
#' sample distribution across them for user inspection.
#'
#' @param data data in wide format as: a data frame, list of imputed
#'   data frames, or mids object
#' @param exposure name of exposure variable
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix
#' @param epochs (optional) data frame of exposure epoch labels and values
#' @param time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param hi_lo_cut list of two numbers indicating quantile values that reflect
#'   high and low values, respectively, for continuous exposure
#' @param ref (optional) string of "-"-separated "l" and "h" values
#'   indicative of a reference exposure history to which to compare comparison,
#'   required if comparison is supplied
#' @param comps (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is supplied
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @return none
#' @export
#' @examples
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
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                time_pts = c(1, 2, 3))
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                time_pts = c(1, 2, 3),
#'                epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
#'                                    values = I(list(c(1, 2), c(3)))))
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                time_pts = c(1, 2, 3),
#'                hi_lo_cut = c(0.6, 0.3))
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                time_pts = c(1, 2, 3),
#'                hi_lo_cut = c(0.6, 0.3),
#'                ref = "l-l-l",
#'                comps = "h-h-h")

eval_hist <- function(data, exposure, tv_confounders, epochs = NULL, time_pts, hi_lo_cut = NULL, ref = NA, comps = NULL, verbose = TRUE){

  exposure_type <- ifelse(inherits(data[, paste0(exposure, '.', time_pts[1])], "numeric"), "continuous", "binary")

  time_varying_wide <- tv_confounders
  time_varying_wide <- sort(time_varying_wide)
  time_varying_wide <- c("ID", time_varying_wide)
  data_wide <- data


  #new will always have exposure main effects (ether exposure time points or epochs)
  # Lists out exposure-epoch combos
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(time_pts),
                         values = time_pts)
    new <- data[, c("ID", paste(exposure, time_pts, sep = "."))]

  } else{

    #new will have cols for epochs
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
        z <- as.numeric(as.character(unlist(data_wide[, which(grepl(paste0(exposure, ".", level), names(data_wide)))])))
        temp <- cbind(temp, z)
      }
      new <- new %>%
        dplyr::mutate(!!new_var := rowMeans(temp, na.rm = TRUE))
    }
  }

  if(verbose){
    cat("Summary of Exposure Main Effects:", "\n")
    print(psych::describe(new %>% select(-c("ID")), fast = TRUE))
    cat("\n")

  }


  tot_hist <- apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                    paste, sep = "", collapse = "-")

  # Assigning history (e.g., h-h-h) based on user-specified hi/lo cutoffs
  if( !is.na(ref) & !is.null(comps)){
    tot_hist <- tot_hist[tot_hist %in% c(ref, comps)]
  }


  epochs$epochs <- as.character(epochs$epochs)


  if (exposure_type == "continuous"){
    if(is.null(hi_lo_cut)){
      # use median as hi/lo split (default)
      new$history <- lapply(seq_len(nrow(new)), function(x) {
        paste(lapply(seq_len(nrow(epochs)), function(y) {
          if (is.na(new[x, y + 1])) {
            return(NA)
          }
          if (new[x, y + 1] >= as.numeric(median((new[, y + 1]), na.rm = TRUE) + 0.001)) { #added 0.001 bc marginaleffects needs different values
            return("h")
          }
          if (new[x, y + 1] <= as.numeric(median((new[, y + 1]), na.rm = TRUE) - 0.001)) {
            return("l")
          }
        }), collapse = "-")
      })
    }
    else{ #user-supplied values
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

      new$history <- lapply(seq_len(nrow(new)), function(x) {
        paste(lapply(seq_len(nrow(epochs)), function(y) {
          if (is.na(new[x, y + 1])) {
            return(NA)
          }
          if (unname(new[x, y + 1]) >= as.numeric(quantile(unname(new[, y + 1]),
                                                           probs = hi_lo_cut[1], na.rm = TRUE))) {
            return("h")
          }
          if (unname(new[x, y + 1]) <= as.numeric(quantile(unname(new[, y + 1]),
                                                           probs =  hi_lo_cut[2], na.rm = TRUE))) {
            return("l")
          }
        }), collapse = "-")
      })
    }
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
    })
  }

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

  if(verbose){
    if(!is.null(hi_lo_cut)){
      cat(paste0("USER ALERT: Out of the total of ", nrow(data_wide), " individuals in the sample, below is the distribution of the ", sum(his_summ$n), " (",
                 round((sum(his_summ$n) / nrow(data_wide)) * 100, 2), "%) that fall into ", nrow(his_summ), " out of the ", length(tot_hist),
                 " the total user-defined exposure histories created from ",
                 hi_lo_cut[2] * 100, "th and ", hi_lo_cut[1] * 100, "th percentile values for low and high levels of exposure ", exposure,
                 ", respectively, across ", paste(epochs$epochs, collapse = ", ")), "\n")
    }
    else{
      cat(paste0("USER ALERT: Out of the total of ", nrow(data_wide), " individuals in the sample, below is the distribution of the ", sum(his_summ$n), " (",
                 round((sum(his_summ$n) / nrow(data_wide)) * 100, 2), "%) that fall into ", nrow(his_summ), " out of the ", length(tot_hist),
                 " the total user-defined exposure histories created from median split values for low and high levels of exposure ", exposure,
                 ", respectively, across ", paste(epochs$epochs, collapse = ", ")), "\n")
    }

    cat("USER ALERT: Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision:", "\n")

    if (nrow(his_summ) != length(tot_hist)) {
      cat(paste0("USER ALERT: There are no individuals in your sample that fall into ",
                 paste(tot_hist[!tot_hist %in% his_summ$history], collapse = " & "),
                 " exposure history/histories. You may wish to consider different high/low cutoffs (for continuous exposures), alternative epochs, or choose a different measure to avoid extrapolation."), "\n")
      cat("\n")
    }

    cat("\n")
    cat(knitr::kable(his_summ, caption = paste0("Summary of user-specified exposure ", exposure,
                                                " histories based on exposure main effects ", paste(epochs$epochs, collapse = ", "),
                                                " containing time points ", paste(epochs$values, collapse = ", "), ":"),
                     format = 'pipe', row.names = F), sep = "\n")
  }

  cat("\n")
}
