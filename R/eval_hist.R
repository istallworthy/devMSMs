#' Visualize distribution of sample across exposure histories
#'
#' Create customized, user-specified exposure histories and tables displaying
#' sample distribution across them for user inspection.
#'
#' @param data data in wide format as: a data frame, list of imputed
#'   data frames, or mids object
#' @param exposure name of exposure variable
#' @param epochs (optional) data frame of exposure epoch labels and values
#' @param time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param hi_lo_cut list of two numbers indicating quantile values that reflect
#'   high and low values, respectively, for continuous exposure
#' @param ref (optional) list of one or more strings of "-"-separated "l" and "h" values
#'   indicative of a reference exposure history to which to compare comparison,
#'   required if comparison is supplied
#' @param comps (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is supplied
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
#'                time_pts = c(1, 2, 3))
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                time_pts = c(1, 2, 3),
#'                epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
#'                                    values = I(list(c(1, 2), c(3)))))
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                time_pts = c(1, 2, 3),
#'                hi_lo_cut = c(0.6, 0.3))
#' h <- eval_hist(data = test,
#'                exposure = "A",
#'                time_pts = c(1, 2, 3),
#'                hi_lo_cut = c(0.6, 0.3),
#'                ref = "l-l-l",
#'                comps = "h-h-h")

eval_hist <- function(data, exposure, time_pts, epochs = NULL, hi_lo_cut = NULL, 
                      ref = NULL, comps = NULL) {
  
  rlang::check_installed("knitr")
  
  
  exposure_type <- if (inherits(data[, paste0(exposure, '.', time_pts[1])], "numeric")) "continuous" else "binary"
  
  data_wide <- data
  
  #new will always have exposure main effects (ether exposure time points or epochs)
  # Lists out exposure-epoch combos
  
  if (is.null(epochs)) { #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(time_pts),
                         values = time_pts)
    new <- data[, c("ID", paste(exposure, time_pts, sep = "."))]
    
  }
  else {
    
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
        z <- as.data.frame(data_wide[, names(data_wide)[grepl(exposure, 
                                                              names(data_wide))]]) #finds exposure vars
        cols <- colnames(z)[as.logical(sapply(strsplit(names(z), "\\."), "[", 2) == as.character(level))]
        cols <- cols[!is.na(cols)]
        z <- as.numeric(as.character(unlist(z[, cols])))
        temp <- cbind(temp, z)
      }
      x <- as.data.frame(rowMeans(temp, na.rm = TRUE))
      colnames(x) <- c(new_var)
      new <- cbind(new, x)
    }
  }
  
  cat("Summary of Exposure Main Effects:", "\n")
  summary(new[, !colnames(new) %in% "ID"])
  cat("\n")
  
  
  tot_hist <- apply(perm2(nrow(epochs), c("l", "h")), 1,
                    paste, sep = "", collapse = "-")
  
  # Assigning history (e.g., h-h-h) based on user-specified hi/lo cutoffs
  
  if (!is.null(ref) && !is.null(comps)) {
    tot_hist <- tot_hist[tot_hist %in% c(ref, comps)]
  }
  
  epochs$epochs <- as.character(epochs$epochs)
  
  if (exposure_type == "continuous") {
    
    if (is.null(hi_lo_cut)) {
      
      # use median as hi/lo split (default)
      
      new$history <- lapply(seq_len(nrow(new)), function(x) {
        
        paste(lapply(seq_len(nrow(epochs)), function(y) {
          
          if (is.na(new[x, y + 1])) {
            return(NA)
          }
          if (new[x, y + 1] >= as.numeric(median((new[, y + 1]),
                                                 na.rm = TRUE) + 0.001)) { #added 0.001 bc marginaleffects needs different values
            return("h")
          }
          if (new[x, y + 1] <= as.numeric(median((new[, y + 1]),
                                                 na.rm = TRUE) - 0.001)) {
            return("l")
          }
        }), collapse = "-")
      })
      
    }
    else { #user-supplied values
      if (length(hi_lo_cut) == 2) {
        
        hi_cutoff <- hi_lo_cut[1]
        lo_cutoff <- hi_lo_cut[2]
        
        if (hi_cutoff > 1 || hi_cutoff < 0) {
          stop('Please select a high cutoff value between 0 and 1',
               call. = FALSE)
        }
        if (lo_cutoff > 1 || lo_cutoff < 0) {
          stop('Please select low cutoff value between 0 and 1',
               call. = FALSE)
        }
      }
      else {
        if (hi_lo_cut > 1 || hi_lo_cut < 0) {
          stop('Please select a hi_lo cutoff value between 0 and 1',
               call. = FALSE)
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
                                                           probs = hi_lo_cut[1],
                                                           na.rm = TRUE))) {
            return("h")
          }
          if (unname(new[x, y + 1]) <= as.numeric(quantile(unname(new[, y + 1]),
                                                           probs =  hi_lo_cut[2],
                                                           na.rm = TRUE))) {
            return("l")
          }
        }), collapse = "-")
      })
    }
  }
  
  else if (exposure_type == "binary") {
    
    new$history <- lapply(seq_len(nrow(new)), function (x) {
      
      paste(lapply(seq_len(nrow(epochs)), function (y) {
        
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
  
  new$history <- unlist(new$history)
  his_summ <- aggregate( ID ~ history,
                         data = new,
                         FUN = length)
  colnames(his_summ) <- c("history", "n")
  
  if (!is.null(ref) && !is.null(comps)) {
    his_summ <- his_summ[his_summ$history %in% c(ref, comps), ]
  }
  
  his_summ <- his_summ[!grepl("NA", his_summ$history), ]
  his_summ <- his_summ[!grepl("NULL", his_summ$history), ]
  
  if (!is.null(hi_lo_cut)) {
    
    cat(sprintf("USER ALERT: Out of the total of %s individuals in the sample, below is the distribution of the %s (%s%%) individuals that fall into %s out of the %s the total user-defined exposure histories created from %sth and %sth percentile values for low and high levels of exposure %s, respectively, across %s. \n",
                nrow(data_wide),
                sum(his_summ$n),
                round((sum(his_summ$n) / nrow(data_wide)) * 100, 2),
                nrow(his_summ),
                length(tot_hist),
                hi_lo_cut[2] * 100,
                hi_lo_cut[1] * 100,
                exposure,
                paste(epochs$epochs, collapse = ", ")))
  }
  else {
    
    cat(sprintf("USER ALERT: Out of the total of %s individuals in the sample, below is the distribution of the %s (%s%%) individuals that fall into %s out of the %s total user-defined exposure histories created from median split values for low and high levels of exposure %s, respectively, across %s. \n",
                nrow(data_wide),
                sum(his_summ$n),
                round((sum(his_summ$n) / nrow(data_wide)) * 100, 2),
                nrow(his_summ),
                length(tot_hist),
                exposure,
                paste(epochs$epochs, collapse = ", ")))
  }
  
  cat("USER ALERT: Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision:", "\n")
  
  if (nrow(his_summ) != length(tot_hist)) {
    
    warning(sprintf("USER ALERT: There are no individuals in your sample that fall into %s exposure history/histories. You may wish to consider different high/low cutoffs (for continuous exposures), alternative epochs, or choose a different measure to avoid extrapolation.\n",
                    paste(tot_hist[!tot_hist %in% his_summ$history], 
                          collapse = " & ")), 
            call. = FALSE)
    cat("\n")
  }
  
  cat("\n")
  cat(knitr::kable(his_summ,
                   caption = sprintf("Summary of user-specified exposure %s histories based on exposure main effects %s containing time points %s:",
                                     exposure, paste(epochs$epochs, 
                                                     collapse = ", "), 
                                     paste(epochs$values, collapse = ", ")),
                   format = 'pipe',
                   row.names = F),
      sep = "\n")
  cat("\n")
  
}
