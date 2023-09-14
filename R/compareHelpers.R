# FUNCTIONS called by compareHistories


#' Finds custom reference values
#'
#' @param d data frame of high and low values per exposure main effect
#' @param reference reference sequence of "h" and/or "l" (e.g., "h-h-h")
#' @return reference values
#' @examples
#' d <- data.frame(e = c("A.1", "A.2", "A.3"),
#'                l = c(0, 0, 0),
#'                h = c(1, 1, 1))
#' r <- get_reference_values(d = d,
#'                           reference = "l-l-l" )
#' r <- get_reference_values(d = d,
#'                           reference = "h-h-h" )

get_reference_values <- function(d, reference) {
  ref_vals <- sapply(seq_len(length(unlist(strsplit(reference, "-")))), function(x) {
    d[x, unlist(strsplit(reference, "-"))[x]]
  })
  ref_vals
}


#' Finds custom comparison values
#'
#' @param d data frame of high and low values per exposure main effect
#' @param comp_histories comparison sequence(s) of "h" and/or "l" (e.g., "h-h-h")
#' @return comparison values
#' @export
#' @examples
#' d <- data.frame(e = c("A.1", "A.2", "A.3"),
#'                l = c(0, 0, 0),
#'                h = c(1, 1, 1))
#' r <- get_comparison_values(d = d,
#'                           comp_histories = "l-l-l" )
#' r <- get_comparison_values(d = d,
#'                           comp_histories = "h-h-h" )
#' r <- get_comparison_values(d = d,
#'                          comp_histories = c("h-h-h", "h-h-l"))
get_comparison_values <- function(d, comp_histories) {
  comp_vals <- sapply(comp_histories, function(comp) {
    sapply(seq_len(length(unlist(strsplit(comp, "-")))), function(x) {
      d[x, unlist(strsplit(comp, "-"))[x]]
    })
  })
  t(comp_vals)
}


#' Create custom contrasts
#'
#' @param d data frame of high and low values per exposure main effect
#' @param reference reference sequence of "h" and/or "l" (e.g., "h-h-h")
#' @param comp_histories comparison sequence(s) of "h" and/or "l" (e.g., "h-h-h")
#' @param exposure name of exposure variable
#' @param preds custom output of marginaleffects::average_predictions()
#'
#' @return contrasts
#' @export
create_custom_contrasts <- function(d, reference, comp_histories, exposure, preds) {
  if (is.na(reference) | is.null(comp_histories)) {
    return(NULL)  # Invalid input, return early
  }

  ref_vals <- get_reference_values(d, reference)
  comp_vals <- get_comparison_values(d, comp_histories)
  cus_comps <- create_custom_comparisons(preds, ref_vals, comp_vals, exposure)

  comps <- lapply(preds, function(y) {
    y |> marginaleffects::hypotheses(cus_comps)
  })

  comps
}


#' Creates custom comparisons
#'
#' @param preds custom output of marginaleffects::average_predictions()
#' @param ref_vals reference values
#' @param comp_vals comparison values
#' @param exposure name of exposure variable
#' @return custom comparisons
#' @export

create_custom_comparisons <- function(preds, ref_vals, comp_vals, exposure) {
  cus_comps <- matrix(ncol = nrow(comp_vals), nrow = nrow(as.data.frame(preds[[1]][[1]])))

  ref_pos <- which(apply(as.data.frame(preds[[1]])[grepl(exposure, colnames(as.data.frame(preds[[1]])))], 1,
                         paste, collapse = ",") == paste0(ref_vals, collapse = ","))
  cus_comps[ref_pos, ] <- -1

  for (x in seq_len(nrow(comp_vals))) {
    c <- paste(comp_vals[x, ], collapse = ",")
    cus_comps[which(apply(as.data.frame(preds[[1]])[grepl(exposure, colnames(as.data.frame(preds[[1]])))], 1,
                          paste, collapse = ",") == paste0(c, collapse = ",")), x] <- 1
  }

  if (nrow(comp_vals) > 1) {
    colnames(cus_comps) <- paste0("(", paste0(paste0(ref_vals, collapse = ","), ") - (",
                                              apply(comp_vals, 1, paste, collapse = ",")), ")")
  }
  else {
    colnames(cus_comps) <- paste0("(", paste0(paste0(ref_vals, collapse = ","), ") - (",
                                              paste0(paste0(comp_vals, collapse = ",")), ")"))
  }

  cus_comps[is.na(cus_comps)] <- 0
  cus_comps
}



#' Add history labels to table
#'
#' @param p table output from marginaleffects::avg_predictions() or hypotheses()
#' @param d data frame of high and low values per exposure main effect
#'
#' @return table with histories labeled
#' @export

add_histories <- function(p, d) {
  if((is.list(p)) & length(p) == 1){
    history <- matrix(data = NA, nrow = nrow(p[[1]]), ncol = 1) # Get histories from the first element
    p <- p[[1]]
  }

  #for preds
  if (sum(d$e %in% colnames(p)) == nrow(d)){
    for (i in seq_len(nrow(p))) {
      vals <- as.data.frame(p)[i, seq_len(nrow(d))]
      history[i] <- as.data.frame(paste(ifelse(round(as.numeric(as.character(vals)), 3) ==
                                                 round(as.numeric(as.character(d$l)), 3), "l", "h"), collapse = "-")) #was getting different rounding values due to floaters (i think), so made char first
    }
    history <- unlist(history)
  }

  if("term" %in% colnames(p)) { #preds_pool, comps
    history <- matrix(data = NA, nrow = nrow(p), ncol = 1) # Get histories from the first element

    if(grepl("\\=", p$term[1])) {
      for (i in seq_len(nrow(p))){
        vals <- as.numeric(sapply(strsplit(unlist(strsplit(as.character(p$term[i]), "\\,")), "="), "[", 2))
        history[i] <- as.data.frame(paste(ifelse(round(as.numeric(as.character(vals)), digits = 3) ==
                                                   round(as.numeric(as.character(d$l)), digits = 3), "l", "h"), collapse = "-"))
      }
      history <- unlist(history)

    }
    else { #comps
      for (i in seq_len(nrow(p))) {
        temp <- as.character(p$term[i])
        pair <- lapply(1:2, function(y) {
          a <- sapply(strsplit(temp, " - "), "[", y)
          his <- lapply(seq_len(nrow(d)), function(z) {
            ifelse(round(as.numeric(as.character(gsub("[^0-9.-]", "", sapply(strsplit(a, "\\,"), "[", z)))), 3) ==
                     round(as.numeric(as.character(d[z, "l"])), 3), "l", "h")
          })
        })
        history[i, 1] <- paste(sapply(pair, paste, collapse = "-"), collapse = " vs ")
      }
    }
  }
  p <- cbind (p, history = history)
  p
}


#' Add dose tally to table
#'
#' @param p table output from marginaleffects::avg_predictions() or hypotheses()
#' @param dose_level "l" or "h" indicating whether low or high doses should be tallied in tables and plots
#' @return table with dose level tally
#' @export

add_dose <- function(p, dose_level) {
  if( length(p$history[1]) == 1 ) {
    if(grepl("vs", p$history[1])) {
      dose_a <- stringr::str_count(sapply(strsplit(p$history, "vs"), "[", 1), dose_level)
      dose_b <- stringr::str_count(sapply(strsplit(p$history, "vs"), "[", 2), dose_level)
      dose_count <- data.frame(dose = gsub(" ", " vs ", paste(dose_a, dose_b)))
    }
    else{
      dose_count <- stringr::str_count(p$history, dose_level)
    }
  }
  if (length(p$history[1]) > 1){
    dose_count <- stringr::str_count(p$history, dose_level)
  }
  p <- cbind (p, dose_count = dose_count)
  p
}


#' Conduct multiple comparison correction
#'
#' @param comps
#' @param reference reference sequence of "h" and/or "l" (e.g., "h-h-h")
#' @param comp_histories comparison sequence(s) of "h" and/or "l" (e.g., "h-h-h")
#' @param method character abbreviation for multiple comparison correction method
#' @return comparison table with corrected p-values
#' @export

perform_multiple_comparison_correction <- function(comps, reference, comp_histories, method) {
  if (nrow(comps) > 1) {
    cat("\n")
    cat(paste0("Conducting multiple comparison correction using the ", method, " method."), "\n")
    cat("\n")
    corr_p <-stats::p.adjust(comps$p.value, method = method)
    comps <- cbind(comps, p.value_corr = corr_p)
  }
  else {
    cat("\n")
    cat(paste0("The user specified comparison only between ", reference, " and a single comparison, ", comp_histories,
               ", so no correction for multiple comparisons will be implemented."), "\n")
  }
  comps
}
