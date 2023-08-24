# FUNCTIONS called by compareHistories

#custom contrasts
get_reference_values <- function(d, reference) {
  ref_vals <- sapply(1:length(unlist(strsplit(reference, "-"))), function(x) {
    d[x,unlist(strsplit(reference, "-"))[x]]
  })
  ref_vals
}


get_comparison_values <- function(d, comp_histories) {
  comp_vals <- sapply(comp_histories, function(comp) {
    sapply(1:length(unlist(strsplit(comp, "-"))), function(x) {
      d[x, unlist(strsplit(comp, "-"))[x]]
    })
  })
  t(comp_vals)
}


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


create_custom_comparisons <- function(preds, ref_vals, comp_vals, exposure) {
  cus_comps <- matrix(ncol = nrow(comp_vals), nrow = nrow(as.data.frame(preds[[1]][[1]])))

  ref_pos <- which(apply(as.data.frame(preds[[1]])[grepl(exposure, colnames(as.data.frame(preds[[1]])))], 1,
                         paste, collapse = ",") == paste0(ref_vals, collapse = ","))
  cus_comps[ref_pos, ] <- -1

  for (x in 1:nrow(comp_vals)) {
    c <- paste(comp_vals[x, ], collapse = ",")
    cus_comps[which(apply(as.data.frame(preds[[1]])[grepl(exposure, colnames(as.data.frame(preds[[1]])))], 1,
                          paste, collapse = ",") == paste0(c, collapse = ",")), x] <- 1
  }

  if (nrow(comp_vals) > 1) {
    colnames(cus_comps) <- paste0("(", paste0(paste0(ref_vals, collapse = ","), ") - (",
                                              apply(comp_vals, 1, paste, collapse = ",")), ")")
  } else {
    colnames(cus_comps) <- paste0("(", paste0(paste0(ref_vals, collapse = ","), ") - (",
                                              paste0(paste0(comp_vals, collapse = ",")), ")"))
  }

  cus_comps[is.na(cus_comps)] <- 0
  cus_comps
}



add_histories <- function(p, d) {
  if( ("list" %in% class(p)) & length(p) == 1){
    history <- matrix(data = NA, nrow = nrow(p[[1]]), ncol = 1) # Get histories from the first element
    p <- p[[1]]
  }

  if (sum(d$e %in% colnames(p)) == nrow(d)){
    for (i in 1:nrow(p)) {
      vals <- as.data.frame(p)[i, 1:nrow(d)]
      history[i] <- as.data.frame(paste(ifelse(round(vals, 3) == round(d$l, 3), "l", "h"), collapse = "-"))
    }
    history <- unlist(history)
  }

  if("term" %in% colnames(p)) { #preds_pool
    history <- matrix(data = NA, nrow = nrow(p), ncol = 1) # Get histories from the first element

    if(grepl("\\=", p$term[1])) {
      for (i in 1:nrow(p)){
        vals <- as.numeric(sapply(strsplit(unlist(strsplit(as.character(p$term[i]), "\\,")), "="), "[",2))
        history[i] <- as.data.frame(paste(ifelse(round(vals, digits = 2) == round(d$l, digits = 2), "l", "h"), collapse = "-"))
      }
      history <- unlist(history)

    } else {
      for (i in 1:nrow(p)) {
        temp <- as.character(p$term[i])
        pair <- lapply(1:2, function(y) {
          a <- sapply(strsplit(temp, " - "), "[", y)
          his <- lapply(1:nrow(d), function(z) {
            ifelse(round(as.numeric(gsub("[^0-9.-]", "",
                                         sapply(strsplit(a, "\\,"), "[", z))), 2) == round(d[z, "l"], 2), "l", "h")
          })
        })
        history[i, 1] <- paste(sapply(pair, paste, collapse = "-"), collapse = " vs ")
      }
    }
  }
  p <- cbind (p, history = history)
  p
}


add_dose <- function(p, dose_level) {
  if( length(p$history[1]) == 1 ) {
    if(grepl("vs", p$history[1])) {
      dose_a <- stringr::str_count(sapply(strsplit(p$history, "vs"), "[", 1), dose_level)
      dose_b <- stringr::str_count(sapply(strsplit(p$history, "vs"), "[", 2), dose_level)
      dose_count <- data.frame(dose = gsub(" ", " vs ", paste(dose_a, dose_b)))
    } else{
      dose_count <- stringr::str_count(p$history, dose_level)
    }
  }
  if (length(p$history[1]) > 1){
    dose_count <- stringr::str_count(p$history, dose_level)
  }
  p <- cbind (p, dose_count = dose_count)
  p
}


perform_multiple_comparison_correction <- function(comps, reference, comp_histories, method) {
  if (any(is.na(reference) & is.na(comp_histories)) | length(comp_histories) > 1) {
    cat("\n")
    cat(paste0("Conducting multiple comparison correction using the ", method, " method."), "\n")
    cat("\n")
    corr_p <-stats::p.adjust(comps$p.value, method = method)
    comps <- cbind(comps, p.value_corr = corr_p)
  } else {
    cat("\n")
    cat(paste0("The user specified comparison only between ", reference, " and a single comparison, ", comp_histories,
               ", so no correction for multiple comparisons will be implemented."), "\n")
  }
  comps
}
