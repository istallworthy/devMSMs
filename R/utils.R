# utils ----
.create_dir_if_needed <- function(dir) {
  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir, recurse = TRUE)
  }
}

# Adaptation of gtools::permutations(2, repeats.allowed = TRUE)
perm2 <- function(r, v) {
  if (!is.numeric(r) || length(r) != 1 || r < 1 || (r %% 1) != 0) {
    stop("bad value of r")
  }
  if (!is.atomic(v)) {
    stop("v is non-atomic")
  }

  v <- unique(sort(v))
  n <- length(v)

  sub <- function(r, v) {
    if (r == 1) {
      return(matrix(v, ncol = 1))
    }

    inner <- Recall(r - 1, v)
    cbind(
      rep(v, rep(nrow(inner), n)), 
      matrix(t(inner), ncol = ncol(inner), nrow = nrow(inner) * n, byrow = TRUE)
    )
  }

  sub(r, v)
}

# Base R rewrite of stringr::str_count(); pattern can only be length 1
.string_count <- function(string, pattern = "") {
  lengths(regmatches(string, gregexpr(pattern, string)))
}

# Checks ----
.check_data <- function(data) {
  is_mids <- inherits(data, "mids")
  is_df <- is.data.frame(data)
  is_list_of_df <- is.list(data) && all(vapply(data, is.data.frame, logical(1L)))
  if (!(is_mids || is_df || is_list_of_df)) {
    stop(
      "Please provide either a 'mids' object, a data frame, or a list of imputed datasets in the `data` argument.",
      call. = FALSE
    )
  }
  if (is_mids) {
    rlang::check_installed("mice")
  }
  return(invisible(NULL))
}

.check_weights <- function(weights) {
  dreamerr::check_arg(weights, "class(devMSM_weights) | list", .message = "The 'weighted' mode of this function requires weights be supplied in the form of output from `createWeights()`.")

  for (i in seq_along(weights)) {
    w <- weights[[i]]
    dreamerr::check_value(
      w, "class(weightitMSM)",
      .arg_name = paste0("weights[[", i, "]]"),
      .message = "Please supply a list of weights output from the `createWeights()` function."
    )
  }

  return(invisible(NULL))
}

# Common functions ----
#' Grabs trailing numbers after the last occurance of `sep` in `vars`
#' @noRd
.extract_time_pts_from_vars <- function(vars, sep = "[\\._]") {
  regex_time_pts <- paste0(sep, "([0-9]+)$")
  sapply(vars, function(var) {
    as.numeric(regmatches(var, regexec(regex_time_pts, var))[[1]][2])
  })
}



#' Gets prefix of variable by numbers after last occurance of `sep` and optionally deleting `sep`
#' @noRd
.remove_time_pts_from_vars <- function(vars, sep = "[\\._]", keep_sep = FALSE) {
  regex_time_pts <- paste0("(", sep, ")", "([0-9]+)$")
  replace <- if (keep_sep) "\\1" else ""
  sapply(
    vars,
    function(var) {
      gsub(regex_time_pts, replace, var)
    },
    USE.NAMES = FALSE
  )
}

# initMSM ----
.create_var_tab <- function(exposure, tv_conf, ti_conf = NULL, sep = "[\\._]") {
  regex_time_pts <- paste0(sep, "([0-9]+)$")

  exposure_time_pts <- .extract_time_pts_from_vars(exposure, sep = sep)
  if (anyNA(exposure_time_pts)) {
    stop("Please supply exposure variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }

  tv_conf_time_pts <- .extract_time_pts_from_vars(tv_conf, sep = sep)
  if (anyNA(tv_conf_time_pts)) {
    stop("Please supply tv_conf variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }

  if (!is.null(ti_conf)) {
    ti_conf_time_pts <- .extract_time_pts_from_vars(ti_conf, sep = sep)
    time_pt_detected <- !is.na(ti_conf_time_pts)
    if (any(time_pt_detected)) {
      warning(sprintf(
        "Potentially time-varying covariate was passed to ti_conf:\n  %s.\nThis is likely because it ends with a seperator and a number. If this variable is exogenous or time-invariant, please ignore.",
        paste(ti_conf[time_pt_detected], collapse = ", ")
      ), call. = FALSE)
    }
  }

  data.frame(
    var = c(exposure, tv_conf, ti_conf),
    type = c(rep("exposure", length(exposure)), rep("tv_conf", length(tv_conf)), 
             rep("ti_conf", length(ti_conf))),
    time = c(exposure_time_pts, tv_conf_time_pts, rep(-1, length(ti_conf)))
  )
}

# assessBalance ----
#' Used in `assessBalance` to define history
#' @noRd
.characterize_exposure <- function(mat, exposure_type, hi_lo_cut = NULL) {
  
  if (!is.null(hi_lo_cut)){
    if (length(hi_lo_cut) == 1) hi_lo_cut = c(hi_lo_cut, hi_lo_cut)
    hi <- lapply(mat, function(x) quantile(x, hi_lo_cut[1])) #IS amended to use hi_lo_cut instead of median
    lo <- lapply(mat, function(x) quantile(x, hi_lo_cut[2]))
  }
  
  if (exposure_type == "continuous") {
    if (!is.null(hi_lo_cut)){
      res <- lapply(1:ncol(mat), function(x) {
        ifelse(mat[, x] <= lo[x], "l", 
               ifelse(mat[, x] > hi[x], "h", NA))
      })
      names(res) <- colnames(mat)
    } else{
      
      res <- lapply(mat, function(x) ifelse(x <= median(x), "l", "h"))
    }
  } else {
    res <- lapply(mat, function(x) ifelse(x == 0, "l", "h"))
  }
  do.call(function(...) paste(..., sep = "-"), res)
}

#' Added by IS: Averages balance stats across imputed datasets: list by time
#' 
#' @noRd
.avg_imp_bal_stats_time <- function(bal_stats) {
  average_sub_lists <- function(lists) {
    summed <- Reduce(`+`, lists)/length(lists)
    summed
  }
  
  # TODO: double check this and maybe cleanup?
  avgs <- lapply(
      lapply(seq_along(bal_stats[[1]]), function(z) {
      # lapply(
        # for a single time point, gets all imps
        lapply(seq_along(bal_stats), function(y) {
          bal_stats[[y]][[z]]$std_bal_stats
        }) # ,
        # # TODO: Absolute value or not? 
        # function(q) abs(q[["std_bal_stats"]]))
    }), 
    average_sub_lists)

  test <- bal_stats[[1]]
  for (i in seq_along(test)) {
    test[[i]]$std_bal_stats <- avgs[[i]]
    test[[i]]$balanced <- ifelse(abs(test[[i]]$std_bal_stats) < test[[i]]$bal_thresh, 1, 0) # recalc balanced
  }

  test
}


# fitModel ----
#' Takes row means of eposure variables in the same epoch
#' 
#' @noRd
.generate_epoch_exposures <- function(data, exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  sp <- split(exposure, epoch)
  res <- lapply(sp, function(epoch_cols) {
    rowMeans(data[, epoch_cols, drop = FALSE], na.rm = TRUE)
  })
  res <- do.call("cbind", res)
  colnames(res) <- paste0(prefix, names(sp))
  return(res)
}

#' Return interactions up to order n of variables in `x`
#' @param x character vector of variable names
#' @param n order of interaction
#' @noRd
.get_interaction_terms <- function(x, n = 2) {
  x <- sprintf("`%s`", x)
  x <- attr(terms(reformulate(x)), "term.labels")
  str <- sprintf("~ (%s)^%s", paste(unique(x), collapse = " + "), n)
  
  attr(terms(as.formula(str)), "term.labels")
}

# compareHistories ----
#' Simple function to get the variable names associated with epochs 
#' @noRd
.get_epoch_var_names <- function(exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  paste0(prefix, unique(epoch))
}

#' Grab values of epoch variables to evaluate in `marginaleffects::predictions`
#' 
#' For continuous variables, grab quantiles of the epoch variables based on hi_lo_cut
#' For binary variables, use 0 and 1.
#' @noRd
.get_avg_predictions_variables <- function(data, epoch_vars, exposure_type, hi_lo_cut) {
  
  if (exposure_type == "continuous") {
    res <- lapply(epoch_vars, function(var_name) {
      qs <- quantile(data[[var_name]], probs = hi_lo_cut, na.rm = TRUE)
      if (length(qs) == 1L) qs <- c(qs - 0.0001, qs + 0.0001)
      return(unname(qs))
    })
  } else {
    res <- lapply(epoch_vars, function(var_name) c(0, 1))
  }
  
  names(res) <- epoch_vars
  res
}

#' Create contrast map using comparison and reference
#' @param histories character vector (from preds$term corresponding to average predictions)
#' @param reference character vector of exposure histories
#' @param comparison character vector of exposure histories
#' @noRd
.create_hypotheses_mat <- function(histories, reference, comparison) {
  # Create a matrix with length(histories) rows.
  # Each column sets the value -1 for the reference group row
  # and the value 1 for the comparison group row. 
  # When multiplied by the predicted values it does 
  # R' b = b_comparison - b_reference
  sub_mats <- lapply(reference, function(ref) {
    n_cols <- length(comparison)
    cus_comps <- matrix(0, ncol = n_cols, nrow = length(histories))
    cus_comps[match(ref, histories), ] <- -1
    # cus_comps[match(comparison, histories), 1:n_cols] <- 1
    cus_comps[cbind(match(comparison, histories), seq_len(n_cols))] <- 1 # changed by IS to appropriately assign 1 to comparison histories
    colnames(cus_comps) <- sprintf("(%s) - (%s)", comparison, ref)
    return(cus_comps)
  })
  
  do.call("cbind", sub_mats)
}

#' Add history labels and dose to table based on `epoch_vars`
#'
#' @param p table output from marginaleffects::avg_predictions() or hypotheses()
#' @param epoch_vars variable names of exposure epoch (in time order)
#' @return table with histories labeled
#' @noRd
.add_histories <- function(p, epoch_vars) {
  # NOTE:
  # 1. This function is based on the observation that there is a low and a high exposure value, so you can classify "h"/"l" by just comparing to the min of a column from `p`
  # 2. For the resulting string to be in the correct order, epoch_vars must be in the correct order. It should be since it is created in `initMSM`.

  # sapply is necessary here because p is wrapped in classes that prevent simple subsetting
  e <- sapply(epoch_vars, function(x) p[[x]])
  classification <- apply(e, 2, function(x) ifelse(x > min(x), "h", "l"))

  # history string, e.g. "l-l-h"
  histories <- apply(classification, 1, function(row) paste(row, collapse = "-"))
  p[["term"]] <- histories
  return(p)
}

#' Calculate number of doses from history string (e.g. "l-l-l")
#' 
#' @examples
#' library(devMSMs)
#' devMSMs:::.dose_from_history(c("l-l-l-h", "l-l-l-l"), dose_level = "h")
#' devMSMs:::.dose_from_history(c("l-l-l-h", "l-l-l-l"), dose_level = "l")
#' @noRd
.dose_from_history <- function(history, dose_level = "h") {
  .string_count(history, dose_level)
}

#' Calculate number of doses from history string (e.g. "l-l-l")
#' 
#' @examples
#' library(devMSMs)
#' devMSMs:::.dose_from_history_contrast(c("(h-h-h) - (l-l-l)"), dose_level = "h")
#' devMSMs:::.dose_from_history_contrast(c("(h-h-h) - (l-l-l)"), dose_level = "l")
#' @noRd
.dose_from_history_contrast <- function(history, dose_level = "h") {
  vapply(
    strsplit(history, " - "),
    function(hist) {
      counts <- .string_count(hist, dose_level)
      sprintf("%s - %s", counts[1], counts[2])
    },
    character(1L)
  )
}



