# utils ----
.create_dir_if_needed <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

#' Adaptation of gtools::permutations(2, repeats.allowed = TRUE)
#' @keywords internal
perm2 <- function(r, v) {
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) != 0) {
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

#' Base R rewrite of stringr::str_count(); pattern can only be length 1
#' @keywords internal
.string_count <- function(string, pattern = "") {
  lengths(regmatches(string, gregexpr(pattern, string)))
}

# Checks ----
.check_data <- function(data) {
  is_mids <- inherits(data, "mids")
  is_df <- is.data.frame(data)
  is_list_of_df <- is.list(data) && all(sapply(data, is.data.frame))
  if (!(is_mids | is_df | is_list_of_df)) {
    stop(
      "Please provide either a 'mids' object, a data frame, or a list of imputed datasets in the `data` argument.",
      call. = FALSE
    )
  }
  if (is_mids) {
    rlang::check_installed("mice", "mitml")
  }
  return(invisible(NULL))
}
.check_weights <- function(weights) {
  dreamerr::check_arg(weights, "class(devMSM_weights) | list", .message = "The 'weighted' mode of this function requires weights be supplied in the form of output from createWeights.")

  lapply(seq_along(weights), function(i) {
    w <- weights[[i]]
    dreamerr::check_value(
      w, "class(weightitMSM)",
      .arg_name = paste0("weights[[", i, "]]"),
      .message = "Please supply a list of weights output from the createWeights function."
    )
  })
  return(invisible(NULL))
}

# Common functions ----
#' Grabs trailing numbers after the last occurance of `sep` in `vars`
#' @keywords internal
.extract_time_pts_from_vars <- function(vars, sep = "[\\._]") {
  regex_time_pts <- paste0(sep, "([0-9]+)$")
  sapply(vars, function(var) {
    as.numeric(regmatches(var, regexec(regex_time_pts, var))[[1]][2])
  })
}

#' Gets prefix of variable by numbers after last occurance of `sep` and optionally deleting `sep`
#' @keywords internal
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
  if (any(is.na(exposure_time_pts))) {
    stop("Please supply exposure variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }

  tv_conf_time_pts <- .extract_time_pts_from_vars(tv_conf, sep = sep)
  if (any(is.na(tv_conf_time_pts))) {
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

  var_tab <- data.frame(
    var = c(exposure, tv_conf, ti_conf),
    type = c(rep("exposure", length(exposure)), rep("tv_conf", length(tv_conf)), rep("ti_conf", length(ti_conf))),
    time = c(exposure_time_pts, tv_conf_time_pts, rep(-1, length(ti_conf)))
  )
  return(var_tab)
}

# assessBalance ----
#' Used in `assessBalance` to define history
#' @keywords internal
.characterize_exposure <- function(mat, exposure_type) {
  if (exposure_type == "continuous") {
    res <- lapply(mat, function(x) ifelse(x <= median(x), "l", "h"))
  } else {
    res <- lapply(mat, function(x) ifelse(x == 0, "l", "h"))
  }
  do.call(function(...) paste(..., sep = "-"), res)
}

#' Averages balance stats across imputed datasets: df of all times
#' 
#' @keywords internal
.avg_imp_bal_stats <- function(bal_stats) {

  temp <- lapply(bal_stats, function(x){
    b <- do.call(rbind, lapply(x, as.data.frame))
  })
  bal_stats <- do.call(rbind, lapply(bal_stats[[1]], as.data.frame))
  bal_stats$std_bal_stats = rowMeans(do.call(cbind, 
                                             lapply(temp, function(x) abs(x[["std_bal_stats"]]))))
  bal_stats$balanced <- ifelse(bal_stats$std_bal_stats < 
                                 bal_stats$bal_thresh, 1, 0) # recalc balanced

  return(bal_stats)
}

#' Averages balance stats across imputed datasets: list by time
#' 
#' @keywords internal
.avg_imp_bal_stats_time <- function(bal_stats) {

  avgs = lapply(1:5, function(y) {
    rowMeans(do.call(cbind, 
                     lapply(
                       lapply(1:2, function(x){
                         x <- bal_stats[[x]]
                         out <- x[[y]]
                         out
                       }), function(x) abs(x[["std_bal_stats"]]))))
  })
  
  test <- bal_stats[[1]]
  bal_stats <- lapply(1:5, function(x){
    test[[x]]$std_bal_stats <- avgs[[x]]
    test[[x]]
  })
  bal_stats
}


# fitModel ----
#' Takes row means of eposure variables in the same epoch
#' 
#' @keywords internal
.generate_epoch_exposures <- function(data, exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  sp <- split(exposure, epoch)
  res <- lapply(sp, function(epoch_cols) {
    rowMeans(data[, epoch_cols, drop = FALSE], na.rm = TRUE)
  })
  res <- do.call(cbind, res)
  colnames(res) <- paste0(prefix, names(sp))
  return(res)
}

#' Return interactions up to order n of variables in `x`
#' @param x character vector of variable names
#' @param n order of interaction
#' @keywords internal
.get_interaction_terms <- function(x, n = 2) {
  x = sprintf("`%s`", x)
  x = attr(terms(reformulate(x)), "term.labels")
  str = sprintf("~ (%s)^%s", paste(unique(x), collapse = " + "), n)
  interactions <- attr(terms(as.formula(str)), "term.labels")
  return(interactions)
}

# compareHistories ----
#' Simple function to get the variable names associated with epochs 
#' 
#' @keywords internal
.get_epoch_var_names <- function(exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  paste0(prefix, unique(epoch))
}

#' Grab values of epoch variables to evaluate in `marginaleffects::predictions`
#' 
#' @description
#' For continuous variables, grab quantiles of the epoch variables based on hi_lo_cut
#' For binary variables, use 0 and 1.
#' 
#' @keywords internal
.get_avg_predictions_variables <- function(data, epoch_vars, exposure_type, hi_lo_cut = NULL) {
  if (is.null(hi_lo_cut)) hi_lo_cut = c(0.25, 0.75)
  if (exposure_type == "continuous") {
    res = lapply(epoch_vars, function(var_name) {
      qs = quantile(data[[var_name]], probs = hi_lo_cut, na.rm = TRUE)
      if (length(qs) == 1) qs = c(qs - 0.0001, qs + 0.0001)
      return(unname(qs))
    })
    names(res) = epoch_vars
    return(res)
  } else {
    res = lapply(epoch_vars, function(var_name) c(0, 1))
    names(res) = epoch_vars
    return(res)
  }
}

#' Create contrast map using comparison and reference
#' @param histories character vector (from preds$term corresponding to average predictions)
#' @param reference character vector of exposure histories
#' @param comparison character vector of exposure histories
#'
#' @keywords internal
.create_hypotheses_mat <- function(histories, reference, comparison) {
  sub_mats <- lapply(reference, function(ref) {
    n_cols <- length(comparison)
    cus_comps <- matrix(0, ncol = n_cols, nrow = length(histories))
    cus_comps[match(ref, histories), ] <- -1
    cus_comps[match(comparison, histories), 1:n_cols] <- 1
    colnames(cus_comps) <- sprintf("(%s) - (%s)", comparison, ref)
    return(cus_comps)
  })
  do.call(cbind, sub_mats)
}

#' Add history labels and dose to table based on `epoch_vars`
#'
#' @param p table output from marginaleffects::avg_predictions() or hypotheses()
#' @param epoch_vars variable names of exposure epoch (in time order)
#' @return table with histories labeled
#'
#' @keywords internal
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
#' 
#' @keywords internal
.dose_from_history <- function(history, dose_level = "h") {
  .string_count(history, dose_level)
}

#' Calculate number of doses from history string (e.g. "l-l-l")
#' 
#' @examples
#' library(devMSMs)
#' devMSMs:::.dose_from_history_contrast(c("(h-h-h) - (l-l-l)"), dose_level = "h")
#' devMSMs:::.dose_from_history_contrast(c("(h-h-h) - (l-l-l)"), dose_level = "l")
#' 
#' @keywords internal
.dose_from_history_contrast <- function(history, dose_level = "h") {
  sapply(
    strsplit(history, " - "),
    function(hist) {
      counts <- .string_count(hist, dose_level)
      sprintf("%s - %s", counts[1], counts[2])
    }
  )
}


#' Finds custom reference values (currently unused)
#'
#' @param epoch_d data frame of high and low values per exposure main effect
#' @param exp_history sequence of "h" and/or "l" (e.g., "h-h-h")
#' @return exposure history values
#'
#' @examples
#' epoch_d <- data.frame(
#'   e = c("A.1", "A.2", "A.3"),
#'   l = c(0, 0, 0),
#'   h = c(1, 1, 1)
#' )
#' r <- devMSMs:::.get_history_values(
#'   epoch_d = epoch_d,
#'   exp_history = "l-l-l"
#' )
#' r <- devMSMs:::.get_history_values(
#'   epoch_d = epoch_d,
#'   exp_history = "h-h-h"
#' )
#'
#' @keywords internal
.get_history_values <- function(epoch_d, exp_history) {
  sapply(exp_history, function(hist) {
    ifelse(strsplit(hist, "-")[[1]] == "l", epoch_d$l, epoch_d$h)
  })
}




