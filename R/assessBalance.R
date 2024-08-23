#' Assesses confounder balancing
#'
#' Draws on functions from the cobalt package to quantify the relations between
#' exposure and confounders at each exposure time point according to the
#' guidelines from Jackson, 2016 on how to assess balance for time-varying
#' exposures.
#'
#' @seealso 
#'  [cobalt] package, <https://cran.r-project.org/web/packages/cobalt/index.html>;
#'  Jackson, 2016 for more on assessing balance for time-varying exposures, 
#'  <https://pubmed.ncbi.nlm.nih.gov/27479649/>
#'
#' @inheritParams devMSM_common_docs
#' @param weights (optional) list of IPTW weights output from [createWeights()]
#' @param balance_thresh (optional) one or two numbers between 0 and 1
#'   indicating a single balancing threshold or thresholds for more and less
#'   important confounders, respectively (default = 0.1)
#' @param imp_conf (optional) list of variable names reflecting important
#'   confounders, required if two balance thresholds are supplied
#'
#' @returns a list containing balance statistics as a dataframe. It is the length
#'  of the number of datasets (1 for a data.frame or the number of imputed datasets)
#'
#' @examples
#' library(devMSMs)
#' data <- data.frame(
#'   ID = 1:50,
#'   A.1 = rnorm(n = 50),
#'   A.2 = rnorm(n = 50),
#'   A.3 = rnorm(n = 50),
#'   B.1 = rnorm(n = 50),
#'   B.2 = rnorm(n = 50),
#'   B.3 = rnorm(n = 50),
#'   C = rnorm(n = 50),
#'   D.3 = rnorm(n = 50)
#' )
#' obj <- initMSM(
#'   data,
#'   exposure = c("A.1", "A.2", "A.3"),
#'   ti_conf = c("C"),
#'   tv_conf = c("B.1", "B.2", "B.3", "D.3")
#' )
#'
#' # Prebalance
#' b <- assessBalance(data = data, obj = obj)
#' print(b)
#'
#' # returns ggplot of balance stats for all exposure variables
#' plots <- plot(b, t = TRUE)
#' # can plot only specific exposure time periods
#' plot(b, t = "A.3")
#' plot(b, t = 3)
#'
#' # Weighted
#' f <- createFormulas(obj, type = "short")
#' w <- createWeights(data = data, formulas = f)
#' bw <- assessBalance(data = data, weights = w)
#' print(bw)
#' plot(bw)
#'
#' @export
assessBalance <- function(
    data, obj, weights = NULL, balance_thresh = NULL, imp_conf = NULL,
    verbose = FALSE, save.out = FALSE) {
  ### Checks ----
  dreamerr::check_arg(verbose, "scalar logical")

	dreamerr::check_arg(save.out, "scalar logical | scalar character")
  
  .check_data(data)
  
  if (is.null(weights)) {
    type <- "prebalance"
    weight_method <- "no weights"
  } else {
    .check_weights(weights)
    type <- "weighted"
    weight_method <- weights[[1]]$method
    obj <- attr(weights, "obj")
  }
  
  if (missing(obj) || !inherits(obj, "devMSM")) {
    stop("`obj` must be output from `initMSM()`", call. = FALSE)
  }
  # dreamerr::check_arg(obj, "class(devMSM)", .message = "`obj` must be output from `initMSM`")
  
  dreamerr::check_set_arg(
    balance_thresh, "numeric vector len(1,2) GT{0} LT{1} NULL{0.1}"
  )
  dreamerr::check_arg(imp_conf, "character vector NULL")
  if (length(balance_thresh) == 2 && is.null(imp_conf)) {
    stop("If you wish to provide different balance threshold for important and less important confounders, please provide a list of important confounders in the 'imp_conf' field.", call. = FALSE)
  }
  if (!is.null(imp_conf) && length(balance_thresh) == 1) {
    stop("If you provide a list of important confounders, please provide a list of two balance thresholds for important and less important confounders, respectively", call. = FALSE)
  }
  
  ### START ----
  if (inherits(data, "mids")) {
    data_type <- "mids"
    m <- data$m
  } else if (inherits(data, "list")) {
    data_type <- "list"
    m <- length(data)
  } else {
    data_type <- "data.frame"
    data <- list(data)
    m <- 1
  }
  
  # calcBalStats ----
  all_bal_stats <- list_of_omitted_histories <- vector("list", m) # store omitted histories for each dataset 
  
  for (k in seq_len(m)) {
    if (data_type == "mids") {
      d <- as.data.frame(mice::complete(data, k))
    } else {
      d <- data[[k]]
    }
    # `weightitMSM` object
    w <- weights[[k]]
    
    res <- calc_bal_stats(
      data = d, obj = obj, weights = w$weights,
      balance_thresh = balance_thresh, imp_conf = imp_conf
    )
    
    if (!is.null(res$omitted_histories)) {
      list_of_omitted_histories[[k]] <- res$omitted_histories
    }
    
    all_bal_stats[[k]] <- res$all_bal_stats
  }

  class(all_bal_stats) <- "devMSM_bal_stats"
  attr(all_bal_stats, "obj") <- obj
  attr(all_bal_stats, "list_of_omitted_histories") <- list_of_omitted_histories

  
  is_weighted <- !(is.null(weights))
  attr(all_bal_stats, "weighted") <- is_weighted
  if (is_weighted) {
    formula_type <- attr(weights, "form_type")
    attr(all_bal_stats, "weight_method") <- attr(weights, "method")
    attr(all_bal_stats, "form_type") <- formula_type
  }
  
  if (verbose) print(all_bal_stats, i = 1, t = TRUE)

  if (isTRUE(save.out) || is.character(save.out)) {
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "balance"))
    .create_dir_if_needed(out_dir)
    
    if (is_weighted){
      ft_str <- sprintf("form_%s", formula_type)
    } else {
      ft_str <- ""
    }

    if (is.character(save.out)) {
      file_name <- save.out
    } else {
      if (is_weighted) { 
        file_name <- sprintf(
          "balance-weighted-exposure_%s-method_%s-%s.rds",
          obj[["exposure_root"]], attr(weights, "method"), ft_str
        )
      } else {
        file_name <- sprintf(
          "prebalance-exposure_%s.rds",
          obj[["exposure_root"]]
        )
      }
    }

    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving balance statistics to `.rds` file. To load, call:\nreadRDS("%s")\n',
      out
    ))
    saveRDS(all_bal_stats, out)
  }
  return(all_bal_stats)
}

#' @rdname assessBalance
#'
#' @param x devMSM_bal_stats object from [assessBalance()]
#' @inheritParams devMSM_common_docs
#' @param ... ignored
#'
#' @export
print.devMSM_bal_stats <- function(x, i = NA, t = TRUE, save.out = FALSE, ...) {

  obj <- attr(x, "obj")
  data_type <- obj[["data_type"]]
  exposure <- obj[["exposure"]]
  exposure_root <- obj[["exposure_root"]]
  is_weighted <- attr(x, "weighted")
  weight_method <- attr(x, "weight_method")
  formula_type <- attr(x, "form_type")
  
  if (data_type %in% c("mids", "list")) data_type <- "imputed"
  else i <- 1L
  
  if (isTRUE(i)) {
    i <- seq_along(x)
  }
  
  if ((length(i) != 1L || !is.na(i)) && !all(i %in% seq_along(x))) {
    stop("`i` must be `TRUE`, `NA`, or an integer indicating which imputation for which to summarize balance.", call. = FALSE)
  }
  
  if (length(i) > 1L) {
    for (j in i) {
      print(x, i = j, t = t, save.out = save.out, ...)
    }
    return(invisible(NULL))
  }
  else if (is.na(i)) {
    all_bal_stats <- .avg_imp_bal_stats_time(x) # TODO: IS notes that this is giving different bal stats than the summary() method does averaging across imps
  } else {
    all_bal_stats <- x[[i]]
  }

  # Process t
  dreamerr::check_arg(t, "integer vector | logical scalar | character vector")
  if (isTRUE(t)) {
    t <- seq_along(all_bal_stats)
    t_str <- "All Exposure Time Points"
  } else if (is.character(t)) {
    t <- match(t, exposure)
    t_str = sprintf("Exposure Time %s", t)
  }
  else {
    t <- match(t, seq_along(exposure))
    t_str = sprintf("Exposure Time %s", t)
  }
  
  if (anyNA(t)) {
    stop("`t` must be `TRUE` or the names or indices of exposures.", call. = FALSE)
  }
  
  if (data_type == "data.frame") {
    tbl_caption <- sprintf("Balance Stats for %s", t_str)
  } else if (!is.na(i)) {
    tbl_caption <- sprintf("Balance Stats for %s for imputation (%s)", t_str, i)
  } else {
    tbl_caption <- sprintf("Balance Stats for %s Averaging Across Imputed Datasets", t_str)
  }
  
  balance_stats <- do.call("rbind", all_bal_stats[t])
  balance_stats <- balance_stats[, c("exposure", "covariate", "std_bal_stats", "bal_thresh", "balanced")]
  b <- tinytable::tt(balance_stats, digits = 3, caption = tbl_caption)

  print(b, "markdown")

  if (isTRUE(save.out) || is.character(save.out)) {
    rlang::check_installed("pandoc")
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "balance"))
    .create_dir_if_needed(out_dir)

    if (data_type == "data.frame") {
      i_str <- ""
    } else if (!is.na(i)) {
      i_str <- sprintf("-imp_%s", i)
    } else {
      i_str <- "-imp_averaged"
    }

    if (is_weighted){
      ft_str <- sprintf("form_%s", formula_type)
    } else {
      ft_str <- ""
    }
    
    if (is.character(save.out)) {
      file_name <- save.out
    } else if (is_weighted) { 
      file_name <- sprintf(
        "bal_stats_table_weighted-exposure_%s-method_%s-form_%s%s-%s.docx",
        exposure_root, weight_method, ft_str, i_str, gsub(" ", "_", t_str)
      )
    } else {
      file_name <- sprintf(
        "bal_stats_table_prebalance-exposure_%s%s-%s.docx",
        exposure_root, i_str, gsub(" ", "_", t_str)
      )
    }

    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving balance statistics table to file:\n%s\n',
      out
    ))
    if (fs::path_ext(out) == "pdf") {
      tinytable::save_tt(
        tinytable::format_tt(b, escape = TRUE),
        output = out, overwrite = TRUE
      )
    } else {
      tinytable::save_tt(b, output = out, overwrite = TRUE)
    }
  }

  return(invisible(b))

  # TODO: Print omitted histories?
  # if (verbose) {
  #   if (data_type == "imputed") {
  #     cat(sprintf(
  #       "USER ALERT: the following history/histories, %s has/have been omitted from balance checking for exposure %s imputation %s at time point %s due to insufficient counts:",
  #       omitted_histories, exposure, k, exposure_time_pt
  #     ))
  #     cat("\n")
  #   } else {
  #     cat(sprintf(
  #       "USER ALERT: the following history/histories, %s has/have been omitted from balance checking for exposure %s at time point %s due to insufficient counts:",
  #       omitted_histories, exposure, exposure_time_pt
  #     ))
  #     cat("\n")
  #   }
  # }
}

#' @rdname assessBalance
#'
#' @param object devMSM_bal_stats object from [assessBalance()]
#' @inheritParams devMSM_common_docs
#' @param ... ignored
#'
#' @export
summary.devMSM_bal_stats <- function(object, i = NA, save.out = FALSE, ...) {
  obj <- attr(object, "obj")
  data_type <- obj[["data_type"]]
  exposure_root <- obj[["exposure_root"]]
  weight_method <- attr(object, "weight_method")
  is_weighted <- attr(object, "weighted")
  formula_type <- attr(object, "form_type")
  
  if (data_type %in% c("mids", "list")) data_type <- "imputed"
  else i <- 1L
  
  if (isTRUE(i)) {
    i <- 1L
  }
  
  if (length(i) != 1L || !(is.na(i) || (i %in% seq_along(object)))) {
    stop("`i` must be `NA` or an integer indicating which imputation for which to summarize balance.", call. = FALSE)
  }
  
  if (is.na(i)) {
    all_bal_stats <- .avg_imp_bal_stats_time(object)
  } else {
    all_bal_stats <- object[[i]]
  }
  
  
  
  ### Summarizing balance ----
  n_covars <- Reduce(`+`, lapply(all_bal_stats, nrow))
  imbalanced_covars <- unlist(lapply(all_bal_stats, function(object) object$covariate[object$balanced == 0]))
  n_imbalanced_covars <- length(imbalanced_covars)
  
  # Early return, no imbalanced
  if (n_imbalanced_covars == 0) {
    msg <- {
      if (data_type == "imputed") {
        if (!is.null(weight_method)) {
          if (!is.na(i)) {
            sprintf("No covariates remain imbalaned for imputation %s using `%s` weighting method.", i, weight_method)
          } else {
            sprintf("Averaging across imputed datasets using %s weighting method.", weight_method)
          }
        } else {
          if (!is.na(i)) {
            sprintf("No covariates are imbalanced for imputation %s.", i)
          } else {
            "No covariates are imbalanced averaging across imputed datasets."
          }
        }
      } else {
        if (!is.null(weight_method)) {
          sprintf("No covariates remain imbalaned using `%s` weighting method.", weight_method)
        } else {
          "No covariates are imbalaned."
        }
      }
    }
    
    message(msg)
  }
  
  imbalanced_std_bal_stats <- unlist(lapply(
    all_bal_stats, function(object) object$std_bal_stats[object$balanced == 0]
  ))
  
  tab_bal_summary <- do.call(
    "rbind",
    lapply(all_bal_stats, function(bal_stats) {
      b <- bal_stats$balanced
      data.frame(
        exposure = bal_stats$exposure[1],
        n_covars = length(b),
        n_imablance_covars = sum(1 - b)
      )
    })
  )
  colnames(tab_bal_summary) <- c("Exposure", "Total # of covariates", "# of imbalanced covariates")
  
  msg <- {
    if (data_type == "imputed") {
      if (!is.null(weight_method)) {
        if (!is.na(i)) {
          sprintf("USER ALERT: For imputation %s using `%s` weighting method:", i, weight_method)
        } else {
          sprintf("USER ALERT: Averaging across imputed datasets using `%s` weighting method:", weight_method)
        }
      } else {
        if (!is.na(i)){
          sprintf("USER ALERT: For imputation %s:", i)
        } else {
          sprintf("USER ALERT: Averaging across imputated datasets:")
        }
      }
    } else {
      if (!is.null(weight_method)) {
        sprintf("USER ALERT: Using `%s` weighting method:", weight_method)
      } else {
        ""
      }
    }
  }
  
  main_msg <- {
    if (n_imbalanced_covars == 0) {
      sprintf(
        "As shown below, %s out of %s (%0.1f%%) covariates across time points remain imbalanced:",
        n_imbalanced_covars,
        n_covars,
        (n_imbalanced_covars / n_covars) * 100
      )
    } else {  
      sprintf(
        "As shown below, %s out of %s (%0.1f%%) covariates across time points remain imbalanced with a remaining median absolute %s of %.2f (max: %.2f):",
        n_imbalanced_covars,
        n_covars,
        (n_imbalanced_covars / n_covars) * 100,
        switch(obj[["exposure_type"]], "continuous" = "correlation", "standardized mean difference"),
        median(abs(imbalanced_std_bal_stats)),
        max(abs(imbalanced_std_bal_stats))
      )
    }}
  
  caption <- {
    if (data_type == "imputed") {
      if (!is.null(weight_method)) {
        if (!is.na(i)) {
          sprintf("Imbalanced Covariates for imputation %s using `%s`", i, weight_method)
        } else {
          sprintf("Imbalanced Covariates Averaging Across Imputed Datasetrs using `%s`", weight_method)
        }
      } else {
        if (!is.na(i)) {
          sprintf("Imbalanced Covariates for imputation %s", i)
        } else {
          "Imbalanced Covariates Averaged Across Imputed Datasets"
        }
      }
    } else {
      if (!is.null(weight_method)) {
        sprintf("Imbalanced Covariates using `%s`", weight_method)
      } else {
        "Imbalanced covariates"
      }
    }
  }
  t <- tinytable::tt(tab_bal_summary, caption = caption)
  
  cat(msg, main_msg)
  print(t, output = "markdown")
  
  if (isTRUE(save.out) || is.character(save.out)) {
    rlang::check_installed("pandoc")
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "balance"))
    .create_dir_if_needed(out_dir)
    
    if (data_type == "data.frame") {
      i_str <- ""
    } else if (!is.na(i)) {
      i_str <- sprintf("-imp_%s", i)
    } else {
      i_str <- "-imp_averaged"
    }
    
    if (is_weighted){
      ft_str <- sprintf("form_%s", formula_type)
    } else {
      ft_str <- ""
    }
    
    if (is.character(save.out)) {
      file_name <- save.out
    } else if (is_weighted) { 
      file_name <- sprintf(
        "summ_bal_stats_table_weighted-exposure_%s-method_%s-%s%s.docx",
        exposure_root, weight_method, ft_str, i_str
      )
    } else {
      file_name <- sprintf(
        "summ_bal_stats_table_prebalance-exposure_%s%s.docx",
        exposure_root, i_str
      )
    }
    
    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving balance statistics table to file:\n%s\n',
      out
    ))
    if (fs::path_ext(out) == "pdf") {
      tinytable::save_tt(
        tinytable::format_tt(t, escape = TRUE),
        output = out, overwrite = TRUE
      )
    } else {
      tinytable::save_tt(t, output = out, overwrite = TRUE)
    }
  }
  
  return(invisible(t))
}


#' @rdname assessBalance
#'
#' @param x devMSM_bal_stats object from `assessBalance`
#' @inheritParams devMSM_common_docs
#' @param ... ignored
#'
#' @export
plot.devMSM_bal_stats <- function(x, i = NA, t = TRUE, save.out = FALSE, ...) {

  obj <- attr(x, "obj")
  data_type <- obj[["data_type"]]
  exposure <- obj[["exposure"]]
  exposure_root <- obj[["exposure_root"]]
  exposure_type <- obj[["exposure_type"]]
  exposure_time <- obj[["exposure_time_pts"]]
  is_weighted <- obj[["weighted"]]
  weight_method <- attr(x, "weight_method")
  is_weighted <- !is.null(weight_method)
  formula_type <- attr(x, "form_type")
  
  
  if (data_type %in% c("mids", "list")) data_type <- "imputed"
  else i <- 1L

  if (isTRUE(i)) {
    i <- seq_along(x)
  }

  ## IS added to print avgs over impts for i = NA
  if ((length(i) != 1L || !is.na(i)) && !all(i %in% seq_along(x))) {
    stop("`i` must be `NA` or an integer indicating which imputation for which to summarize balance.", call. = FALSE)
  }
  
  if (length(i) > 1L) {
    for (j in i) {
      print(plot(x, i = j, t = t, save.out = save.out, ...))
    }
    return(invisible(NULL))
  }
  else if (is.na(i)) {
    all_bal_stats <- .avg_imp_bal_stats_time(x) 
  } else {
    all_bal_stats <- x[[i]]
  }

  k <- switch(data_type, "imputed" = i, 0)
  
  # Process t
  dreamerr::check_arg(t, "integer vector | logical scalar | character vector")
  if (isTRUE(t)) {
    t <- seq_along(all_bal_stats)
  } else if (is.character(t)) {
    t <- match(t, exposure)
  }
  else {
    t <- match(t, seq_along(exposure))
  }
  
  if (anyNA(t)) {
    stop("`t` must be `TRUE` or the names or indices of exposures.", call. = FALSE)
  }
  
  balance_stats <- do.call("rbind", all_bal_stats[t])
  
  # Loop through selected exposure variables
  lp <- make_love_plot(
    obj = obj,
    balance_stats = balance_stats,
    k = k,
    t = t,
    weight_method = attr(x, "weight_method")
  )

  if (isTRUE(save.out) || is.character(save.out)) {
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "balance", "plots"))
    .create_dir_if_needed(out_dir)
    
    # timing title
    if (length(t) > 1){
      t_t <- "t=all"
    } else {
      t_t <- sprintf("t=%s", exposure_time[t])
    }
    
    if (is_weighted){
      ft_str <- sprintf("form_%s", formula_type)
    } else {
      ft_str <- ""
    }

    if (data_type == "data.frame") {
      i_str <- ""
    } else if (!is.na(i)) {
      i_str <- sprintf("-imp_%s", i)
    } else {
      i_str <- "-imp_averaged"
    }

    if (is.character(save.out)) {
      file_names <- save.out
    } else if (is_weighted) { 
      file_names <- sprintf(
        "weighted-exposure_%s-method_%s-%s%s-%s.png",
        exposure_root, weight_method, ft_str, i_str, t_t
      )
    } else {
      file_names <- sprintf(
        "prebalance-exposure_%s%s-%s.png",
        exposure_root, i_str, t_t
      )
    }
    
    out <- fs::path_join(c(out_dir, file_names))
    cat(sprintf(
      "\nSaving love-plot of balance statistics:\n%s\n",
      out
    ))

    ggplot2::ggsave(out, plot = lp, height = 8, width = 14)
  }

  return(lp)
}
