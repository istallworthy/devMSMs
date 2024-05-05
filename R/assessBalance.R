#' Assesses confounder balancing
#'
#' Draws on functions from the cobalt package to quantify the relations between
#' exposure and confounders at each exposure time point according to the
#' guidelines from Jackson, 2016 on how to assess balance for time-varying
#' exposures.
#'
#' @seealso {[cobalt] package,
#'   <https://cran.r-project.org/web/packages/cobalt/index.html>}
#' @seealso {Jackson, 2016 for more on assessing balance for time-varying
#'   exposures, <https://pubmed.ncbi.nlm.nih.gov/27479649/>}
#'
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param obj initialized MSM object from `initMSM`
#' @param weights list of IPTW weights output from `createWeights`,
#'   required for type 'weighted'
#' @param balance_thresh (optional) one or two numbers between 0 and 1
#'   indicating a single balancing threshold or thresholds for more and less
#'   important confounders, respectively (default = 0.1)
#' @param imp_conf (optional) list of variable names reflecting important
#'   confounders, required if two balance thresholds are supplied
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param home_dir (optional) path to home directory
#'   (required if save.out = TRUE)
#' @returns a data frame of balance statistics
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
#' f <- createFormulas(obj, type = "short")
#' 
#' # Prebalance
#' b <- assessBalance(data = data, obj = obj)
#' print(b)
#' 
#' # returns list of ggplots for each exposure variable
#' plots <- plot(b) 
#' plots[[1]]
#' plots[[2]]
#' plots[[3]]
#' 
#' 
#' # Weighted
#' w <- createWeights(data = data, obj = obj, formulas = f)
#' bw <- assessBalance(data = data, obj = obj, weights = w)
#' print(bw)
#' 
#' @export
assessBalance <- function(
    data, obj, weights = NULL, balance_thresh = NULL, imp_conf = NULL,
    verbose = FALSE, save.out = FALSE, home_dir = NULL) {
  ### Checks ----
  dreamerr::check_arg(verbose, save.out, "scalar logical")
  if (save.out) {
    dreamerr::check_arg_plus(home_dir, "path dir")
    rlang::check_installed(c("knitr", "kableExtra", "stargazer"))
  }

  .check_data(data)
  if (!inherits(obj, "devMSM")) {
    stop("`obj` must be output from `initMSM`", call. = FALSE)
  }
  # dreamerr::check_arg(obj, "class(devMSM)", .message = "`obj` must be output from `initMSM`")

  if (is.null(weights)) {
    type <- "prebalance"
    weights_method <- "no weights"
  } else {
    .check_weights(weights)
    type <- "weighted"
    weights_method <- weights[[1]]$method
  }

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

  if (save.out) {
    .create_dir_if_needed(file.path(home_dir, "balance"))
    .create_dir_if_needed(file.path(home_dir, "balance", type))
    .create_dir_if_needed(file.path(home_dir, "balance", type, "plots"))
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
  all_bal_stats <- lapply(seq_len(m), function(k) {
    if (data_type == "mids") {
      d <- as.data.frame(mice::complete(data, k))
    } else {
      d <- data[[k]]
    }
    w <- weights[[k]]$weights

    return(calcBalStats(
      data = d, obj = obj, weights = w,
      balance_thresh = balance_thresh, imp_conf = imp_conf
    ))
  })

  class(all_bal_stats) <- c("devMSM_bal_stats", "list")
  attr(all_bal_stats, "data_type") <- data_type
  attr(all_bal_stats, "exposure_type") <- attr(obj, "exposure_type")

  is_weighted <- !(is.null(weights))
  attr(all_bal_stats, "weighted") <- is_weighted
  if (is_weighted) {
    attr(all_bal_stats, "weight_method") <- attr(weights, "method")
  }
  return(all_bal_stats)
}

#' @export
print.devMSM_bal_stats <- function(x, i = 1, ...) {
  all_bal_stats <- x[[i]]
  data_type <- attr(x, "data_type")
  weight_method <- attr(x, "weight_method")

  ### Summarizing balance ----
  vars <- unique(unlist(lapply(all_bal_stats, "[[", "covariate")))
  n_covars <- Reduce(`+`, lapply(all_bal_stats, nrow))
  imbalanced_covars <- unlist(lapply(all_bal_stats, function(x) x$covariate[x$balanced == 0]))
  n_imbalanced_covars <- length(imbalanced_covars)

  # Early return, no imbalanced
  if (n_imbalanced_covars == 0) {
    msg <- if (data_type == "imputed") {
      if (!is.null(weight_method)) {
        sprintf("No covariates remain imbalaned for imputation %s using `%s` weighting method.", i, weight_method)
      } else {
        sprintf("No covariates are imbalaned for imputation %s.", i)
      }
    } else {
      if (!is.null(weight_method)) {
        sprintf("No covariates remain imbalaned using `%s` weighting method.", weight_method)
      } else {
        sprintf("No covariates are imbalaned.")
      }
    }
    message(msg)

    return(invisible())
  }

  # TODO: This will break with factor variables e.g. `F.1_b`
  # remaining_imbalanced_domains <- .remove_time_pts_from_vars(unlist(imbalanced_covars))
  imbalanced_std_bal_stats <- unlist(lapply(
    all_bal_stats, function(x) x$std_bal_stats[x$balanced == 0]
  ))
  tab_bal_summary <- do.call(
    rbind,
    lapply(all_bal_stats, function(bal_stats) {
      b <- bal_stats$balanced
      data.frame(
        exposure = bal_stats$exposure[1],
        n_covars = length(b),
        n_imablance_covars = sum(1 - b)
      )
    })
  )

  msg <- if (data_type == "imputed") {
    if (!is.null(weight_method)) {
      sprintf("USER ALERT: For imputation %s using `%s` weighting method:", i, weight_method)
    } else {
      sprintf("USER ALERT: For imputation %s:", i)
    }
  } else {
    if (!is.null(weight_method)) {
      sprintf("USER ALERT: Using `%s` weighting method:", weight_method)
    } else {
      sprintf("USER ALERT:")
    }
  }
  caption <- if (data_type == "imputed") {
    if (!is.null(weight_method)) {
      sprintf("Imbalanced Covariates for imputation %s using `%s`", i, weight_method)
    } else {
      sprintf("Imbalanced Covariates for imputation %s", i)
    }
  } else {
    if (!is.null(weight_method)) {
      sprintf("Imbalanced Covariates using `%s`", weight_method)
    } else {
      sprintf("Imbalanced covariates")
    }
  }

  # TODO:
  # , corresponding to %s out of %s domains,
  # remaining_imbalanced_domains,
  # total_domains,
  message(msg)
  cat(sprintf(
    "As shown below, %s out of %s (%0.1f%%) covariates across time points remain imbalanced with a remaining median absolute value correlation/std mean difference of %.2f (range = %.2f - %.2f):\n",
    n_imbalanced_covars,
    n_covars,
    (n_imbalanced_covars / n_covars) * 100,
    median(abs(imbalanced_std_bal_stats)),
    min(imbalanced_std_bal_stats),
    max(imbalanced_std_bal_stats)
  ))
  colnames(tab_bal_summary) <- c("Exposure", "Total # of covariates", "# of imbalanced covariates")
  print(tinytable::tt(tab_bal_summary, caption = caption), output = "markdown")
}

#' @export
plot.devMSM_bal_stats <- function(x, i = 1, verbose = TRUE, save.out = FALSE, ...) {
  all_bal_stats <- x[[i]]
  data_type <- attr(x, "data_type")
  weight_method <- attr(x, "weight_method")

  exposure_type <- attr(x, "exposure_type")
  x_lab <- if (exposure_type == "continuous") {
    "Correlation with Exposure"
  } else {
    "Standardized Mean Difference Between Exposures"
  }
  k <- ifelse(data_type == "data.frame", 0, i)

  # Loop through exposure variables
  lps <- lapply(all_bal_stats, function(balance_stats) {
    lp <- make_love_plot(
      balance_stats = balance_stats,
      exposure_type = exposure_type,
      k = k,
      weight_method = attr(x, "weight_method")
    )
  })

  for (lp in lps) {
    if (verbose) { 
      plot(lp)
    } 
    
    # TODO: Save plots
    # if (save.out) {
    #   filename <- file.path(
    #     home_dir,
    #     "balance",
    #     folder,
    #     "plots",
    #     sprintf(
    #       "%s_imp_%s_%s_%s_%s_summary_balance_plot.jpeg",
    #       form_name,
    #       k,
    #       exposure,
    #       exposure_time_pt,
    #       weight_method
    #     )
    #   )
    # 
    #   suppressMessages(
    #     ggplot2::ggsave(lp, filename = filename, width = 6, height = 8)
    #   )
    # }
  }

  return(lps)
}
