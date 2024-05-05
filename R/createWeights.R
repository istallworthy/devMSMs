#' Creates IPTW balancing weights
#'
#' Creates IPTW balancing weights at each user-specified exposure time point
#' using balancing formulas that relate exposure at each time point to all
#' relevant confounders.
#'
#' @export
#' @seealso {[WeightIt::weightitMSM()],
#'   <https://ngreifer.github.io/WeightIt/reference/weightitMSM.html>}
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param obj initialized MSM object from `initMSM`
#' @param formulas list of balancing formulas at each time point output from
#'   `createFormulas()`
#' @param method (optional) character string of weightitMSM() balancing method
#'   abbreviation (default is Covariate Balancing Propensity Score "cbps")
#' @param ... pass custom arguments to `weightitMSM()`
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   FALSE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is FALSE)
#' @param home_dir path to home directory (required if 'save.out' = TRUE)
#' @return list of IPTW balancing weights
#' @export
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
#' w <- createWeights(data = data, obj = obj, formulas = f)
#' print(w)
#' plot(w)
#'
#' # Methods from `WeightIt::weightitMSM`
#' @examplesIf requireNamespace("CBPS", quietly = TRUE)
#' w <- createWeights(data = data, obj = obj, formulas = f, method = "cbps")
#'
#' @examplesIf requireNamespace("gbm", quietly = TRUE)
#' w <- createWeights(data = data, obj = obj, formulas = f, method = "gbm")
#'
#' @examplesIf requireNamespace("dbarts", quietly = TRUE)
#' w <- createWeights(data = data, obj = obj, formulas = f, method = "bart")
#'
#' @examplesIf requireNamespace("SuperLearner", quietly = TRUE)
#' w <- createWeights(data = data, obj = obj, formulas = f, method = "super")
#'
createWeights <- function(
    data, obj, formulas,
    method = c("glm", "gbm", "bart", "super", "cbps"),
    verbose = FALSE, save.out = FALSE, home_dir = NULL,
    ...) {
  ### Checks ----
  dreamerr::check_arg(verbose, save.out, "scalar logical")
  if (save.out) dreamerr::check_arg_plus(home_dir, "path dir")

  .check_data(data)
  dreamerr::check_arg(formulas, "class(devMSM_formulas)")
  method <- match.arg(method, several.ok = FALSE)
  form_type <- attr(formulas, "type")

  if (save.out) {
    .create_dir_if_needed(file.path(home_dir, "weights"))

    folder_dir <- file.path(home_dir, "weights", form_type)
    .create_dir_if_needed(folder_dir)

    plots_dir <- file.path(home_dir, "weights", form_type, "plots")
    .create_dir_if_needed(plots_dir)
  }

  ### Create weights ----
  dots <- list(...)

  calculate_weights <- function(formulas, data, method, verbose, dots) {
    # TODO: Write `check_formula_data_for_na` to only look at variables in the formula
    if (anyNA(data)) {
      stop("This code requires complete data. Consider imputation if missingness < 20% and is reasonably Missing at Random (MAR).", call. = FALSE)
    }

    args <- list(
      formula.list = formulas,
      data = data,
      method = method,
      verbose = verbose, 
      include.obj = TRUE
      # stabilize = TRUE,
    )
    
    # TODO: Discuss defaults
    custom_args <- switch(method,
      "super" = list(SL.library = "SL.glm"),
      "glm" = list(use.kernel = TRUE),
      "gbm" = list(use.kernel = TRUE, criterion = "p.mean"), # required? even tho doc says there is default
      "cbps" = list(use.kernel = TRUE, over = FALSE),
      "bart" = list(use.kernel = TRUE, over = FALSE)
    )
    args <- modifyList(args, custom_args)

    # overwrite options with captured `...` from `createWeights`
    args <- modifyList(args, dots, keep.null = TRUE)

    # do.call(WeightIt::weightitMSM, args)
    do.call(WeightIt::weightitMSM, args = args)
  }

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

  weights <- lapply(seq_len(m), function(k) {
    if (data_type == "mids") {
      d <- as.data.frame(mice::complete(data, k))
    } else {
      d <- data[[k]]
    }
    calculate_weights(formulas = formulas, data = d, method = method, verbose = verbose, dots = dots)
  })

  ### Output ----
  class(weights) <- c("devMSM_weights", "list")
  attr(weights, "data_type") <- data_type
  attr(weights, "method") <- method
  attr(weights, "form_type") <- form_type

  return(weights)
}


# TODO: Print for each x? `all = FALSE` default, but print them
#' @export
print.devMSM_weights <- function(x, i = 1, ...) {
  w <- x[[i]]$weights

  trim <- attr(x, "trim")
  trim_lower <- attr(x, "trim")
  data_type <- attr(x, "data_type")
  method <- attr(x, "method")

  trim_str <- .get_trim_str(trim, trim_lower)
  if (trim_str != "") trim_str <- paste0(trim_str, ", ")

  list_str <- ""
  if (data_type == "mids" | data_type == "list") {
    list_str <- sprintf("imputation %s and ", i)
  }

  msg <- sprintf(
    "For %sthe \`%s\` weighting method, %sthe median weight value is %s (SD = %s; range = %s-%s).",
    list_str,
    method,
    trim_str,
    round(median(w), 2),
    round(sd(w), 2),
    round(min(w), 2),
    round(max(w))
  )

  cat(msg)
  return(invisible(NULL))
}

#' @export
plot.devMSM_weights <- function(x, i = 1, ...) {
  w <- x[[i]]$weights
  trim <- attr(x, "trim")
  trim_lower <- attr(x, "trim")
  method <- attr(x, "method")

  trim_str <- .get_trim_str(trim, trim_lower)

  if (is.null(trim)) {
    title <- sprintf("Distribution of %s weights", method)
  } else {
    title <- sprintf(
      "Distribution of %s weights (%s)",
      method, trim_str
    )
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(aes(x = w), color = "black", bins = 15) +
    ggplot2::labs(
      title = title, x = "Weights"
    )
  return(p)
}

.get_trim_str <- function(trim, trim_lower) {
  trim_str <- ""
  if (!is.null(trim)) {
    if (trim >= 1) {
      if (trim_lower == TRUE) {
        trim_str <- sprintf("after trimming the top and bottom %s weights", trim)
      } else {
        trim_str <- sprintf("after trimming the top %s weights", trim)
      }
    } else {
      if (trim_lower == TRUE) {
        trim_str <- sprintf("after trimming between %sth and %sth quantiles", round((1 - trim) * 1000) / 10, round(trim * 1000) / 10)
      } else {
        trim_str <- sprintf("after trimming at %sth quantile", round(trim * 1000) / 10)
      }
    }
  }

  return(trim_str)
}
