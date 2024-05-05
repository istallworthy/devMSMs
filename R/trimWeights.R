#' Trim IPTW balancing weights
#'
#' Trims IPTW balancing weights with heavy right tails by populating all weight
#' values above a given quantile with the weight value of that quantile.
#'
#' @seealso [WeightIt::trim()],
#'   <https://search.r-project.org/CRAN/refmans/WeightIt/html/trim.html> which
#'   this function wraps
#' @param weights list of IPTW weights output from createWeights()
#' @inheritParams WeightIt::trim
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @return list of model output with trimmed weights
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
#' tw <- trimWeights(w)
#' print(tw)
#' plot(tw)
#' 
#' trimWeights(w, at = 0.975, lower = TRUE)
#' 
#' 
trimWeights <- function(weights, at = 0, lower = FALSE, verbose = FALSE, save.out = FALSE, home_dir = NULL) {
  ### Checks ----
  dreamerr::check_arg(verbose, save.out, "scalar logical")
  if (save.out) dreamerr::check_arg_plus(home_dir, "path dir")
  
  dreamerr::check_arg(at, "scalar numeric GT(0.5) LT(1) | scalar integer GT(1)")
  dreamerr::check_arg(lower, "scalar logical")
  
  if (!is.null(weights)) .check_weights(weights)
  if (save.out) {
    .create_dir_if_needed(file.path(home_dir, "weights"))
    .create_dir_if_needed(file.path(home_dir, "weights", "values"))
  }

  # No verbose?
  trim_weights <- lapply(weights, WeightIt::trim, at = at, lower = lower)
  
  class(trim_weights) <- c("devMSM_weights", "list")
  attr(trim_weights, "data_type") <- attr(weights, "data_type")
  attr(trim_weights, "method") <- attr(weights, "method")
  attr(trim_weights, "trim") <- attr(trim_weights[[1]][["weights"]], "trim")
  attr(trim_weights, "trim.lower") <- attr(trim_weights[[1]][["weights"]], "trim.lower")

  # TODO: 
  # if (save.out) {
  #   saveRDS(
  #     trim_weights,
  #     file.path(
  #       home_dir, "weights", "values",
  #       sprintf(
  #         "%s-%s_%s_%s_weights_trim.rds",
  #         exposure, outcome, weights[[1]]$method,
  #         quantile
  #       )
  #     )
  #   )
  # }

  if (verbose) print(trim_weights)
  return(trim_weights)
}
