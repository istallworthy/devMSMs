#' Trim IPTW balancing weights, if needed
#'
#' Trims IPTW balancing weights with heavy right tails by populating all weight
#' values above a given quantile with the weight value of that quantile.
#'
#' @seealso [WeightIt::trim()],
#'   <https://search.r-project.org/CRAN/refmans/WeightIt/html/trim.html> which
#'   this function wraps
#' 
#' @inheritParams devMSM_common_docs
#' @inheritParams WeightIt::trim
#' 
#' @return a list containing [WeightIt::weightitMSM()] output. It is the length 
#'  of the number of datasets (1 for a data.frame or the number of imputed datasets).
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
#' f <- createFormulas(obj, type = "short")
#' 
#' w <- createWeights(data = data, obj = obj, formulas = f)
#' tw <- trimWeights(obj, w)
#' print(tw)
#' plot(tw)
#' 
#' trimWeights(obj, w, at = 0.975, lower = TRUE)
#' 
#' @export
trimWeights <- function(obj, weights, at = 0, lower = FALSE, verbose = FALSE, save.out = FALSE) {
  ### Checks ----
  dreamerr::check_arg(verbose, "scalar logical")

	dreamerr::check_arg(save.out, "scalar logical | scalar character")
  dreamerr::check_arg(at, "scalar numeric GT(0.5) LT(1) | scalar integer GT(1)")
  dreamerr::check_arg(lower, "scalar logical")
  
  if (!is.null(weights)) .check_weights(weights)

  # No verbose?
  trim_weights <- lapply(weights, WeightIt::trim, at = at, lower = lower)
  
  class(trim_weights) <- c("devMSM_weights", "list")
  attr(trim_weights, "obj") <- obj
  attr(trim_weights, "method") <- attr(weights, "method")
  attr(trim_weights, "form_type") <- attr(weights, "form_type")
  attr(trim_weights, "trim") <- attr(trim_weights[[1]][["weights"]], "trim")
  attr(trim_weights, "trim_lower") <- attr(trim_weights[[1]][["weights"]], "trim.lower")

  if (verbose) print(trim_weights)
  
  if (save.out == TRUE || is.character(save.out)) {
    home_dir <- attr(obj, "home_dir")
    out_dir <- fs::path_join(c(home_dir, "weights"))
    .create_dir_if_needed(out_dir)

    if (is.character(save.out)) {
      file_name <- save.out
    } else {
      file_name <- sprintf(
        "type_%s-exposure_%s-method_%s-trim_at_%s-lower_%s.rds",
        attr(weights, "form_type"), 
        attr(obj, "exposure_root"),  
        attr(weights, "method"), 
        at, 
        tolower(lower)
      )
    }
    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving (trimmed) weights to `.rds` file. To load, call:\nreadRDS("%s")\n',
      out
    ))
    saveRDS(weights, out)
  }

  return(trim_weights)
}
