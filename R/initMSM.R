# TODO: add mice support to `initMSM`

#' Initial step in `devMSMs` workflow
#'
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param exposure names of exposure variables with ".timepoint" suffix
#' @param tv_conf list of time-varying confounders with ".timepoint"
#'   suffix, should include exposure and outcome variables (at least
#'   time-varying exposure variables required here)
#' @param ti_conf list of time invariant confounders (at least one required)
#' 
#' @details
#' By `.timepoint` suffix, we mean that time-varying and exposure variable 
#' names must end with either `.#` or a `_#`. This allows us to extract the 
#' time-period when the variable is measured and will allow us to properly
#' create the formulae (omitting future mediators)
#'
#' @return object of class `devMSM` that contains the initialized information.
#' @export
initMSM <- function(data, exposure, tv_conf, ti_conf) {
  dreamerr::check_arg(tv_conf, exposure, "vector character len(2, )")
  dreamerr::check_arg(ti_conf, "vector character len(1, ) | NULL")
  dreamerr::check_value(
    reformulate(exposure), "formula var(data)",
    .data = data, .arg_name = "exposure"
  )
  dreamerr::check_value(
    reformulate(tv_conf), "formula var(data)",
    .data = data, .arg_name = "tv_conf"
  )
  dreamerr::check_value(
    reformulate(ti_conf), "formula var(data)",
    .data = data, .arg_name = "ti_conf"
  )

  if (any(exposure %in% tv_conf)) {
    stop("`exposure` must not be in `tv_conf`.", call. = FALSE)
  }

  # checking exposure type
  lapply(seq_along(exposure), function(i) {
    x = data[[exposure[i]]]
    if (is.integer(x) && !all(unique(x) %in% c(1, 0))) {
      stop("Please make sure your exposure levels are 1s and 0s for integer exposures.", call. = FALSE)
    }
    if (!is.numeric(x) && !is.integer(x)) {
      stop("Please provide an exposure in numeric or integer form.", call. = FALSE)
    }
  })
  exposure_type <- ifelse(
    is.numeric(data[[exposure[1]]]), 
    "continuous", "binary"
  )

  # Create var_tab containing all the information needed
  var_tab = .create_var_tab(exposure = exposure, ti_conf = ti_conf, tv_conf = tv_conf)
  exposure_time_pts <- var_tab$time[var_tab$type == "exposure"]
  rownames(var_tab) <- NULL

  # return initialized object
  obj <- list()
  class(obj) <- c("devMSM", "list")
  attr(obj, "var_tab") <- var_tab
  attr(obj, "exposure") <- exposure
  attr(obj, "exposure_time_pts") <- exposure_time_pts
  # TODO: Can this change by exposure?
  attr(obj, "exposure_type") <- exposure_type

  return(obj)
}

#' @export 
print.devMSM = function(x, ...) {
  cat(sprintf(
    "Exposure (%s): %s\n", 
    attr(x, "exposure_type"),
    paste(attr(x, "exposure"), collapse = ", ")
  ))
  cat(sprintf("Variable and their encodings:\n"))
  print(attr(x, "var_tab"), row.names = FALSE)
}


