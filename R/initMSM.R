#' Initial step in `devMSMs` workflow
#'
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or `mids` object from the `mice` package
#' @param exposure names of exposure variables with ".timepoint" suffix
#' @param epoch (optional) group set of exposure variables into categories.
#'  Provide a character vector corresponding to the category for each exposure
#' @param tv_conf list of time-varying confounders with ".timepoint"
#'   suffix, should include exposure and outcome variables (at least
#'   time-varying exposure variables required here)
#' @param ti_conf list of time invariant confounders. Can be left as NULL for none.
#' @param concur_conf (optional) list of variable names reflecting time-varying
#'   confounders to retain in formulas contemporaneously (default is none)
#' @param home_dir (optional) directory for saving output. Either an absolute path or a relative path with respect to `getwd()`
#' @param sep (optional) The seperator between the variable and the time period. The
#'   variable names will be split by the last occurance of `sep` with the
#'   second string containing the time. This uses regex notation, so `.` must
#'   be `\\.`
#'
#' @details
#' By `.timepoint` suffix, we mean that time-varying and exposure variable
#' names must end with either `.#` or a `_#`. This allows us to extract the
#' time-period when the variable is measured and will allow us to properly
#' create the formulae (omitting future mediators)
#'
#' @return object of class `devMSM` that contains the initialized information.
#'
#' @examples
#' data <- data.frame(
#'   A.1 = rnorm(n = 50),
#'   A.2 = rnorm(n = 50),
#'   A.3 = rnorm(n = 50),
#'   B.1 = rnorm(n = 50),
#'   B.2 = rnorm(n = 50),
#'   B.3 = rnorm(n = 50),
#'   D.3 = rnorm(n = 50),
#'   L.1 = sample(c(0, 1), size = 50, replace = TRUE),
#'   C   = rnorm(n = 50)
#' )
#' obj <- initMSM(
#'   data = data,
#'   exposure = c("A.1", "A.2", "A.3"),
#'   tv_conf = c("B.1", "B.2", "B.3", "D.3"),
#'   ti_conf = "C"
#' )
#'
#' obj
#'
#' @export
initMSM <- function(data, exposure, epoch = NULL, tv_conf, ti_conf = NULL,
                    concur_conf = NULL, home_dir = NULL, sep = "[\\._]") {
  if (inherits(data, "mids")) {
    d <- data[[1]]
    data_type <- "mids"
  } else if (is.data.frame(data)) {
    d <- data
    data_type <- "data.frame"
  } else if (is.list(data)) {
    d <- data[[1]]
    data_type <- "list"
  }
  else {
    stop("`data` must be a data.frame, a list of data.frames, or a `mids` object (from `mice::mice()`).",
         call. = FALSE)
  }

  dreamerr::check_arg(exposure, "vector character len(1, )")
  dreamerr::check_value(
    reformulate(exposure), "formula var(data)",
    .data = d, .arg_name = "exposure"
  )
  dreamerr::check_set_arg(
    epoch, "NULL{exposure} | vector character len(data)",
    .data = exposure
  )
  dreamerr::check_arg(tv_conf, "vector character len(1, )")
  dreamerr::check_arg(ti_conf, "vector character | NULL")
  dreamerr::check_value(
    reformulate(tv_conf), "formula var(data)",
    .data = d, .arg_name = "tv_conf"
  )
  if (!is.null(ti_conf)) {
    dreamerr::check_value(
      reformulate(ti_conf), "formula var(data)",
      .data = d, .arg_name = "ti_conf"
    )
  }
  if (!is.null(concur_conf)) {
    dreamerr::check_value(
      reformulate(concur_conf), "formula var(data)",
      .data = d, .arg_name = "concur_conf"
    )
    concur_in_tv <- concur_conf %in% tv_conf
    if (!all(concur_in_tv)) {
      stop(sprintf("The following `concur_conf` do not show up in `tv_conf` but need to:%s",
                   paste(concur_conf[!concur_in_tv], collapse = ", ")), call. = FALSE)
    }
  }
  if (any(exposure %in% tv_conf)) {
    stop("`exposure` must not be in `tv_conf`.", call. = FALSE)
  }

  dreamerr::check_arg_plus(home_dir, "path dir | NULL")
  # Make absolute path
  cwd <- getwd()
  if (is.null(home_dir)) {
    home_dir <- cwd
  } else if (!fs::is_absolute_path(home_dir)) {
    home_dir <- fs::path_join(c(cwd, home_dir))
  }
  home_dir <- fs::path_norm(home_dir)
  .create_dir_if_needed(home_dir)

  # checking exposure type
  for (e in exposure) {
    x <- d[[e]]
    if (is.integer(x) && !all(na.omit(unique(x)) %in% c(1, 0))) {
      stop("Please make sure your exposure levels are 1s and 0s for integer exposures.", call. = FALSE)
    }
    if (!is.numeric(x) && !is.integer(x)) {
      stop("Please provide an exposure in numeric or integer form.", call. = FALSE)
    }
  }
  exp1 <- d[[exposure[1]]]
  exposure_type <- if (is.logical(exp1) || all(exp1 %in% c(0, 1), na.rm = TRUE)) "binary" else "continuous"

  # Create var_tab containing all the information needed
  var_tab <- .create_var_tab(exposure = exposure, ti_conf = ti_conf, tv_conf = tv_conf, sep = sep)
  exposure_time_pts <- var_tab$time[var_tab$type == "exposure"]
  rownames(var_tab) <- NULL

  sort <- order(exposure_time_pts)
  exposure <- exposure[sort]
  exposure_time_pts <- exposure_time_pts[sort]
  epoch <- epoch[sort]

  if (!is.null(concur_conf)) {
    which_concur <- (var_tab$var %in% concur_conf)
    var_tab$concur_conf <- which_concur
  }

  # return initialized object
  obj <- list(
    data_type = data_type,
    exposure = exposure,
    epoch = epoch,
    exposure_root = .remove_time_pts_from_vars(exposure[1], sep = sep),
    exposure_time_pts = exposure_time_pts,
    exposure_type = exposure_type,
    sep = sep,
    home_dir = home_dir,
    var_tab = var_tab
  )
  
  class(obj) <- "devMSM"

  return(obj)
}

#' @export
print.devMSM <- function(x, ...) {
  exposure <- x[["exposure"]]
  epoch <- x[["epoch"]]
  cat(sprintf(
    "Exposure (%s): %s\n",
    x[["exposure_type"]],
    paste(exposure, collapse = ", ")
  ))
  if (!all(exposure == epoch)) {
    cat(sprintf(
      "Corresponding epoch: %s\n",
      paste(epoch, collapse = ", ")
    ))
  }
  cat(sprintf("Variable and their encodings:\n"))
  print(x[["var_tab"]], row.names = FALSE)
  
  invisible(x)
}
