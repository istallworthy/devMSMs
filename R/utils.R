.extract_time_pts_from_vars <- function(vars, regex_time_pts = "[\\.\\_]([0-9]+)$") {
  sapply(vars, function(var) {
    as.numeric(regmatches(var, regexec(regex_time_pts, var))[[1]][2])
  })
}
.remove_time_pts_from_vars <- function(vars, regex_time_pts = "[\\.\\_]([0-9]+)$") {
  sapply(vars, function(var) {
    gsub(regex_time_pts, "", var)
  }, USE.NAMES = FALSE)
}
.create_var_tab <- function(exposure, ti_conf, tv_conf) {
  dreamerr::check_arg(tv_conf, exposure, "vector character len(2, )")
  dreamerr::check_arg(ti_conf, "vector character len(1, ) | NULL")

  exposure_time_pts = .extract_time_pts_from_vars(exposure)
  if (any(is.na(exposure_time_pts))) {
    stop("Please supply exposure variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }
  
  tv_conf_time_pts = .extract_time_pts_from_vars(tv_conf)
  if (any(is.na(tv_conf_time_pts))) {
    stop("Please supply tv_conf variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }

  ti_conf_time_pts = .extract_time_pts_from_vars(ti_conf)
  time_pt_detected = !is.na(ti_conf_time_pts)
  if (any(time_pt_detected)) {
    warning(sprintf(
      "Potentially time-varying covariate was passed to ti_conf:\n  %s.\nThis is likely because it ends with a seperator and a number. If this variable is exogenous or time-invariant, please ignore.",
      paste(ti_conf[time_pt_detected], sep = ", ")
    ))
  }
  
  var_tab = rbind(
    data.frame(var = exposure, type = "exposure", time = exposure_time_pts),
    data.frame(var = tv_conf,  type = "tv_conf",  time = tv_conf_time_pts),
    data.frame(var = ti_conf,  type = "ti_conf",  time = -1)
  ) 
  return(var_tab)
}

.create_dir_if_needed <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# Used in `assessBalance` to define history
.characterize_exposure = function(mat, exposure_type) {
  if (exposure_type == "continuous") {
    res = lapply(mat, function(x) as.numeric(x <= median(x)))
  } else {
    res = lapply(mat, function(x) as.numeric(x == 0))
  }
  do.call(function(...) paste(..., sep = "-"), res)
}

# Checks
.check_data <- function(data) {
  is_mids = inherits(data, "mids")
  is_df = is.data.frame(data)
  is_list_of_df = is.list(data) && all(sapply(data, is.data.frame))
  if (!(is_mids | is_df | is_list_of_df)) {
    stop(
      "Please provide either a 'mids' object, a data frame, or a list of imputed datasets.",
      call. = FALSE
    )
  }
  if (is_mids) rlang::check_installed("mice")
}

.check_weights <- function(weights) {
  # dreamerr::check_arg(weights, "class(devMSM_weights) | list", .message = "The 'weighted' mode of this function requires weights be supplied in the form of output from createWeights.")

  lapply(seq_along(weights), function(i) {
    w <- weights[[i]]
    # dreamerr::check_value(
    #   w, "class(weightitMSM)",
    #   .arg_name = paste0("weights[[", i, "]]"),
    #   .message = "Please supply a list of weights output from the createWeights function."
    # )
  })
}





# Currently unused
check_interactions <- function(tv_conf, ti_conf, regex_time_pts = "([0-9]+)$") {
  idx_tv_interactions = which(grepl("\\:", tv_conf))
  idx_ti_interactions = which(grepl("\\:", ti_conf))
  sapply(
    strsplit(tv_conf[idx_tv_interactions], "\\:"), 
    function(x) {
      in_tv_conf = x %in% tv_conf
      in_ti_conf = x %in% ti_conf

      if (!any(in_tv_conf)) {
        stop(
          paste0("The term `", paste0(x, collapse = ":"), "` is invalid. At least one variable in a time-varying confounder interaction must be a time-varying confounder."), 
          call. = FALSE
        )
      }
      if (!all(in_tv_conf | in_ti_conf)) {
        stop(
          paste0("The term `", paste0(x, collapse = ":"), "` is invalid. At least one of the terms is not present by itself in `tv_conf` or `ti_conf`."), 
          call. = FALSE
        )
      }
    }
  )

  sapply(
    strsplit(ti_conf[idx_ti_interactions], "\\:"), function(x) {
      in_tv_conf = x %in% tv_conf
      in_ti_conf = x %in% ti_conf

      if (any(in_tv_conf)) {
        stop(
          paste0("The term `", paste0(x, collapse = ":"), "` is invalid. It contains a time-varying confounder and should be put in `tv_conf`."), 
          call. = FALSE
        )
      }
      if (!all(in_ti_conf)) {
        stop(
          paste0("The term `", paste0(x, collapse = ":"), "` is invalid. All terms must be found in `ti_conf`."), 
          call. = FALSE
        )
      }
    }
  )
}
process_interactions <- function(vars, regex_time_pts = "([0-9]+)$") {
  sapply(strsplit(vars, "\\:"), function(var) {
      # makes tv confounder last always
      # reorder terms: ti confounder first, then in order of exposure period
      var_time_pts = .extract_time_pts_from_vars(var, regex_time_pts)
      var = c(
        var[is.na(var_time_pts)], var[!is.na(var_time_pts)][order(var_time_pts)]
      )
      var = na.omit(var)
      paste(var, collapse = ":")
    }
  )
}
