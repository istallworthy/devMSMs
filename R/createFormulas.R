#' Create balancing formulas
#'
#' Creates balancing formulas relating exposure to all relevant time-varying and
#' time invariant confounders at each exposure time point to be used to create
#' IPTW weights.
#'
#' @param obj initialized MSM object from `initMSM`
#' @param type type of formula to create from 'full' (includes all lagged
#'   time-varying confounders), 'short' (includes time-varying confounders at
#'   t-1 lag only), or 'update' (adds to 'short' formulas any imbalanced
#'   time-varying confounders at lags great than t-1)
#' @param concur_conf (optional) list of variable names reflecting time-varying
#'   confounders to retain in formulas contemporaneously (default is none)
#' @param keep_conf (optional) list of variable names reflecting confounders to
#'   always retain in formulas (default depends on type)
#' @param bal_stats list of balance statistics from assessBalance(), required
#'   for 'update' type
#' @param custom (optional) custom list of formulas at each exposure time point
#'   (default is to create automatically according to type)
#' @param verbose (optional) TRUE or FALSE indicator for user output
#'   (default is TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is FALSE)
#' @param home_dir path to home directory (required if save.out = TRUE)
#'
#' @return list of balancing formulas at each exposure time point
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
#' # Full Formulas
#' f <- createFormulas(obj, type = "full")
#'
#' # Short Formulas
#' f <- createFormulas(obj, type = "short")
#'
#' # Update Formulas
#' w <- createWeights(data = data, obj = obj, formulas = f)
#' b <- assessBalance(data = data, obj = obj, weights = w)
#' f <- createFormulas(obj, type = "update", bal_stats = b)
#'
#' @export
createFormulas <- function(
    obj, type = c("full", "short", "update"), custom = NULL,
    concur_conf = NULL, keep_conf = NULL, bal_stats = NULL,
    verbose = FALSE, save.out = FALSE, home_dir = NULL) {
  
  ### Checks ----
  dreamerr::check_arg(verbose, save.out, "scalar logical")
  if (save.out) dreamerr::check_arg_plus(home_dir, "path dir")

  # Extract var_tab and get variables
  var_tab <- attr(obj, "var_tab")
  ti_conf <- var_tab$var[var_tab$type == "ti_conf"]
  ti_conf_time <- var_tab$time[var_tab$type == "ti_conf"]
  tv_conf <- var_tab$var[var_tab$type == "tv_conf"]
  tv_conf_time <- var_tab$time[var_tab$type == "tv_conf"]
  exposure <- var_tab$var[var_tab$type == "exposure"]
  exposure_time <- var_tab$time[var_tab$type == "exposure"]

  # Check concur_conf and keep_conf
  dreamerr::check_arg(concur_conf, keep_conf, "vector character | NULL")
  if (!all(concur_conf %in% tv_conf)) {
    stop("The variables in the concur_conf field must be included in the tv_conf.", call. = FALSE)
  }
  if (any(exposure %in% concur_conf)) {
    stop("Do not include the exposure concurrently as a confounder. Please revise the concur_conf field.", call. = FALSE)
  }
  if (!all(keep_conf %in% tv_conf)) {
    stop("The variables in the keep_conf field must be included in the tv_conf.", call. = FALSE)
  }
  if (any(exposure %in% keep_conf)) {
    stop("Do not include the exposure concurrently as a confounder. Please revise the concur_conf field.", call. = FALSE)
  }

  type <- match.arg(type, several.ok = FALSE)
  if (type != "update" && !is.null(bal_stats)) {
    stop("Please only provide balance statistics for the type 'update'.", call. = FALSE)
  }
  if (type == "update") {
    # TODO: Check the class
    dreamerr::check_arg(
      bal_stats, "class(devMSM_bal_stats)",
      .message = "Please provide a data frame of balance statistics from the assessBalance function."
    )
  }

  # TODO: custom formulas passed in
  ### For custom formulas, check and then pass through ----
  # if (!is.null(custom_fmla)) {
  #   check_arg(
  #     custom_fmla, "list len(data)",
  #     .data = exposure,
  #     .message = "If you wish to supply custom formulas, please provide a list of formulas for each exposure time point."
  #   )
  #
  #   # Very minorly clean custom formulas
  #   # TODO: Check that data check works in formula check
  #   forms <- lapply(seq_along(custom_fmla), function(i) {
  #     form <- custom_fmla[[i]]
  #
  #     dreamerr::check_value(
  #       form, "ts formula var(data)",
  #       .data = data,
  #       .arg_name = paste0("custom[[", i, "]]")
  #     )
  #
  #     # NOTE: This works because we have checked for `ts formula` above
  #     form_outcome <- deparse1(form[[2]])
  #     outcome_time_pt <- .extract_time_from_vars(form_outcome)
  #     if (!(form_outcome %in% exposure)) {
  #       stop(sprintf(
  #         "%s\nThe LHS of the formula is not one of the exposure var: %s.",
  #         deparse1(form), paste(exposure, sep = ",")
  #       ), call. = FALSE)
  #     }
  #
  #     form_rhs <- strsplit(deparse1(form[[3]]), "\\s+\\+\\s+")[[1]]
  #     # Just-in-case, still sort based on time-varying covariates in interactions
  #     form_rhs <- process_interactions(form_rhs)
  #     covars_time <- stats::na.omit(.extract_time_from_vars(form_rhs))
  #
  #     if (length(covars_time) > 0 && any(covars_time > outcome_time_pt)) {
  #       stop(sprintf(
  #         "%s\nBalancing formula cannot contain time-varying confounders measured after the exposure time point for that formula.",
  #         deparse1(form)
  #       ), call. = FALSE)
  #     }
  #
  #     form_rhs_vars <- all.vars(reformulate(form_rhs))
  #     miss <- !(form_rhs_vars %in% c(tv_conf, ti_conf))
  #     if (any(miss)) {
  #       warning(sprintf(
  #         "Please make sure all variables in your custom formulas are included as either time-varying or time invariant confounders.\nThe following variables are not: %s", paste(form_rhs_vars[miss], sep = ", ")
  #       ), call. = FALSE)
  #     }
  #
  #     f <- reformulate(form_rhs, response = form_outcome)
  #     class(f) <- c("devMSM_formulas", "formula")
  #     attr(f, "type") <- type
  #     attr(f, "time") <- exposure_time[i]
  #     attr(f, "exposure") <- exposure
  #     attr(f, "outcome") <- outcome
  #     attr(f, "formula_msg") <- sprintf(
  #       "At time point %s, the %s formula for %s and outcome %s is: \n%s",
  #       exposure_time[i], type, exposure[i], outcome, deparse1(f)
  #     )
  #     return(f)
  #   })
  #
  #   # TODO: Figure out this error
  #   # if (any(
  #   #   as.logical(unlist(lapply(1:length(covars), function(x) {
  #   #     any(na.omit(as.numeric(sapply(strsplit(unlist(covars[x]), "\\."), "[", 2))) ==
  #   #     as.numeric(sapply(strsplit(names(covars[x]), "-"), "[", 2)))
  #   #   })))
  #   # )) {
  #   #   warning(sprintf("Be sure that time-varying confounders included in balancing formulas contemporaneously are not mediators."),#     call. = FALSE
  #   #   )
  #   # }
  #
  #   names(forms) <- paste(paste0(type, "_form"), exposure_time, sep = "-")
  #   if (verbose) {
  #     message("The user-supplied custom balancing formula for each exposure time point are below: ")
  #     lapply(forms, print)
  #   }
  #
  #   return(forms)
  # }

  ### Create formulas  ----
  forms <- list()
  exposure_times <- c()
  exposure_vars <- c()
  msgs <- c()
  update_msgs <- c()
  for (i in seq_along(exposure_time)) {
    i_time <- exposure_time[i]
    im1_time <- if (i == 1) 0 else exposure_time[i - 1]

    ti_include <- ti_conf
    if (type == "full") {
      tv_include <- tv_conf[tv_conf_time < i_time]
      tv_include <- c(tv_include, exposure[exposure_time < i_time])
    } else if (type == "short" | type == "update") {
      tv_include <- tv_conf[tv_conf_time >= im1_time & tv_conf_time < i_time]
      tv_include <- c(tv_include, exposure[exposure_time >= im1_time & exposure_time < i_time])
    }

    if (type == "update") {
      update_include <- bal_stats$covariate[
        bal_stats$balanced == 0 &
          bal_stats$exp_time == i_time &
          as.numeric(bal_stats$covar_time) < im1_time &
          as.numeric(bal_stats$covar_time) > 0
      ]
      tv_include <- c(tv_include, update_include)

      # if (length(update_include) > 0) {
      #   update_msg <- sprintf(
      #     "For %s at exposure time point %s the following covariate(s) will be added to the short balancing formula: %s \n\n",
      #     exposure, time, paste(new, collapse = ", ")
      #   )
      # } else {
      #   update_msg <- sprintf(
      #     "For %s at exposure time point %s no time-varying confounders at additional lags were added. \n\n",
      #     exposure, time
      #   )
      # }
    }

    # adding in any user-specified concurrent confounders (default is only lagged)
    concur_time <- .extract_time_pts_from_vars(concur_conf)
    concur_include <- concur_conf[concur_time == i_time]
    keep_time <- .extract_time_pts_from_vars(keep_conf)
    keep_include <- keep_conf[keep_time == i_time]

    vars_to_include <- unique(c(ti_include, tv_include, keep_include, concur_include))

    # Creates formula for the given exposure time point
    f <- reformulate(vars_to_include, response = exposure[i])

    if (!exists("update_msg")) update_msg <- NULL
    exposure_times <- c(exposure_times, i_time)
    exposure_vars <- c(exposure_vars, exposure[i])
    msgs <- c(
      msgs,
      sprintf(
        "At time point %s, the %s formula for %s is: \n%s",
        i_time, type, exposure[i], deparse1(f)
      )
    )
    update_msgs <- c(update_msgs, update_msg)

    # Assigns the form to forms list
    forms[[i]] <- f
  }

  ### Output ----
  class(forms) <- c("devMSM_formulas", "list")
  attr(forms, "type") <- type
  attr(forms, "exposure_times") <- exposure_times
  attr(forms, "exposure_vars") <- exposure_vars
  attr(forms, "msgs") <- msgs
  attr(forms, "update_msgs") <- update_msgs

  if (verbose) print(forms, short = FALSE)
  if (save.out) {
    forms_dir <- file.path(home_dir, "formulas", type)
    if (!dir.exists(forms_dir)) dir.create(forms_dir, recursive = TRUE)
    msgs <- sapply(forms, function(f) attr(f, "formula_msg"))
    forms_csv_file <- file.path(forms_dir, sprintf(
      "%s_%s-%s_%s_balancing_formulas.csv",
      type, exposure, type
    ))
    writeLines(msgs, con = forms_csv_file)
  }

  return(forms)
}


#' @export
print.devMSM_formulas <- function(x, ...) {
  user_alert <- switch(attr(x, "type"),
    "full" = "USER ALERT: Please manually inspect the full balancing formula below:",
    "short" = "USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:",
    "update" = "USER ALERT: Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:",
    "custom" = ""
  )
  
  lapply(seq_along(x), function(i) {
    msg <- attr(x, "msgs")[[i]]
    message(user_alert)

    if (i < length(x)) msg <- paste0(msg, "\n\n")
    cat(msg)

    return(NULL)
  })
}
