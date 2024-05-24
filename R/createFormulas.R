#' Create balancing formulas
#'
#' Creates balancing formulas relating exposure to all relevant time-varying and
#' time invariant confounders at each exposure time point to be used to create
#' IPTW weights.
#'
#' @inheritParams devMSM_common_docs
#' @param type type of formula to create from 'full' (includes all lagged
#'   time-varying confounders), 'short' (includes time-varying confounders at
#'   t-1 lag only), or 'update' (adds to 'short' formulas any imbalanced
#'   time-varying confounders at lags great than t-1)
#' @param keep_conf (optional) list of variable names reflecting confounders to
#'   always retain in formulas (default depends on type)
#' @param custom (optional) custom list of formulas at each exposure time point
#'   (default is to create automatically according to type)
#'
#' @return a list containing balancing formulas. It is the length of the number of exposure variables.
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
#' print(f)
#'
#' # Short Formulas
#' f <- createFormulas(obj, type = "short")
#' print(f)
#'
#' # Update Formulas
#' w <- createWeights(data = data, obj = obj, formulas = f)
#' b <- assessBalance(data = data, obj = obj, weights = w)
#' f <- createFormulas(obj, type = "update", bal_stats = b)
#' print(f)
#'
#' @export
createFormulas <- function(
    obj, type = c("full", "short", "update"), custom = NULL,
    keep_conf = NULL, bal_stats = NULL,
    verbose = FALSE, save.out = FALSE) {
  
  ### Checks ----
  dreamerr::check_arg(verbose, save.out, "scalar logical")
  
  home_dir <- attr(obj, "home_dir")
  if (is.null(home_dir) && save.out){
    stop("Please provide a home directory in the MSM object to save.", call. = FALSE)
  } else if (save.out) dreamerr::check_arg_plus(home_dir, "path dir")

  # Extract var_tab and get variables
  var_tab <- attr(obj, "var_tab")
  ti_conf <- var_tab$var[var_tab$type == "ti_conf"]
  ti_conf_time <- var_tab$time[var_tab$type == "ti_conf"]
  tv_conf <- var_tab$var[var_tab$type == "tv_conf"]
  tv_conf_time <- var_tab$time[var_tab$type == "tv_conf"]
  if (!is.null(var_tab$concur_conf)) { # added to fix bug if no concur_conf
    tv_conf_time <- tv_conf_time - 0.01 * var_tab$concur_conf[var_tab$type == "tv_conf"]
  }
  exposure <- var_tab$var[var_tab$type == "exposure"]
  exposure_time <- var_tab$time[var_tab$type == "exposure"]
  
  data_type <- attr(bal_stats, "data_type")
  

  # Check keep_conf
  dreamerr::check_arg(keep_conf, "vector character | NULL")
  if (!all(keep_conf %in% tv_conf)) {
    stop("The variables in the keep_conf field must be included in the tv_conf.", call. = FALSE)
  }
  if (any(exposure %in% keep_conf)) {
    stop("Do not include the exposure concurrently as a confounder. Please revise the keep_conf field.", call. = FALSE)
  }

  type <- match.arg(type, several.ok = FALSE)
  if (type != "update" && !is.null(bal_stats)) {
    stop("Please only provide balance statistics for the type 'update'.", call. = FALSE)
  }
  if (type == "update" & is.null(bal_stats)) {
    stop("You set the `type` argument to `'update'`, but did not provide `bal_stats`.", call. = FALSE)
  }
  if (type == "update") {
    dreamerr::check_arg(
      bal_stats, "class(devMSM_bal_stats)",
      .message = "Please provide a data frame of balance statistics from the assessBalance function."
    )
  }

  ### For custom formulas, check and then pass through ----
  if (!is.null(custom)) {
    dreamerr::check_arg(
      custom, "list len(data)",
      .data = exposure,
      .message = "If you wish to supply custom formulas, please provide a list of formulas for each exposure time point."
    )

    # Very minorly clean custom formulas
    processed <- lapply(seq_along(custom), function(i) {
      
      form <- custom[[i]]
      dreamerr::check_value(
        form, "ts formula",
        .arg_name = paste0("custom[[", i, "]]")
      )

      ## Check exposure variable (LHS)
      exposure_var <- rlang::expr_text(rlang::f_lhs(form))
      if (!(exposure_var == exposure[i])) {
        stop(sprintf(
          "%s\nThe LHS of the formula is not the correct exposure var: %s.",
          deparse1(form), paste(exposure[[i]], sep = ",")
        ), call. = FALSE)
      }
      
      ## Check covariates
      included_covs <- all.vars(form)
      covars_time <- .extract_time_pts_from_vars(included_covs, sep = attr(obj, "sep"))
      if (any(!is.na(covars_time)) && any(covars_time > exposure_time[i], na.rm = TRUE)) {
        # TODO: warn or stop?
        stop(sprintf(
          "%s\nBalancing formula cannot contain time-varying confounders measured after the exposure time point for that formula.",
          deparse1(form)
        ), call. = FALSE)
      }
      
      miss <- !(included_covs %in% c(tv_conf, ti_conf, exposure))
      if (any(miss)) {
        warning(sprintf(
          "Please make sure all variables in your custom formulas are included as either exposure variables, time-varying confounders, or time invariant confounders.\nThe following variables are not: %s", paste(included_covs[miss], sep = ", ")
        ), call. = FALSE)
      }

      # TODO: Process interactions? Seems not needed

      return(list(
        formula = form,
        exposure_time = exposure_time[i],
        exposure_var = exposure_var,
        msg = sprintf(
          "At time point %s, the custom formula for %s is: \n%s",
          exposure_time[i], exposure[i], deparse1(form)
        )
      ))
    })

    forms <- lapply(processed, function(x) x$formula)
    class(forms) <- c("devMSM_formulas", "list")
    attr(forms, "type") <- "custom"
    attr(forms, "exposure_times") <- sapply(processed, function(x) x$exposure_time)
    attr(forms, "exposure_vars") <- sapply(processed, function(x) x$exposure_var)
    attr(forms, "msgs") <- sapply(processed, function(x) x$msg)

    if (verbose) print(forms)
    return(forms)
  }

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
      # TODO: Check this <--IS: this need to be looking at balance stats averaged across imputed data, w/ imputed data; not working w/ imputed data
      ## IS added
      if (data_type == "data.frame") {  # non-imputed data
        bal_stats <- bal_stats[[1]]
        bal_stats <- do.call(rbind, lapply(bal_stats, as.data.frame))
      } else if (data_type == "list" || data_type == "mids") { # imputed data --needs averaging
        bal_stats <- .avg_imp_bal_stats(bal_stats)
      }

      update_include <- bal_stats$covariate[
         bal_stats$balanced == 0 &
          bal_stats$exposure_time == i_time &
          as.numeric(bal_stats$covar_time) < im1_time &
          as.numeric(bal_stats$covar_time) > 0
      ]
      tv_include <- c(tv_include, update_include)
      
    }

    # adding in any user-specified concurrent confounders (default is only lagged)
    keep_time <- .extract_time_pts_from_vars(keep_conf, sep = attr(obj, "sep"))
    keep_include <- keep_conf[keep_time < i_time]

    vars_to_include <- unique(c(ti_include, tv_include, keep_include))

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

  if (verbose) print(forms)
  return(forms)
}

#' @rdname createFormulas
#'
#' @param x devMSM_formulas object from [createFormulas()]
#' @param ... ignored
#' @export
print.devMSM_formulas <- function(x, ...) {
  user_alert <- switch(attr(x, "type"),
    "full" = "USER ALERT: Please manually inspect the full balancing formula below:",
    "short" = "USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:",
    "update" = "USER ALERT: Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:",
    "custom" = "USER ALERT: Please manually inspect the slightly cleaned custom formula below: "
  )

  lapply(seq_along(x), function(i) {
    msg <- attr(x, "msgs")[[i]]
    message(user_alert)

    if (i < length(x)) msg <- paste0(msg, "\n\n")
    cat(msg)

    return(NULL)
  })
  return(invisible(NULL))
}
