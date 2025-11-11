#' Fit outcome model
#'
#' Fits weighted marginal outcome model as a generalized linear model of the
#' user's choosing, relating exposure main effects to outcome using IPTW
#' weights.
#' @seealso [WeightIt::glm_weightit()] for more on family/link specifications.
#'
#' @param outcome name of outcome variable with ".timepoint" suffix. 
#'   See [initMSM()] for details on suffix
#' @param model character indicating one of the following outcome models:
#'  * "m0" (exposure main effects)
#'  * "m1" (exposure main effects & covariates)
#'  * "m2" (exposure main effects & their interactions)
#'  * "m3" (exposure main effects, their interactions, & covariates)
#' @param int_order integer specification of highest order exposure main effects
#'   interaction, required for interaction models ("m2", "m3")
#' @param covariates list of characters reflecting variable names of covariates,
#'   required for covariate models ("m1", "m3")
#' @param family (optional) family function specification for [WeightIt::glm_weightit()] model.
#'   Note that this should be specified as as a function, not a character string unless you have
#'   a multinomial outcome, in which case set this to "multinomial", or if you have an ordinal 
#'   outcome with more than 2 levels, in which case set this to "ordinal". (default is gaussian)
#' @param link (optional) link function specification for [WeightIt::glm_weightit()] model.
#'
#' @return list containing [WeightIt::glm_weightit()] model output. It is the length 
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
#' f <- createFormulas(obj, type = "short")
#' w <- createWeights(data = data, formulas = f)
#' 
#' fit_m0 <- fitModel(
#'   data = data, weights = w, 
#'   outcome = "D.3", model = "m0"
#' )
#' print(fit_m0)
#' 
#' fit_m1 <- fitModel(
#'   data = data, weights = w, 
#'   outcome = "D.3", model = "m1", 
#'   covariates = c("C")
#' )
#' print(fit_m1)
#' 
#' fit_m2 <- fitModel(
#'   data = data, weights = w, 
#'   outcome = "D.3", model = "m2", 
#'   int_order = 2
#' )
#' print(fit_m2)
#' 
#' fit_m3 <- fitModel(
#'   data = data, weights = w, 
#'   outcome = "D.3", model = "m3",
#'   int_order = 2, covariates = c("C")
#' )
#' print(fit_m3)
#' 
#' data <- data.frame(
#'   ID = 1:50,
#'   A.1 = rnorm(n = 50),
#'   A.2 = rnorm(n = 50),
#'   A.3 = rnorm(n = 50),
#'   B.1 = rnorm(n = 50),
#'   B.2 = rnorm(n = 50),
#'   B.3 = rnorm(n = 50),
#'   C = rnorm(n = 50),
#'   D.3 = c(rep(c("A", "B", "C"), 16), "A", "B")
#' )
#' obj <- initMSM(
#'   data,
#'   exposure = c("A.1", "A.2", "A.3"),
#'   ti_conf = c("C"),
#'   tv_conf = c("B.1", "B.2", "B.3", "D.3")
#' )
#' f <- createFormulas(obj, type = "short")
#' w <- createWeights(data = data, formulas = f)
#' 
#' fit_m0 <- fitModel(
#'   data = data, weights = w, 
#'   outcome = "D.3", model = "m0", family = "multinomial"
#' )
#' print(fit_m0)
#' 
#'
#' @export
fitModel <- function(
    data, obj, weights = NULL, outcome,
    model = "m0", int_order = NA, covariates = NULL,
    family = gaussian(), link = "identity",
    verbose = FALSE, save.out = FALSE) {
  ### Checks ----
  dreamerr::check_arg(verbose, "scalar logical")
  dreamerr::check_arg(save.out, "scalar logical | scalar character")
  
  .check_data(data)
  if (!is.null(weights)) {
    .check_weights(weights)
    obj <- attr(weights, "obj")
  } else if (is.null(weights)) {
    warning("ALERT: You are fitting an unweighted model that does not use IPTW to adjust for confounding!", call. = FALSE)
  }
  
  if (missing(obj) || !inherits(obj, "devMSM")) {
    stop("`obj` must be output from `initMSM()`.", call. = FALSE)
  }
  
  var_tab <- obj[["var_tab"]]
  exposure <- obj[["exposure"]]
  exposure_root <- obj[["exposure_root"]]
  exposure_time_pts <- obj[["exposure_time_pts"]]
  epoch <- obj[["epoch"]]
  sep <- obj[["sep"]]
  
  dreamerr::check_arg(outcome, "character scalar MBT")
  if (is.data.frame(data)) { 
    dreamerr::check_value(
      reformulate(outcome), "formula var(data)",
      .data = data, .arg_name = "outcome"
    )
  } else {
    dreamerr::check_value(
      reformulate(outcome), "formula var(data)",
      .data = data[[1]], .arg_name = "outcome"
    )
  }
  outcome_time_pt <- .extract_time_pts_from_vars(outcome, sep = sep)
  if (is.na(outcome_time_pt)) {
    stop("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your data", call. = FALSE)
  }
  if (outcome_time_pt < max(exposure_time_pts)) {
    stop("Please supply an outcome variable and time point that is equal to or greater than the last exposure time point.", call. = FALSE)
  }
  
  model <- match.arg(model, c("m0", "m1", "m2", "m3"))
  if (model %in% c("m1", "m3")) {
    if (is.null(covariates)) {
      stop("Please provide a list of covariates as characters if you select a covariate model.", call. = FALSE)
    }
    dreamerr::check_arg(covariates, "vector character")
  }
  else if (!is.null(covariates)) {
    warning(sprintf("Covariates are ignored when `model = \"%s\"`.", model),
            call. = FALSE, immediate. = TRUE)
    covariates <- NULL
  }
  
  if (model %in% c("m2", "m3")) {
    if (is.na(int_order)) {
      stop("Please provide an integer interaction order if you select a model with interactions.", call. = FALSE)
    }
    
    dreamerr::check_arg(int_order, "scalar integer")
    
    if (int_order > length(epoch)) {
      stop("Please provide an interaction order equal to or less than the total number of epoch/time points.", call. = FALSE)
    }
  }
  else if (!all(is.na(int_order))) {
    warning(sprintf("Interaction order is ignored when `model = \"%s\"`.", model),
            call. = FALSE, immediate. = TRUE)
    int_order <- NA
  }
  
  ## Check that covariates are not post-exposure and are in var_tab
  covars_time <- .extract_time_pts_from_vars(covariates, sep = obj[["sep"]])
  post_baseline_time <- covars_time > min(exposure_time_pts)
  if (any(!is.na(covars_time)) && any(post_baseline_time, na.rm = TRUE)) {
    warning(sprintf(
      "Covariates should not contain time-varying confounders measured after the first treatment time point.\nProblematic covariates: %s",
      paste(dQuote(na.omit(covariates[post_baseline_time]), FALSE), collapse = ", ")
    ), call. = FALSE)
  }
  
  # getting null comparisons for LHT
  null_model <- if (model %in% c("m0", "m2")) "int" else "covs"
  
  ### START ----
  if (inherits(data, "mids")) {
    data_type <- "mids"
    m <- data$m
  } else if (is.data.frame(data)) {
    data_type <- "data.frame"
    data <- list(data)
    m <- 1
  } else {
    data_type <- "list"
    m <- length(data)
  }
  
  # fit_glm ----
  all_fits <- all_fits_null <- vector("list", m)
  
  for (k in seq_len(m)) {
    d <- {
      if (data_type != "mids") data[[k]]
      else as.data.frame(mice::complete(data, k))
    }
    
    # `weightitMSM` object
    w <- weights[[k]]
    
    all_fits[[k]] <- .fit_glm(
      data = d, weights = w,
      outcome = outcome, exposure = exposure, epoch = epoch,
      model = model, int_order = int_order, covariates = covariates,
      family = family, link = link, sep = sep
    )
    
    all_fits_null[[k]] <- .fit_glm(
      data = d, weights = w,
      outcome = outcome, exposure = exposure, epoch = epoch,
      model = null_model, int_order = int_order, covariates = covariates,
      family = family, link = link, sep = sep
    )
  }
  
  # Test of null ----
  if (data_type %in% c("mids", "list")) {
    rlang::check_installed("mice")
    null_test <- mice::D2(all_fits, all_fits_null)
  } else {
    null_test <- anova(all_fits[[1]], all_fits_null[[1]])
  }
  
  class(all_fits) <- "devMSM_models"
  attr(all_fits, "outcome") <- outcome
  attr(all_fits, "model") <- model
  attr(all_fits, "null_test") <- null_test
  attr(all_fits, "obj") <- obj
  
  if (verbose) print(all_fits, i = 1)
  
  if (isTRUE(save.out) || is.character(save.out)) {
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "models"))
    .create_dir_if_needed(out_dir)
    
    if (is.character(save.out)) {
      file_name <- save.out
    } else {
      file_name <- sprintf(
        "outcome_%s-exposure_%s-model_%s.rds",
        gsub("\\.", "_", outcome), 
        exposure_root, 
        model
      )
    }
    
    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving model fits to `.rds` file. To load, call:\nreadRDS("%s")\n',
      out
    ))
    saveRDS(all_fits, out)
  }
  
  return(all_fits)
}

#' @rdname fitModel
#'
#' @param x devMSM_models object from `fitModel`
#' @inheritParams devMSM_common_docs
#' 
#' @param ... ignored
#' @export
print.devMSM_models <- function(x, i = NA, save.out = FALSE, ...) {
  model <- attr(x, "model")
  outcome <- attr(x, "outcome")
  null_test <- attr(x, "null_test")
  obj <- attr(x, "obj")
  data_type <- obj[["data_type"]]
  exposure_root <- obj[["exposure_root"]]
  model_class <- class(x[[1]]) # added by IS 11/11/25
  
  if (data_type %in% c("mids", "list")) data_type <- "imputed"
  else i <- 1L
  
  if (isTRUE(i)) {
    i <- seq_along(x)
  }
  
  if (length(i) > 1) {
    for (j in i) {
      print(x, i = j, save.out = save.out, ...)
    }
    return(invisible(x))
  }
  
  msg_null_test <- c(
    "Please inspect the Wald test to determine if the exposures collectively predict significant variation in the outcome compared to a model without exposure terms.",
    "\n",
    "We strongly suggest only conducting history comparisons if the test is significant.",
    "\n"
  )
  
  if (is.na(i)) {
    rlang::check_installed("mice")
    xi <- mice::pool(x, dfcom = Inf)
  }
  else {
    xi <- x[[i]]
  }
  
  # added by IS to accommodate multiple outcome levels 11/11/25
  if(any(grepl("multinom_weightit", model_class))){
    t <- modelsummary::modelsummary(
      xi,
      statistic = c("CI" = "[{conf.low}, {conf.high}]", "p" = "{p.value}"),
      shape = term + response ~ model + statistic,
      gof_map = c("nobs"),
      output = "tinytable"
    )
  } else{
    
    t <- modelsummary::modelsummary(
      xi,
      statistic = c("CI" = "[{conf.low}, {conf.high}]", "p" = "{p.value}"),
      shape = term ~ model + statistic,
      gof_map = c("nobs"),
      output = "tinytable"
    )
  }
  
  if (data_type == "imputed") {
    if (is.na(i)) {
      msg_model_sum <- sprintf(
        "\n\nThe marginal model, %s, pooled across imputed datasets is summarized below:\n",
        model
      )
    }
    else {
      msg_model_sum <- sprintf(
        "\n\nThe marginal model, %s, for imputed dataset %s is summarized below:\n",
        model, i
      )
    }
    
  } else {
    msg_model_sum <- sprintf(
      "\n\nThe marginal model, %s, is summarized below:\n",
      model
    )
  }
  message(msg_null_test)
  print(null_test)
  cat(msg_model_sum)
  print(t, "markdown")
  
  if (isTRUE(save.out) || is.character(save.out)) {
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "models"))
    .create_dir_if_needed(out_dir)
    
    if (data_type == "imputed") {
      if (is.na(i)){
        i_str <- sprintf("-all_avg")
      } else{
        i_str <- sprintf("-imp_%s", i)}
    } else {
      i_str <- ""
    }
    
    if (is.character(save.out)) {
      file_name <- save.out
    } else {
      file_name <- sprintf(
        "fit_model_summary-outcome_%s-exposure_%s-model_%s%s.docx",
        gsub("\\.", "_", outcome), 
        exposure_root, 
        model, i_str
      )
    }
    
    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving model statistics table to file:\n%s\n',
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

# Internal call
.fit_glm <- function(
    data, weights = NULL, outcome, exposure, epoch,
    model, int_order, family, link, covariates, sep) {
  
  # Create epoch-averages of exposure variables
  epoch_vars <- exposure
  if (any(exposure != epoch)) {
    epoch_df <- .generate_epoch_exposures(data, exposure, epoch, sep)
    epoch_vars <- colnames(epoch_df)
    missing_columns <- !(colnames(epoch_df) %in% colnames(data))
    data <- cbind(data, epoch_df[, missing_columns])
  }
  
  if (model == "m2") {
    interactions <- .get_interaction_terms(epoch_vars, int_order)
  } else if (model == "m3") { # added to interact exp w/ covars
    interactions <- lapply(covariates, function(x) {
      .get_interaction_terms(c(epoch_vars, x), int_order)
    })
    interactions <- unique(unlist(interactions))
  }
  
  covs <- switch(model,
                 "int" = "1",
                 "covs" = covariates,
                 "m0" = epoch_vars,
                 "m1" = c(epoch_vars, covariates),
                 "m2" = c(epoch_vars, interactions),
                 "m3" = c(epoch_vars, interactions)) #, covariates))
  
  # IS added 10/29/25
  if(!class(family) %in% ("function")){
    if(any(family %in% "multinomial")){ # multinomial outcome
      WeightIt::multinom_weightit(
        formula = reformulate(covs, response = outcome),
        data = data,
        link = link,
        weightit = weights
      )
    }else if(any(family %in% "ordinal")){ # ordinal outcome with more than 2 levels
      WeightIt::ordinal_weightit(
        formula = reformulate(covs, response = outcome),
        data = data,
        link = link,
        weightit = weights
      )
    }else{
      WeightIt::glm_weightit(
        formula = reformulate(covs, response = outcome),
        data = data,
        family = family,
        link = link,
        weightit = weights
      )
    }} else{ # send function to glm_weightit
      WeightIt::glm_weightit(
        formula = reformulate(covs, response = outcome),
        data = data,
        family = family,
        link = link,
        weightit = weights)
    }
}