#' Estimate, compare, and visualize exposure histories
#'
#' Takes fitted model output to created predicted values for user-specified
#' histories (pooling for imputed data), before conducting contrast comparisons
#' (pooling for imputed data), correcting for multiple comparisons, and then
#' plotting results.
#' @seealso
#'  [marginaleffects::avg_predictions()],
#'  <https://cran.r-project.org/web/packages/marginaleffects/marginaleffects.pdf>;
#'  [marginaleffects::hypotheses()],
#'  <https://cran.r-project.org/web/packages/marginaleffects/marginaleffects.pdf>;
#'  [stats::p.adjust()],
#'  <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust>;
#'
#' @inheritParams devMSM_common_docs
#' @param hi_lo_cut list of two numbers indicating quantile values
#'   that reflect high and low values, respectively, for continuous exposure
#' @param dose_level (optional) "l" or "h" indicating whether low or high doses
#'   should be tallied in tables and plots (default is high "h")
#' @param reference lists of one or more strings of "-"-separated "l"
#'   and "h" values indicative of a reference exposure history to which to
#'   compare comparison, required if comparison is supplied
#' @param comparison (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is supplied
#' @param mc_comp_method (optional) character abbreviation for multiple
#'   comparison correction method for stats::p.adjust, default is
#'   Benjamini-Hochburg ("BH")
#'
#' @return list containing two dataframes: `preds` with predictions from
#'  [marginaleffects::avg_predictions()] containing average expected outcome
#'  for different exposure histories and `comps` with contrasts from
#'  [marginaleffects::comparisons()] comparing different exposure history
#'
#' @examples
#' library(devMSMs)
#' set.seed(123)
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
#' fit <- fitModel(
#'   data = data, weights = w,
#'   outcome = "D.3", model = "m0"
#' )
#' 
#' comp <- compareHistories(
#'   fit = fit,
#'   hi_lo_cut = c(0.3, 0.6)
#' )
#' print(comp)
#' plot(comp)
#' summary(comp, "preds")
#' summary(comp, "comps")
#' 
#' comp2 <- compareHistories(
#'   fit = fit,
#'   hi_lo_cut = c(0.3, 0.6),
#'   reference = "l-l-l",
#'   comparison = c("h-h-h", "h-h-l")
#' )
#' print(comp2)
#' plot(comp2)
#' summary(comp2, "preds")
#' summary(comp2, "comps")
#'
#' @export
compareHistories <- function(
    fit,
    hi_lo_cut,
    dose_level = "h",
    reference = NULL, comparison = NULL,
    mc_comp_method = "BH",
    verbose = FALSE, save.out = FALSE) {
  ### Checks ----
  dreamerr::check_arg(verbose, "scalar logical")

	dreamerr::check_arg(save.out, "scalar logical | scalar character")
  if (verbose) {
    rlang::check_installed("tinytable")
  }

  dreamerr::check_arg(fit, "class(devMSM_models)")

  # Get objects from `obj`
  obj <- attr(fit, "obj")
  exposure <- obj[["exposure"]]
  exposure_time_pts <- obj[["exposure_time_pts"]]
  exposure_type <- obj[["exposure_type"]]
  exposure_root <- obj[["exposure_root"]]
  epoch <- obj[["epoch"]]
  sep <- obj[["sep"]]
  outcome <- attr(fit, "outcome")
  model <- attr(fit, "model")

  n_epoch <- length(unique(epoch))
  exposure_levels <- apply(
    perm2(n_epoch, c("l", "h")),
    1,
    function(z) {
      paste(z, sep = "", collapse = "-")
    }
  )

  dreamerr::check_arg(reference, "character vector | NULL")
  dreamerr::check_arg(comparison, "character vector | NULL")
  if ((is.null(comparison) && !is.null(reference)) || (!is.null(comparison) && is.null(reference))) {
    stop("If either comparison or reference is specified, then both must be specified", call. = FALSE)
  }
  
  is_invalid_comparison <- !(comparison %in% exposure_levels)
  if (any(is_invalid_comparison)) {
    stop(sprintf("The following elements of `comparison` are invalid: %s", 
                 paste0(comparison[is_invalid_comparison], collapse = ", ")), call. = FALSE)
  }
  
  is_invalid_reference <- !(reference %in% exposure_levels)
  if (any(is_invalid_reference)) {
    stop(sprintf("The following elements of `reference` are invalid: %s", 
                 paste0(reference[is_invalid_reference], collapse = ", ")), call. = FALSE)
  }

  # Multiple comparison methods
  dreamerr::check_arg(mc_comp_method, "scalar character")
  mc_comp_method <- match.arg(mc_comp_method, stats::p.adjust.methods)

  if (exposure_type == "continuous") {
    dreamerr::check_arg(hi_lo_cut, "vector numeric len(2) GE{0} LE{1}")
    hi_lo_cut <- sort(hi_lo_cut)
  } else {
    hi_lo_cut <- NULL
  }

  # Dose_level
  dose_level <- match.arg(dose_level, c("h", "l"))

  # STEP 1 ----
  # Define variables for average partial effects
  # getting data to use for determining hi/lo values: should have any epochs created / used in the model
  data <- do.call("rbind", lapply(fit, function(z) z[["data"]]))

  epoch_vars <- exposure
  if (!all(exposure == epoch)) {
    epoch_vars <- .get_epoch_var_names(exposure, epoch, sep = sep)
  }

  mat <- fit[[1]][["data"]][, epoch_vars, drop = FALSE]
  epoch_history <- .characterize_exposure(mat, exposure_type, hi_lo_cut)

  prediction_vars <- .get_avg_predictions_variables(
    data = data, epoch_vars = epoch_vars,
    exposure_type = exposure_type, hi_lo_cut = hi_lo_cut
  )
  epoch_d <- data.frame(
    e = names(prediction_vars),
    l = sapply(prediction_vars, function(x) x[1]),
    h = sapply(prediction_vars, function(x) x[2]),
    row.names = NULL
  )
  epoch_d$z <- prediction_vars


  # STEP 2 ----
  # Estimated marginal predictions for each history
  preds <- lapply(fit, function(m) { # Goes through different fitted model
    d <- m$data

    p <- marginaleffects::avg_predictions(
      m,
      newdata = d, variables = prediction_vars,
      type = "response"
    )
    
    .add_histories(p, epoch_vars)
  })

  # STEP 3 ----
  # Conduct history comparisons
  if (is.null(reference) && is.null(comparison)) {
    # All pairwise comparisons if no ref/comparisons were specified by user
    # Pairwise comparision, but we want term to be nicely labeled by histories
    comps <- lapply(preds, function(y) {
      h <- marginaleffects::hypotheses(
        model = as.data.frame(y), vcov = vcov(y), hypothesis = "pairwise"
      )
      class(h) <- c("comp_custom", class(h))
      return(h)
    })
  } else {
    hypothesis <- .create_hypotheses_mat(preds[[1]]$term, reference, comparison)
    
    # IS added to remove same-same history comparisons 
    nonid <- colnames(hypothesis)[sapply(strsplit(colnames(hypothesis), " - "), "[", 1) != 
                                    sapply(strsplit(colnames(hypothesis), " - "), "[", 2)]
    hypothesis <- hypothesis[, nonid, drop = FALSE]
    
    if (length(hypothesis) == 0L){
      stop("Please a comparison history that differs from the reference history.", call. = FALSE)
    }

    # for normal outcome variables
    comps <- lapply(preds, function(y) {
      h <- marginaleffects::hypotheses(y, hypothesis = hypothesis)
      class(h) <- c("comp_custom", class(h))
      return(h)
    })
    

  }

  # STEP 4 ----
  # pooling predicted values and contrasts for imputed data
  is_pooled <- length(preds) > 1
  if (is_pooled) { # IMPUTED DATA
    # TODO: Delete this code when `modelsummary` goes onto CRAN & require that version in DESCRIPTION 
    # `modelsummary (>= 1.4.4),`
    # TODO: should be fixed by modelsummary: https://github.com/vincentarelbundock/modelsummary/commit/5f13fe03683016ae92e5ffdd4b8b6b402614409e
    # 
    # assign(
    #   "glance_custom_internal.glm_weightit", 
    #   function(model, vcov_type = vcov_type, gof = gof) {
    #     out = NextMethod("glance", model)
    #     if ("logLik" %in% colnames(out)) {
    #       out$logLik = as.numeric(out$logLik)
    #     }
    #     return(out)
    #   },
    #   envir = globalenv()
    # )
    # on.exit({ rm(glance_custom_internal.glm_weightit, envir = globalenv()) })

    rlang::check_installed("mice")
    preds <- summary(mice::pool(preds, dfcom = Inf), conf.int = TRUE)
    comps <- summary(mice::pool(comps, dfcom = Inf), conf.int = TRUE)
    names(preds)[7:8] <- names(comps)[7:8] <- c("conf.low", "conf.high")
  } else { # REGULAR DATA
    comps <- as.data.frame(comps[[1]])
    preds <- as.data.frame(preds[[1]])
  }

  preds$term <- as.character(preds$term)
  preds$dose <- .dose_from_history(preds$term, dose_level = dose_level)
  preds <- preds[, c("term", setdiff(colnames(preds), "term"))]

  comps$term <- as.character(comps$term)
  comps$dose <- .dose_from_history_contrast(comps$term, dose_level = dose_level)
  comps$p.value_corr <- stats::p.adjust(comps$p.value, method = mc_comp_method)

  # If the user specified reference and comparison groups, subset pred_pool for inspection and plotting
  if (!is.null(reference) && !is.null(comparison)) {
    idx_keep <- preds$term %in% c(reference, comparison)
    preds <- preds[idx_keep, , drop = FALSE]
  }

  # Return ----
  exp_lab <- .remove_time_pts_from_vars(exposure[1], sep = sep)

  res <- list(preds = preds, comps = comps)
  class(res) <- "devMSM_comparisons"
  attr(res, "obj") <- obj
  attr(res, "outcome") <- outcome
  attr(res, "hi_lo_cut") <- hi_lo_cut
  attr(res, "is_pooled") <- is_pooled
  attr(res, "epoch_d") <- epoch_d
  attr(res, "epoch_history") <- epoch_history
  attr(res, "mc_comp_method") <- mc_comp_method
  attr(res, "exp_lab") <- exp_lab
  
  if (!is.null(reference) && !is.null(comparison)) {
    attr(res, "reference") <- reference
    attr(res, "comparison") <- comparison
  }

  if (verbose) print(res)

  if (isTRUE(save.out) || is.character(save.out)) {
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "histories"))
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
    saveRDS(res, out)
  }

  return(res)
}


#' @rdname compareHistories
#'
#' @param x devMSM_histories object from [compareHistories()]
#' @param ... ignored
#'
#' @export
print.devMSM_comparisons <- function(x, save.out = FALSE, ...) {
  preds <- x$preds
  comps <- x$comps
  hi_lo_cut <- attr(x, "hi_lo_cut")
  epoch_history <- attr(x, "epoch_history")
  mc_comp_method <- attr(x, "mc_comp_method")
  
  outcome <- attr(x, "outcome")
  obj <- attr(x, "obj")
  epoch <- obj[["epoch"]]
  exposure_root <- obj[["exposure_root"]]

  reference <- NULL
  comparison <- NULL
  if (!is.null(attr(x, "reference"))) {
    reference <- attr(x, "reference")
  }
  if (!is.null(attr(x, "comparison"))) {
    comparison <- attr(x, "comparison")
  }

  columns_to_drop <- c("statistic", "s.value", "df", "p.value", "dose")
  preds_tab <- tinytable::tt(preds[, setdiff(colnames(preds), columns_to_drop)])
  preds_tab <- tinytable::format_tt(preds_tab, digits = 2)

  columns_to_drop <- c("statistic", "s.value", "df", "dose")
  comps_tab <- tinytable::tt(comps[, setdiff(colnames(comps), columns_to_drop)])

  cat("Summary of Exposure Main Effects:\n")
  .print_eval_hist(epoch_history, epoch, hi_lo_cut, reference, comparison)
  cat("Below are the pooled average predictions by user-specified history:")
  print(preds_tab, "markdown")
  cat(sprintf("\nConducting multiple comparison correction for all pairings between comparison histories and each reference history using the %s method. \n", mc_comp_method))
  cat("\n")
  print(comps_tab, "markdown")

  if (isTRUE(save.out) || is.character(save.out)) {
    rlang::check_installed("pandoc")
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "histories"))
    .create_dir_if_needed(out_dir)

    if (is.character(save.out)) {
      file_name = save.out
    } else {
      file_name <- sprintf(
        "comparisons_table-outcome_%s-exposure_%s.docx", 
        gsub("\\.", "\\_", outcome), 
        exposure_root
      )
    }
    
    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving comparisons table to file:\n%s\n',
      out
    ))
    if (fs::path_ext(out) == "pdf") {
      tinytable::save_tt(
        tinytable::format_tt(comps_tab, escape = TRUE),
        output = out, overwrite = TRUE
      )
    } else {
      tinytable::save_tt(comps_tab, output = out, overwrite = TRUE)
    }
  }

  return(invisible(comps_tab))
}


#' @rdname compareHistories
#'
#' @inheritParams devMSM_common_docs
#' @param x devMSM_histories object from [compareHistories()]
#' @param colors (optional) character specifying Brewer palette or list of
#'   colors (n(epochs)+1) for plotting (default is "Dark2" palette)
#' @param exp_lab (optional) character label for exposure variable in plots
#'   (default is variable name)
#' @param out_lab (optional) character label for outcome variable in plots
#'   (default is variable name)
#' @param ... ignored
#'
#' @export
plot.devMSM_comparisons <- function(x, colors = "Dark2", exp_lab = NULL, out_lab = NULL, save.out = FALSE, ...) {
  dreamerr::check_arg(exp_lab, out_lab, "NULL | character scalar")
  if (is.null(out_lab)) out_lab <- attr(x, "outcome")
  if (is.null(exp_lab)) exp_lab <- attr(x, "exp_lab")

  outcome <- attr(x, "outcome")
  obj <- attr(x, "obj")
  exposure <- obj[["exposure"]]
  exposure_root <- obj[["exposure_root"]]
  epoch <- obj[["epoch"]]

  n_epoch <- length(unique(epoch))
  dreamerr::check_arg(colors, "character vector")
  if (length(colors) > 1 && length(colors) != n_epoch + 1) {
    stop(
      sprintf(
        "Please provide either %s different colors, a Brewer color palette, or leave this entry blank. \n",
        n_epoch + 1
      ),
      call. = FALSE
    )
  }

  preds <- x$preds
  preds$history <- preds$term
  preds$low_ci <- preds$estimate + (qnorm(.025) * preds$std.error)
  preds$high_ci <- preds$estimate + (qnorm(.975) * preds$std.error)
  preds <- preds[order(preds$dose), , drop = FALSE]
  preds$dose <- as.factor(preds$dose)

  if (length(colors) > 1) { # If user input a list of colors
    color_scale <- ggplot2::scale_color_manual(values = colors)
  } else {
    color_scale <- ggplot2::scale_colour_brewer(palette = colors)
  }
  x_lims <- c(
    min(preds$low_ci) - 1 * sd(preds$low_ci),
    max(preds$high_ci) + 1 * sd(preds$high_ci)
  )

  p <- ggplot2::ggplot(data = preds, ggplot2::aes(
    x = .data$estimate,
    y = .data$history,
    color = .data$dose
  )) +
    ggplot2::geom_point(size = 5) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = .data$low_ci,
        xmax = .data$high_ci
      ),
      height = 0.6
    ) +
    color_scale +
    ggplot2::scale_y_discrete(
      limits = preds$history,
      expand = c(0, 0.2)
    ) +
    ggplot2::scale_x_continuous(limits = x_lims) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Dosage")) +
    ggplot2::labs(
      x = paste0("Average Predicted ", out_lab, " Value"),
      y = paste0(exp_lab, " Exposure History")
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 18),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )

  if (isTRUE(save.out) || is.character(save.out)) {
    home_dir <- obj[["home_dir"]]
    out_dir <- fs::path_join(c(home_dir, "histories", "plots"))
    .create_dir_if_needed(out_dir)

    if (is.character(save.out)) {
      file_name = save.out
    } else {
      file_name <- sprintf(
        "comparisons_plot-outcome_%s-exposure_%s.jpeg", 
        gsub("\\.", "\\_", outcome), 
        exposure_root
      )
    }
    
    out <- fs::path_join(c(out_dir, file_name))
    cat(sprintf(
      '\nSaving comparisons plot to `jpeg` file:\n%s\n',
      out
    ))
    ggplot2::ggsave(out, plot = p, height = 8, width = 14)
  }

  return(p)
}

#' @rdname compareHistories
#'
#' @param object devMSM_histories object from [compareHistories()]
#' @param type Either "preds" or "comps" corresponding to the
#'  results of [marginaleffects::avg_predictions()] at low and high dosages or
#'  [marginaleffects::avg_comparisons()] respectively
#' @param ... ignored
#'
#' @export
summary.devMSM_comparisons <- function(object, type = "comps", ...) {
  dreamerr::check_arg(type, "scalar character")
  type <- match.arg(type, c("preds", "comps"))
  return(object[[type]])
}

#' Visualize distribution of sample across exposure histories
#'
#' Create customized, user-specified exposure histories and tables displaying
#' sample distribution across them for user inspection.
#'
#' @param epoch_history character vector of length `nrow(data)` containing epochi histories for 
#'  each unit 
#' @param epoch character vector of epoch values
#' 
#' @noRd
.print_eval_hist <- function(epoch_history, epoch, hi_lo_cut, reference = NULL, comparison = NULL) {
  
  if (length(hi_lo_cut) == 1L) hi_lo_cut <- c(hi_lo_cut, hi_lo_cut)
  else hi_lo_cut <- sort(hi_lo_cut)
  
  epoch_history_tab <- data.frame(table(epoch_history, useNA = "ifany"))
  colnames(epoch_history_tab) <- c("epoch_history", "n")
  n_total <- sum(epoch_history_tab$n)
  n_total_hist <- nrow(epoch_history_tab)
  
  if (!is.null(reference) && !is.null(comparison)) {
    keep_idx <- epoch_history_tab$epoch_history %in% c(reference, comparison)
    if (!any(keep_idx)) {
      warning("There are no participants in your sample with the reference/comparison histories you specified, using high/low cutoff (if applicable).", call. = FALSE)
    }
    epoch_history_tab <- epoch_history_tab[keep_idx, , drop = FALSE]
    if (any(!c(reference, comparison) %in% epoch_history_tab$epoch_history)){
      #NG: maybe should be a warning?
      warning(sprintf("There are no participants in your sample in the following histories: %s. 
                   Please revise your reference/comparison histories and/or the high/low cutoffs, if applicable.",
                   paste(c(reference, comparison)[!c(reference,comparison) %in% epoch_history_tab$epoch_history], collapse = ", ")),
           call. = FALSE)
    }
  }
  n_included <- sum(epoch_history_tab$n)
  n_included_hist <- nrow(epoch_history_tab)
  
  message(sprintf(
    "USER ALERT: Out of the total of %s individuals in the sample, below is the distribution of the %s (%.0f%%) individuals that fall into %s user-selected exposure histories (out of the %s total) created from %sth and %sth percentile values for low and high levels of exposure-epoch %s. \n",
    n_total,
    n_included,
    n_included / n_total * 100,
    n_included_hist,
    n_total_hist,
    hi_lo_cut[1] * 100,
    hi_lo_cut[2] * 100,
    paste(unique(epoch), collapse = ", ")
  ))
  
  message("USER ALERT: Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision:\n")
  
  tab_caption <- sprintf(
    "Summary of user-selected exposure histories based on exposure main effects %s:",
    paste(unique(epoch), collapse = ", ")
  )
  tab <- tinytable::tt(epoch_history_tab, caption = tab_caption)  
  print(tab, "markdown")
  cat("\n")
  
  return(invisible(NULL))
}
