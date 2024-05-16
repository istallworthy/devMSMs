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
#' @param hi_lo_cut (optional) list of two numbers indicating quantile values
#'   that reflect high and low values, respectively, for continuous exposure
#'   (default is median split)
#' @param reference (optional) list sof one or more strings of "-"-separated "l"
#'   and "h" values indicative of a reference exposure history to which to
#'   compare comparison, required if comparison is supplied
#' @param comparison (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is supplied
#' @param mc_comp_method (optional) character abbreviation for multiple
#'   comparison correction method for stats::p.adjust, default is
#'   Benjamini-Hochburg ("BH")
#' @param dose_level (optional) "l" or "h" indicating whether low or high doses
#'   should be tallied in tables and plots (default is high "h")
#'
#' @return list containing two dataframes: `preds` with predictions from
#'  [marginaleffects::avg_predictions()] containing average expected outcome
#'  for different exposure histories and `comps` with contrasts from
#'  [marginaleffects::comparisons()] comparing different exposure history
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
#' w <- createWeights(data = data, obj = obj, formulas = f)
#' fit <- fitModel(
#'   data = data, obj = obj, weights = w, 
#'   outcome = "D.3", model = "m0"
#' )
#' 
#' comp = compareHistories(
#'   obj, fit = fit,
#'   hi_lo_cut = c(0.3, 0.6)
#' )
#' print(comp)
#' plot(comp)
#' summary(comp, "preds")
#' summary(comp, "comps")
#' 
#' comp2 = compareHistories(
#'   obj, fit = fit, 
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
    obj, fit,
    hi_lo_cut = c(0.3, 0.6), dose_level = c("h", "l"),
    reference = NULL, comparison = NULL,
    mc_comp_method = stats::p.adjust.methods,
    verbose = FALSE, save.out = FALSE, home_dir) {
  ### Checks ----
  dreamerr::check_arg(verbose, save.out, "scalar logical")
  if (save.out) {
    dreamerr::check_arg_plus(home_dir, "path dir")
    .create_dir_if_needed(file.path(home_dir, "histories"))
    .create_dir_if_needed(file.path(home_dir, "plots"))
  }
  if (verbose) {
    rlang::check_installed("tinytable")
  }

  dreamerr::check_arg(fit, "class(devMSM_models)")

  # Get objects from `obj`
  exposure <- attr(obj, "exposure")
  exposure_time_pts <- attr(obj, "exposure_time_pts")
  exposure_type <- attr(obj, "exposure_type")
  epoch <- attr(obj, "epoch")
  sep <- attr(obj, "sep")

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
  if ((is.null(comparison) && !is.null(reference)) || (!is.null(comparison) && is.null(reference) == 1)) {
    stop("If either comparison or reference is specified, then both must be specified", call. = FALSE)
  }
  is_invalid_comparison <- !(comparison %in% exposure_levels)
  if (any(is_invalid_comparison)) {
    stop(sprintf("The following elements of `comparison` are invalid: %s", paste0(comparison[is_invalid_comparison], collapse = ", ")), call. = FALSE)
  }
  is_invalid_reference <- !(reference %in% exposure_levels)
  if (any(is_invalid_reference)) {
    stop(sprintf("The following elements of `reference` are invalid: %s", paste0(reference[is_invalid_reference], collapse = ", ")), call. = FALSE)
  }

  # Multiple comparison methods
  mc_comp_method <- match.arg(mc_comp_method)
  dreamerr::check_arg(mc_comp_method, "scalar character | NULL")
  if (is.null(mc_comp_method)) {
    mc_comp_method <- "BH"
  }
  if (!(mc_comp_method %in% stats::p.adjust.methods)) {
    stop("Please provide a single valid character string abbreviation for a multiple comparison method compatible with the stats::p.adjust() function.", call. = FALSE)
  }

  dreamerr::check_arg(hi_lo_cut, "NULL | vector numeric len(2) GE{0} LE{1}")
  hi_lo_cut <- c(min(hi_lo_cut), max(hi_lo_cut)) # sort

  # Dose_level
  dose_level <- match.arg(dose_level, c("h", "l"))


  # STEP 1 ----
  # Define variables for average partial effects

  # TODO: This uses a ton of data
  # getting data to use for determining hi/lo values: should have any epochs created / used in the model
  data <- do.call(rbind, lapply(fit, function(z) z[["data"]]))

  epoch_vars <- exposure
  if (any(exposure != epoch)) {
    epoch_vars <- .get_epoch_var_names(exposure, epoch, sep = sep)
  }

  # TODO: Error out if epoch_history has no values in `table(epoch_history)`?
  mat <- fit[[1]][["data"]][, epoch_vars]
  epoch_history <- .characterize_exposure(mat, exposure_type)

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
    w <- m$weights

    p <- marginaleffects::avg_predictions(
      m,
      newdata = d, variables = prediction_vars, wts = w
    )
    p <- .add_histories(p, epoch_vars)
    p
  })


  # STEP 3 ----
  # Conduct history comparisons
  # All pairwise comparisons if no ref/comparisons were specified by user
  if (is.null(reference) && is.null(comparison)) {
    # Pairwise comparision, but we want term to be nicely labeled by histories
    comps <- lapply(preds, function(y) {
      marginaleffects::hypotheses(
        model = as.data.frame(y), vcov = vcov(y), hypothesis = "pairwise"
      )
    })
  } else {
    hypothesis <- .create_hypotheses_mat(preds[[1]]$term, reference, comparison)
    comps <- lapply(preds, function(y) {
      marginaleffects::hypotheses(y, hypothesis = hypothesis)
    })
  }


  # STEP 4 ----
  # pooling predicted values and contrasts for imputed data
  is_pooled <- length(preds) > 1
  if (is_pooled) { # IMPUTED DATA
    rlang::check_installed("mice")
    # TODO: with multiply imputed data, broom::glimpse() is bugging out
    preds <- summary(mice::pool(preds), digits = 3)
    comps <- summary(mice::pool(comps))
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
  class(res) <- c("devMSM_comparisons", "list")
  attr(res, "outcome") <- attr(fit, "outcome")
  attr(res, "hi_lo_cut") <- hi_lo_cut
  attr(res, "is_pooled") <- is_pooled
  attr(res, "epoch_d") <- epoch_d
  attr(res, "epoch_history") <- epoch_history
  attr(res, "exposure") <- exposure
  attr(res, "epoch") <- epoch
  attr(res, "mc_comp_method") <- mc_comp_method
  attr(res, "exp_lab") <- exp_lab
  if (!is.null(reference) && !is.null(comparison)) {
    attr(res, "reference") <- reference
    attr(res, "comparison") <- comparison
  }

  if (verbose) print(res)
  return(res)
}


#' @rdname compareHistories
#'
#' @param x devMSM_histories object from [compareHistories()]
#' @param ... ignored
#'
#' @export
print.devMSM_comparisons <- function(x, ...) {
  preds <- x$preds
  comps <- x$comps
  epoch <- attr(x, "epoch")
  hi_lo_cut <- attr(x, "hi_lo_cut")
  epoch_history <- attr(x, "epoch_history")
  mc_comp_method <- attr(x, "mc_comp_method")

  reference <- NULL
  comparison <- NULL
  if (!is.null(attr(x, "reference"))) {
    reference <- attr(x, "reference")
  }
  if (!is.null(attr(x, "comparison"))) {
    comparison <- attr(x, "comparison")
  }

  columns_to_drop <- c("statistic", "s.value", "conf.low", "conf.high", "dose")
  preds_tab <- tinytable::tt(preds[, setdiff(colnames(preds), columns_to_drop)])
  preds_tab <- tinytable::format_tt(preds_tab, digits = 2)

  columns_to_drop <- c("statistic", "s.value", "conf.low", "conf.high", "dose")
  comps_tab <- tinytable::tt(comps[, setdiff(colnames(comps), columns_to_drop)])

  cat("Summary of Exposure Main Effects:\n")
  print_eval_hist(epoch_history, epoch, hi_lo_cut, reference, comparison)
  cat("Below are the pooled average predictions by user-specified history:")
  print(preds_tab, "markdown")
  cat("\n")
  cat(sprintf("Conducting multiple comparison correction for all pairings between comparison histories and each refernece history using the %s method. \n", mc_comp_method))
  cat("\n")
  print(comps_tab, "markdown")
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
plot.devMSM_comparisons <- function(x, colors = "Dark2", exp_lab = NULL, out_lab = NULL, save.out = FALSE, home_dir = NULL, ...) {
  dreamerr::check_arg(exp_lab, out_lab, "NULL | character scalar")
  if (is.null(out_lab)) out_lab <- attr(x, "outcome")
  if (is.null(exp_lab)) exp_lab <- attr(x, "exp_lab")

  exposure <- attr(x, "exposure")
  epoch <- attr(x, "epoch")
  outcome <- attr(x, "outcome")

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
  preds$low_ci <- preds$estimate - (1.96 * preds$std.error)
  preds$high_ci <- preds$estimate + (1.96 * preds$std.error)
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
      text = ggplot2::element_text(size = 14),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )

  if (save.out) {
    ggplot2::ggsave(
      file.path(
        home_dir, "plots",
        sprintf(
          "%s-%s.jpeg",
          exposure, outcome
        )
      ),
      plot = p
    )
    cat("\n")
    cat("See the '/plots/' folder for graphical representations of results.")
  }

  return(p)
}

#' @rdname compareHistories
#'
#' @param object devMSM_histories object from [compareHistories()]
#' @param type Either "preds" or "comps" corresponding to the
#'  results of [marginaleffects::avg_predictions()] at low and high dosages or
#'  [marginaleffects::comparisons()] respectively
#' @param ... ignored
#'
#' @export
summary.devMSM_comparisons <- function(object, type = c("preds", "comps"), ...) {
  type <- match.arg(type)
  return(object[[type]])
}
