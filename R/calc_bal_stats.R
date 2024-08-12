#' Calculate balance stats based on Jackson paper
#'
#' Calculate weighted or unweighted standardized balance statistics for a given
#' exposure time point, using all relevant confounders. Draws on Jackson, 2016
#' approaches to assessing balance for time-varying exposures by weighting
#' statistics based on sample distribution in exposure histories.
#'
#' @param data data in wide format as: a data frame, path to folder of imputed
#'   .csv files, or mids object
#' @param obj initialized MSM object from `initMSM`
#' @param weights (optional) list of IPTW weights output from createWeights
#' @param balance_thresh (optional) one or two numbers between 0 and 1
#'   indicating a single balancingn threshold or thresholds for more and less
#'   important confounders, respectively
#' @param imp_conf (optional) list of variable names reflecting important
#'   confounders (required if two balance thresholds are provided)
#' @return data frame of balance statistics
#'
#' @noRd
calc_bal_stats <- function(data, obj, weights = NULL, balance_thresh = NULL, imp_conf = NULL) {
  # NOTE:
  # 1. Checks are already done in assessBalance
  # 2. This assumes data is a `data.frame` (and not list/mids) and weights is a `weightitMSM` object (and not devMSM_weights). This is done in `assessBalance`

  # TODO: Calculate effective sample size with `cobalt:ESS`
  ### Calculate balance stats ----
  exposure <- obj[["exposure"]]
  exposure_time_pts <- obj[["exposure_time_pts"]]
  exposure_type <- obj[["exposure_type"]]
  cobalt_cov_fun <- switch(exposure_type,
    "continuous" = cobalt::col_w_cov,
    "binary" = cobalt::col_w_smd
  )
  is_weighted <- !is.null(weights)

  var_tab <- obj[["var_tab"]]
  ti_conf <- var_tab$var[var_tab$type == "ti_conf"]
  ti_conf_time <- var_tab$time[var_tab$type == "ti_conf"]
  tv_conf <- var_tab$var[var_tab$type == "tv_conf"]
  tv_conf_time <- var_tab$time[var_tab$type == "tv_conf"]
  exposure <- var_tab$var[var_tab$type == "exposure"]
  exposure_time <- var_tab$time[var_tab$type == "exposure"]
  sep <- obj[["sep"]]

  # creating initial data frames
  # data frame with all sampling weights for all exposures at all exposure time points for all histories
  all_prop_weights <- all_bal_stats <- vector("list", length(exposure_time_pts))

  for (z in seq_along(exposure_time_pts)) {
    exposure_name <- exposure[z]
    exposure_time_pt <- exposure_time_pts[z]
    lagged_time_pts <- exposure_time_pts[exposure_time_pts < exposure_time_pt]
    lagged_exposure_names <- exposure[exposure_time_pts < exposure_time_pt]

    # All lagged variable (var_time_pt < exposure_time_pt)
    vars <- c(
      ti_conf,
      tv_conf[tv_conf_time < exposure_time_pt],
      exposure[exposure_time_pts < exposure_time_pt]
    )

    var_is_factorable <- vapply(vars, function(v) {
      is.character(data[[v]]) || is.factor(data[[v]])
    }, logical(1L))
    
    if (any(var_is_factorable)) {
      bal_vars <- cobalt::splitfactor(data[, vars, drop = FALSE])

      orig_vars <- colnames(bal_vars)
      for (j in seq_along(orig_vars)) {
        orig_var_name = attr(bal_vars[, j], "split.var")
        if (!is.null(orig_var_name)) {
          orig_vars[j] = orig_var_name
        }
      } 

    } else {
      bal_vars <- data[, vars, drop = FALSE]
      orig_vars <- colnames(data[, vars, drop = FALSE]) # added by IS for case of no factors 
    }

    # Cacluate `bal_stats`

    # GET BALANCE STATISTICS FOR T=1 (when there is no history to weight on)
    if (length(lagged_time_pts) == 0) {
      # GETTING BALANCE STATS FOR T=1 W/ NO EXPOSURE HISTORY
      # (ok to standardize immediately)
      std_bal_stats <- cobalt_cov_fun(
        bal_vars, data[[exposure_name]],
        weights = weights, std = TRUE
      )
      bal_stats <- data.frame(
        covariate = names(std_bal_stats),
        orig_covariate = orig_vars,
        std_bal_stats = as.numeric(std_bal_stats)
      )
      prop_weights <- prop_sum <- NULL
      omitted_histories <- c()
    }

    # GET BALANCE STATISTICS FOR T>1 (when there is a history to weight on)
    if (length(lagged_time_pts) > 0) {
      # ASSIGNING HISTORIES FOR EXP TIME POINTS T>1
      # creating proportion weights based on proportion of individuals in a given exposure history
      history <- .characterize_exposure(
        data[, lagged_exposure_names, drop = FALSE], exposure_type
      )
      prop_sum <- as.data.frame(table(history))
      prop_sum$prop <- prop_sum$Freq / length(history)

      # Removing any histories that only have 1 or 0 person contributing (cannot calc bal stats)
      omitted_histories <- c()
      if (any(prop_sum$Freq <= 1)) {
        omitted_histories <- as.character(
          prop_sum$history[prop_sum$Freq == 1 | prop_sum$Freq == 0]
        )
      } # ends hist exc

      # finding balance by history
      # if weighted, use IPTW weights from weightitmsm and weight by history
      history_bal_stats <- sapply(
        setdiff(prop_sum$history, omitted_histories),
        function(h) {
          which_idx <- (history == h)
          cobalt_cov_fun(
            bal_vars[which_idx, , drop = FALSE],
            data[[exposure_name]][which_idx],
            weights = weights[which_idx],
            std = FALSE # raw difference
          )
        }
      )

      # getting weighted mean across histories (cols), weighting by proportion of those w/ that same history
      weighted_bal_stats <- history_bal_stats %*% prop_sum$prop

      if (exposure_type == "continuous") {
        # unweighted covar sd *
        s1 <- cobalt::col_w_sd(bal_vars)

        # exposure SD at that time pt
        s2 <- cobalt::col_w_sd(data[[exposure_name]])
        s <- s1 * s2
      } else if (exposure_type == "binary") {
        # dividing by pool SD estimate (unadjusted)
        s_0 <- cobalt::col_w_sd(bal_vars, subset = data[[exposure_name]] == 0)
        s_1 <- cobalt::col_w_sd(bal_vars, subset = data[[exposure_name]] == 1)
        s <- sqrt((s_0^2 + s_1^2) / 2)
      }

      # For a weighted_bal_stat of 0, make std stat also 0 so as not to throw an error
      std_bal_stats <- weighted_bal_stats / s
      std_bal_stats[is.nan(std_bal_stats)] <- 0

      bal_stats <- data.frame(
        covariate = rownames(std_bal_stats),
        orig_covariate = orig_vars,
        std_bal_stats = as.numeric(std_bal_stats)
      )
    }

    # adds custom bal thresh info
    if (!is.null(imp_conf)) {
      bal_stats$bal_thresh <- ifelse(
        bal_stats$covariate %in% imp_conf,
        balance_thresh[1], balance_thresh[2]
      )
    } else {
      bal_stats$bal_thresh <- balance_thresh
    }
    bal_stats$balanced <- ifelse(
      abs(bal_stats$std_bal_stats) < bal_stats$bal_thresh,
      1, 0
    )

    bal_stats$exposure <- exposure_name
    bal_stats$exposure_time <- exposure_time_pt
    
    bal_stats$covar_time <- as.numeric(.extract_time_pts_from_vars(bal_stats$orig_covariate, sep = sep))

    # collect proportions for histories at this time point
    all_prop_weights[[z]] <- prop_sum
    all_bal_stats[[z]] <- bal_stats

  } # Ends exposure_time_pts loop
  
  if (length(omitted_histories) == 0) {
    res <- list(all_bal_stats = all_bal_stats, omitted_histories = NULL)
  } else {
    res <- list(all_bal_stats = all_bal_stats, omitted_histories = omitted_histories)
  }
  
  return(res)
}



