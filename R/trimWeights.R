
#' Trim IPTW balancing weights
#'
#' Trims IPTW balancing weights with heavy right tails by populating all weight
#' values above a given quantile with the weight value of that quantile.
#'
#' @seealso [WeightIt::trim()],
#'   <https://search.r-project.org/CRAN/refmans/WeightIt/html/trim.html> which
#'   this function wraps
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param weights list of IPTW weights output from createWeights()
#' @param quantile (optional) numeric value between 0 and 1 of quantile value at
#'   which to trim weights (default is 0.95)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return list of model output with trimmed weights
#' @export
#' @examples
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
#'
#' test <- data.frame(ID = 1:50,
#'                    A.1 = rnorm(n = 50),
#'                    A.2 = rnorm(n = 50),
#'                    A.3 = rnorm(n = 50),
#'                    B.1 = rnorm(n = 50),
#'                    B.2 = rnorm(n = 50),
#'                    B.3 = rnorm(n = 50),
#'                    C = rnorm(n = 50),
#'                    D.3 = rnorm(n = 50))
#' test[, c("A.1", "A.2", "A.3")] <- lapply(test[, c("A.1", "A.2", "A.3")], as.numeric)
#'
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' t <- trimWeights(exposure = "A",
#'                  outcome = "D.3",
#'                  weights = w,
#'                  save.out = FALSE)
#' t <- trimWeights(exposure = "A",
#'                  outcome = "D.3",
#'                  weights = w,
#'                  quantile = 0.75,
#'                  save.out = FALSE)



trimWeights <- function(home_dir, exposure, outcome, weights, quantile = 0.95, verbose = TRUE, save.out = TRUE) {

  if (save.out) {
    if (missing(home_dir)) {
      stop ("Please supply a home directory.",
           call. = FALSE)
    }
    else if (!is.character(home_dir)) {
      stop ("Please provide a valid home directory path as a string if you wish to save output locally.",
            call. = FALSE)
    }
    else if (!dir.exists(home_dir)) {
      stop ("Please provide a valid home directory path if you wish to save output locally.",
            call. = FALSE)
    }
  }

  if (missing(exposure)) {
    stop ("Please supply a single exposure.",
          call. = FALSE)
  }
  else if (!is.character(exposure) || length(exposure) != 1) {
    stop ("Please supply a single exposure as a character.",
          call. = FALSE)
  }

  if (missing(weights)) {
    stop ("Please supply a list of IPTW weights to trim.",
          call. = FALSE)
  }
  else if (!is.list(weights) || is.data.frame(weights)) {
    stop ("Please supply a list of weights output from the createWeights function.",
         call. = FALSE)
  }
  else if (is.list(weights) && !is.data.frame(weights)) {
    if (sum(sapply(weights, function(x) {
      inherits(x, "weightitMSM")})) != length(weights)) {
      stop ("Please supply a list of weights output from the createWeights function.",
            call. = FALSE)
    }
  }

  if (missing(outcome)) {
    stop ("Please supply a single outcome.",
          call. = FALSE)
  }
  else if (!is.character(outcome) || length(outcome) != 1) {
    stop ("Please supply a single outcome as a character.",
          call. = FALSE)
  }

  if (!is.numeric(quantile) || length(quantile) != 1) {
    stop ('Please provide a single numeric quantile value between 0 and 1.',
          call. = FALSE)
  }
  else if (quantile > 1 || quantile < 0) {
    stop ('Please provide a quantile value between 0 and 1.',
          call. = FALSE)
  }

  if (!is.logical(verbose)) {
    stop ("Please set verbose to either TRUE or FALSE.",
          call. = FALSE)
  }
  else if (length(verbose) != 1) {
    stop ("Please provide a single TRUE or FALSE value to verbose.",
          call. = FALSE)
  }

  if (!is.logical(save.out)) {
    stop ("Please set save.out to either TRUE or FALSE.",
          call. = FALSE)
  }
  else if (length(save.out) != 1) {
    stop ("Please provide a single TRUE or FALSE value to save.out.",
          call. = FALSE)
  }


  if (save.out) {
    weights_dir <- file.path(home_dir, "weights")
    if (!dir.exists(weights_dir)) {
      dir.create(weights_dir)
    }
    values_dir <- file.path(home_dir, "weights", "values")
    if (!dir.exists(values_dir)) {
      dir.create(values_dir)
    }
    hist_dir <- file.path(home_dir, "weights", "histograms")
    if (!dir.exists(hist_dir)) {
      dir.create(hist_dir)
    }
  }

  #imputed data

  if (is.null(names(weights))) {

    trim_weights <- lapply(seq_len(length(weights)), function(x) {
      w <- weights[[x]]
      t <- WeightIt::trim(w$weights, at = quantile)

      if (verbose) {
        cat('\n')
        cat(sprintf("For imputation %s and the %s-%s relation, following trimming at the %s quantile, the median weight value is
                  %s (SD= %s; range= %s-%s). \n",
                    x,
                    exposure,
                    outcome,
                    quantile,
                    round(median(t), 2),
                    round(sd(t), 2),
                    round(min(t), 2),
                    round(max(t))))
        cat('\n')
      }

      # Save histogram of new weights

      p <- ggplot2::ggplot(as.data.frame(t), ggplot2::aes(x = t)) +
        ggplot2::geom_histogram(color = 'black',
                                bins = 15) +
        ggplot2::ggtitle(sprintf("Weights trimmed at the %s th value",
                                 quantile))

      if (verbose) {
        print(p)
      }

      if (save.out) {
        ggplot2::ggsave(sprintf("Hist_%s-%s_%s_weights_trim_%s_imp_%s.png",
                                exposure, outcome, weights[[x]]$method, quantile, x),
                        path = sprintf("%s/weights/histograms/",
                                       home_dir),
                        plot = p,
                        height = 8, width = 14)
      }

      is.na(w$weights) <- TRUE
      w$weights <- t
      w

    } )
  }

  # df

  else if ( !is.null(names(weights)) ) {
    trim_weights <- lapply(1, function(x) {
      w <- weights[[1]]
      t <- WeightIt::trim(w$weights, at = quantile)

      if (verbose) {
        cat('\n')
        cat(sprintf("For the %s-%s relation, following trimming at the %s quantile, the median weight value is
                  %s (SD= %s; range= %s-%s). \n",
                    exposure,
                    outcome,
                    quantile,
                    round(median(t), 2),
                    round(sd(t), 2),
                    round(min(t), 2),
                    round(max(t))))
        cat('\n')
      }

      # Save histogram of new weights
      p <- ggplot2::ggplot(as.data.frame(t), ggplot2::aes(x = t)) +
        ggplot2::geom_histogram(color = 'black',
                                bins = 15) +
        ggplot2::ggtitle(sprintf("Weights trimmed at the %sth value",
                                 quantile))

      if (verbose) {
        print(p)
      }

      if (save.out) {

        ggplot2::ggsave(sprintf("Hist_%s-%s_%s_weights_trim_%s.png",
                                exposure, outcome, weights[[x]]$method, quantile),
                        path = sprintf("%s/weights/histograms/",
                                       home_dir),
                        plot = p,
                        height = 8,
                        width = 14)
      }

      w$weights <- NA
      w$weights <- t
      w
    } )
    names(trim_weights) <- "0"
  }

  if (save.out) {
    saveRDS(trim_weights,
            sprintf("%s/weights/values/%s-%s_%s_weights_trim.rds",
                    home_dir, exposure, outcome,weights[[1]]$method ))
  }

  trim_weights
}



