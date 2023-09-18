
#' Trim IPTW balancing weights
#'
#' Trims IPTW balancing weights with heavy right tails by populating all weight
#' values above a given quantile with the weight value of that quantile.
#'
#' @seealso {[WeightIt::trim()],
#'   <https://search.r-project.org/CRAN/refmans/WeightIt/html/trim.html>}
#' @param home_dir path to home directory
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
#'                    tv_confounders = c("B.1", "B.2", "B.3"),
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



trimWeights <- function(home_dir, exposure, outcome, weights, quantile = 0.95, verbose = TRUE, save.out = TRUE){

  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if(!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.", call. = FALSE)
    }
  }

  if (missing(weights)){
    stop("Please supply a list of IPTW weights to trim.", call. = FALSE)
  }
  if(!is.numeric(quantile)){
    stop('Please sprovide a numeric quantile value between 0 and 1.', call. = FALSE)
  }
  else if (quantile > 1 || quantile < 0) {
    stop('Please provide a quantile value between 0 and 1.', call. = FALSE)
  }

  if (!inherits(weights, "list")){
    stop("Please supply a list of weights output from the createWeights function.", call. = FALSE)
  }

  if(save.out){
    # creating directories
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

    trim_weights <- lapply(seq_len(length(weights)), function(x){
      w <- weights[[x]]
      t <- WeightIt::trim(w$weights, at = quantile)

      if (verbose){
      cat('\n')
      cat(paste0("USER ALERT: For imputation ", x, " and the ", exposure, "-", outcome, " relation, following trimming at the ",
                     quantile, " quantile, the median weight value is ", round(median(t), 2) ,
                     " (SD= ", round(sd(t), 2), "; range= ", round(min(t), 2), "-", round(max(t), 2), ")."), "\n")
      cat('\n')
      }

      # Save histogram of new weights
      p <- ggplot2::ggplot(as.data.frame(t), aes(x = t)) +
        ggplot2::geom_histogram(color = 'black', bins = 15)

      if(verbose){
        print(p)
      }

      if(save.out){
        ggplot2::ggsave(paste("Hist_", exposure, "-", outcome, "_", weights[[x]]$method, "_weights_trim_",
                              quantile, "_imp_", x, ".png", sep = ""), path = paste0(home_dir, "/weights/histograms/"), plot = p,
                        height = 8, width = 14)
      }

      # w$weights <- NA
      is.na(w$weights) <- TRUE
      w$weights <- t
      w

    })
  }

  # df
  else if ( !is.null(names(weights)) ){
    trim_weights <- lapply(1, function(x){
      w <- weights[[1]]
      t <- WeightIt::trim(w$weights, at = quantile)

      if (verbose){
      cat('\n')
      cat(paste0("USER ALERT: For the ", exposure, "-", outcome, " relation, following trimming at the ",
                     quantile, " quantile, the median weight value is ", round(median(t), 2) ,
                     " (SD= ", round(sd(t), 2), "; range= ", round(min(t), 2), "-", round(max(t), 2), ")."), "\n")
      cat('\n')
      }

      # Save histogram of new weights
      p <- ggplot2::ggplot(as.data.frame(t), aes(x = t)) +
        ggplot2::geom_histogram(color = 'black', bins = 15) +
        ggtitle(paste0("Weights trimmed at the ", quantile, "th value"))

      if(verbose){
        print(p)
      }

      if(save.out){
        ggplot2::ggsave(paste("Hist_", exposure, "-", outcome,"_",weights[[x]]$method, "_weights_trim_", quantile, ".png", sep = ""), plot = p,
                        path = paste0(home_dir, "/weights/histograms/"), height = 8, width = 14)
      }

        w$weights <- NA
        w$weights <- t
        w
    })
    names(trim_weights) <- "0"
  }

  if(save.out){
    # Save truncated weight data
    saveRDS(trim_weights, paste0(home_dir, "/weights/values/", exposure, "-", outcome, "_", weights[[1]]$method, "_weights_trim.rds"))
  }

  trim_weights
}



