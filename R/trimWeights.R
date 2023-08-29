
#' Trim IPTW balancing weights
#'
#' @seealso {[WeightIt::trim()], <url1>}
#' @param home_dir path to home directory
#' @param weights list of IPTW weights output from createWeights()
#' @param quantile (optional) numeric value between 0 and 1 of quantile value at which to trim weights (default is 0.95)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is TRUE)
#'
#' @return list of model output with trimmed weights
#' @export
#'
#' @examples
trimWeights <- function(home_dir, weights, quantile = 0.95, verbose = TRUE){

  if (missing(home_dir)){
    stop("Please supply a home directory.", call. = FALSE)
  }
  if (missing(weights)){
    stop("Please supply a list of IPTW weights.", call. = FALSE)
  }

  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.", call. = FALSE)
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
      ggplot2::ggplot(as.data.frame(t), aes(x = t)) +
        ggplot2::geom_histogram(color = 'black', bins = 15)
        ggplot2::ggsave(paste("Hist_", exposure, "-", outcome, "_", weights[[x]]$method, "_weights_trim_",
                              quantile, "_imp_", x, ".png", sep = ""), path = paste0(home_dir, "/weights/histograms/"),
                        height = 8, width = 14)

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
      ggplot2::ggplot(as.data.frame(t), aes(x = t)) +
        ggplot2::geom_histogram(color = 'black', bins = 15)
        ggplot2::ggsave(paste("Hist_", exposure, "-", outcome,"_",weights[[x]]$method, "_weights_trim_", quantile, ".png", sep = ""),
                        path = paste0(home_dir, "/weights/histograms/"), height = 8, width = 14)

        w$weights <- NA
        w$weights <- t
        w
    })
    names(trim_weights) <- "0"
  }

  # Save truncated weight data
  saveRDS(trim_weights, paste0(home_dir, "/weights/values/", exposure, "-", outcome, "_", weights[[1]]$method, "_weights_trim.rds"))

  trim_weights
}



