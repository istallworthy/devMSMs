#Code to 'cut off' weights at 90th percentile and populate all of those above at 90th percentile value

#' truncateWeights
#'
#' @return data_for_model_with_weights_cutoff
#' @param data_for_model_with_weights output from condenseWeights
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param percentile_cutoff percentile value for cutting off and replacing heavy tails of weights
#' @export
#' @importFrom dplyr mutate
#' @importFrom ggplot2
#' @examples
#'
truncateWeights <-function(data_for_model_with_weights, exposures, percentile_cutoff=0.90){

  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  #loops through all treatments and creates tx_weight_cutoff variable: replaces any weights above the 90th quantile with the 90th quantile value
  for (t in 1:length(exposures)){
    exposure=exposures[t]
    new_var=paste0(exposure, "_weight_cutoff")
    data_for_model_with_weights_cutoff=data_for_model_with_weights%>%
      dplyr::mutate(!!new_var := ifelse(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")]<
                                          quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")], probs=percentile_cutoff),
                                        data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")],
                                        quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")], probs=percentile_cutoff)
      ) )

    #print histogram of new weights
    ggplot2::ggplot(data=as.data.frame(data_for_model_with_weights_cutoff), aes(x = as.numeric(data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==paste0(exposure, "_weight_cutoff")]))) +
      geom_histogram(color = 'black', bins = 15)
    ggplot2::ggsave(paste("Hist_", exposure, "_weights_cutoff_", percentile_cutoff, ".png", sep=""), path=paste0(home_dir, "final weights"), height=8, width=14)
    message("Histograms of cut off weights per exposure have been saved to the 'final weights' folder")

  }

  return(data_for_model_with_weights_cutoff)

}
