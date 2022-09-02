#' Truncates weights
#'
#'Code to 'cut off' weights at 90th percentile and populate all of those above at 90th percentile value to avoid heavy tails that can bias results
#'
#' @return data_for_model_with_weights_cutoff
#' @param data_for_model_with_weights output from condenseWeights
#' @param home_dir path to home directory for the project
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param percentile_cutoff percentile value for cutting off and replacing heavy tails of weights
#' @export
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @seealso [condenseWeights()]
#' @examples truncateWeights(data_for_model_with_weights, home_dir, exposures, percentile_cutoff=0.90)
#'
truncateWeights <-function(data_for_model_with_weights, home_dir, exposures, percentile_cutoff=0.90){

  data_for_model_with_weights_cutoff=data_for_model_with_weights
  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  #loops through all treatments and creates tx_weight_cutoff variable: replaces any weights above the 90th quantile with the 90th quantile value
  for (t in 1:length(exposures)){
    exposure=exposures[t]
    new_var=paste0(exposure, "_weight_cutoff")
    cutoff_weights=data_for_model_with_weights%>%
      dplyr::mutate(!!new_var := ifelse(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")]<
                                          quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")], probs=percentile_cutoff),
                                        data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")],
                                        quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "_weight")], probs=percentile_cutoff)
      ) )%>%dplyr::select(all_of(new_var))

    data_for_model_with_weights_cutoff=cbind(data_for_model_with_weights_cutoff, cutoff_weights)


    #print histogram of new weights
    ggplot2::ggplot(data=as.data.frame(data_for_model_with_weights_cutoff), aes(x = as.numeric(data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==paste0(exposure, "_weight_cutoff")]))) +
      geom_histogram(color = 'black', bins = 15)
    ggplot2::ggsave(paste("Hist_", exposure, "_weights_cutoff_", percentile_cutoff, ".png", sep=""), path=paste0(home_dir, "final weights"), height=8, width=14)

    print(paste0("Histogram of weights cut off at ", percentile_cutoff, " for exposure ", exposure, " has been saved to the 'final weights' folder"))

  }

  write.csv(data_for_model_with_weights_cutoff,paste0(home_dir, "final weights/ data_for_model_with_weights_cutoff_", percentile_cutoff, ".csv") )
  print("USER ALERT: final cutoff weights have now been saved merged into datast in 'final weights' folder")
  return(data_for_model_with_weights_cutoff)

}
