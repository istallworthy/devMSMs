#Code to 'cut off' weights at 90th percentile and populate all of those above at 90th percentile value

#' Title
#'
#' @return
#' @export
#' @importFrom dplyr mutate
#' @examples
#'
truncateWeights <-function(exposures, percentile_cutoff=0.90){

  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  #loops through all treatments and creates tx_weight_cutoff variable: replaces any weights above the 90th quantile with the 90th quantile value
  for (t in 1:length(exposures)){
    exposure=exposures[t]
    new_var=paste0(exposure, "_weight_cutoff")
    data_for_model=data_for_model%>%dplyr::mutate(!!new_var := ifelse(
      data_for_model[,which(grepl(paste0(exposure, "_weight"), names(data_for_model)))]<
        quantile(data_for_model[,which(grepl(paste0(exposure, "_weight"), names(data_for_model)))], probs=percentile_cutoff),
      data_for_model[,which(grepl(paste0(exposure, "_weight"), names(data_for_model)))],
      quantile(data_for_model[,which(grepl(paste0(exposure, "_weight"), names(data_for_model)))], probs=percentile_cutoff)
    ))

    return(data_for_model)

  }

#plot

}
