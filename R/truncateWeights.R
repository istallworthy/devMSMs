#' Truncates weights
#'
#'Code to 'cut off' weights at 90th percentile and populate all of those above at 90th percentile value to avoid heavy tails that can bias results
#'
#' @return data_for_model_with_weights_cutoff
#' @param data_for_model_with_weights output from condenseWeights
#' @param object msm object that contains all relevant user inputs
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @seealso [condenseWeights()]
#' @examples truncateWeights(data_for_model_with_weights, home_dir, exposures)
#'
truncateWeights <-function(object, data_for_model_with_weights){

  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))

  if (weights_percentile_cutoff>1 | weights_percentile_cutoff<0){
    stop('Please select a percentile cutoff between 0 and 1 in the msmObject')}


  #creating list of all cutoff values to cycle through
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  data_for_model_with_weights_cutoff=data_for_model_with_weights
  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  #loops through all treatments and creates tx_weight_cutoff variable: replaces any weights above the 90th quantile with the 90th quantile value

  for (c in 1:length(all_cutoffs)){ #cycling through all cutoff values
    cutoff=all_cutoffs[c]

    for (z in seq(length(outcomes))){ #cycling through all outcomes
      outcome=outcomes[z]

      for (t in 1:length(exposures)){ #cycling through all exposures
        exposure=exposures[t]

        new_var=paste0(exposure, "-", outcome,"_", cutoff, "_weight_cutoff") #labeling new weight column
        cutoff_weights=data_for_model_with_weights%>% #conducting perecentile cutoff and top coding
          dplyr::mutate(!!new_var := ifelse(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "-", outcome, "_weight")]<
                                              quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "-", outcome, "_weight")], probs=cutoff),
                                            data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure,"-", outcome,"_weight")],
                                            quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "-", outcome, "_weight")], probs=cutoff)
          ) )%>%dplyr::select(all_of(new_var))

        if(cutoff==weights_percentile_cutoff){ #prints only values for user-specified cutoff
          cat(paste0("For ", exposure, "-", outcome, ", using the user-specified cutoff percentile of ", cutoff, ", the median weight is ", median(as.numeric(unlist(cutoff_weights))) ,
                     " (SD= ",sd(as.numeric(unlist(cutoff_weights))), "; range= ", min(as.numeric(unlist(cutoff_weights))), " - ", max(as.numeric(unlist(cutoff_weights))), ")"), "\n")
        }

        #binding on new cutoff weight
        data_for_model_with_weights_cutoff=cbind(data_for_model_with_weights_cutoff, cutoff_weights)


        #print histogram of new weights by cutoff value
        ggplot2::ggplot(data=as.data.frame(data_for_model_with_weights_cutoff), aes(x = as.numeric(data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff")]))) +
          geom_histogram(color = 'black', bins = 15)
        ggplot2::ggsave(paste("Hist_", exposure, "-", outcome, "_weights_cutoff_", cutoff, ".png", sep=""), path=paste0(home_dir, "final weights/histograms/"), height=8, width=14)

        cat(paste0("A histogram of weights cut off at ", cutoff, " for exposure ", exposure, " on ", outcome, " has been saved to the 'final weights/histograms/' folder"),"\n")
        cat("\n")

      }
    }
  }

  # browser()

  write.csv(data_for_model_with_weights_cutoff,paste0(home_dir, "final weights/data_for_model_with_weights_cutoff_", paste(all_cutoffs, sep=",", collapse=" "), ".csv") )
  cat("\n")
  cat("USER ALERT: final cutoff weights using the user-specified cutoff values and 2 other values for subsequent sensiivity analyses have now each been saved as a dataset in 'final weights' folder","\n")
  return(data_for_model_with_weights_cutoff)


}
