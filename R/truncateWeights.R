#' Truncates weights
#'
#'Code to 'cut off' weights at 90th percentile and populate all of those above at 90th percentile value to avoid heavy tails that can bias results
#'
#' @return dat_w_t
#' @param data_for_model_with_weights output from condenseWeights
#' @param object msm object that contains all relevant user inputs
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_histogram
#' @importFrom WeightIt trim
#' @seealso [condenseWeights()]
#' @examples truncateWeights(object, data_for_model_with_weights)
#'
truncateWeights <-function(object, data_for_model_with_weights){

  home_dir=object$home_dir
  exposure=object$exposure
  exposure_epochs=object$exposure_epochs
  outcome=object$outcome
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))
  m=object$m
  ID=object$ID
  weights_method=object$weights_method

  if (weights_percentile_cutoff>1 | weights_percentile_cutoff<0){
    stop('Please select a percentile cutoff between 0 and 1 in the msmObject')}

  dat_w=data_for_model_with_weights

  #creating list of all cutoff values to cycle through (user-specified plus 2 others for sens checks)
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  data_for_model_with_weights_cutoff=list()
  ids=as.data.frame(dat_w[1])
  ids=as.data.frame(ids[,ID])
  ids=ids[!duplicated(ids),1]
  truncated_weights=data.frame(expand.grid(ids, 1:m))
  colnames(truncated_weights)<-c(ID,".imp")


  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  all_trunc_weights=data.frame()
  data_all=data_for_model_with_weights

  dat_w_t<- lapply(all_cutoffs, function(c){ #cycling through cutoff value
    cat("\n")
    cat(paste0("** Truncation Value ", c, " **"), "\n")

    lapply(1:length(dat_w), function(t){ #cycling through imputations
      cutoff=c
      k=t
      data=dat_w[[t]]
      name=paste0(cutoff, "_weight_cutoff") #labeling new weight column
      cutoff_weights=WeightIt::trim(data[,paste0("weights")], cutoff, lower=F)

      cat(paste0("For ", exposure, "-", outcome, ", using the ", weights_method, " weights method and user-specified cutoff percentile of ", cutoff, ", the median weight is ", round(median(as.numeric(unlist(cutoff_weights))),2) ,
                 " (SD= ",round(sd(as.numeric(unlist(cutoff_weights))),2), "; range= ", round(min(as.numeric(unlist(cutoff_weights))),2), " - ", round(max(as.numeric(unlist(cutoff_weights))),2), ")"), "\n")
      # cat("\n")

      cutoff_weights=data.frame(x=cutoff_weights)
      colnames(cutoff_weights)=name
      new=data.frame(ID=data[,ID])

      #adds exposure epochs
      #calculates the mean value for each exposure for each exposure epoch
      for (e in 1:nrow(exposure_epochs)){
        epoch=exposure_epochs[e,1]
        temp=data.frame(row.names=1:nrow(data))
        new_var=paste0(exposure, "_", epoch)
        #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
        for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
          level=as.numeric(unlist(exposure_epochs[e,2]))[l]
          z=data[,which(grepl(paste0(exposure, ".", level), names(data)))]
          temp=cbind(temp, z)
        }
        #adds a new variable of the exposure averaged within epoch
        new=new%>%dplyr::mutate(.imp=k,
                                !!new_var :=rowMeans(temp, na.rm=T))
        new[,new_var]=as.numeric(new[,new_var])
      }
      new=cbind(new, cutoff_weights)
      data=cbind(data, new)

      #print histogram of new weights
      ggplot2::ggplot(data=as.data.frame(data), aes(x = as.numeric(data[,name]))) +
        ggplot2::geom_histogram(color = 'black', bins = 15)
      ggplot2::ggsave(paste("Hist_", exposure, "-", outcome, "_weights_cutoff_", cutoff, "_", weights_method, "_imp_", k, ".png", sep=""),
                      path=paste0(home_dir, "final weights/histograms/"), height=8, width=14)
      cat(paste0("A histogram with weights truncated at ", cutoff, " with ", weights_method, " for ", "imputation ", k, " have now been saved in the 'final weights/histograms/' folder."), "\n")
      write.csv(data,paste0(home_dir, "final weights/", exposure, "-", outcome, "_weights_cutoff_", cutoff, "_", weights_method, "_imp_", k,".csv") )

      data
    })
  })
  names(dat_w_t)=all_cutoffs

  saveRDS(dat_w_t, paste0(home_dir, "final weights/values/",  exposure, "-", outcome, "_", weights_method, "_dat_w_t.rds"))

  cat("\n")
  cat("USER ALERT: final truncated weights (using the user-specified cutoff values and 2 other values for subsequent sensiivity analyses) have now each been saved as a dataset in 'final weights' folder","\n")
  cat(paste0("A histogram of weights for exposure ", exposure, " on ", outcome, ", all imputations has been saved to the 'final weights/histograms/' folder"),"\n")
  cat("\n")
  cat("\n")

  return(dat_w_t)
}
