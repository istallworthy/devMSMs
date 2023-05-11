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
#' @examples truncateWeights(data_for_model_with_weights, home_dir, exposure)
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

  if (weights_percentile_cutoff>1 | weights_percentile_cutoff<0){
    stop('Please select a percentile cutoff between 0 and 1 in the msmObject')}

  # dat_w=readRDS(paste0(home_dir, "original weights/dat_w.rds", sep=""))
  dat_w=data_for_model_with_weights

  #creating list of all cutoff values to cycle through
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  data_for_model_with_weights_cutoff=list()
  # ids=as.data.frame(mice::complete(imp_data_w,1)[,ID])
  ids=as.data.frame(dat_w[1])
  ids=as.data.frame(ids[,ID])
  ids=ids[!duplicated(ids),1]
  truncated_weights=data.frame(expand.grid(ids, 1:m))
  colnames(truncated_weights)<-c(ID,".imp")


  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  #loops through all treatments and creates tx_weight_cutoff variable: replaces any weights above the 90th quantile with the 90th quantile value

  all_trunc_weights=data.frame()

 data_all=data_for_model_with_weights


  # for (k in 1:m){ #cycling through imputations
  dat_w_t<- lapply(all_cutoffs, function(c){
    lapply(1:length(dat_w), function(t){
      cutoff=c
      k=t
      data=dat_w[[t]]
      # browser()

      name=paste0(cutoff, "_weight_cutoff") #labeling new weight column
      # name=paste0(cutoff, "_weight_cutoff") #labeling new weight column
      cutoff_weights=WeightIt::trim(data[,paste0("weights")], cutoff, lower=F)

      # if(cutoff==weights_percentile_cutoff){ #prints only values for user-specified cutoff
        cat(paste0("For ", exposure, "-", outcome, ", using the user-specified cutoff percentile of ", cutoff, ", the median weight is ", round(median(as.numeric(unlist(cutoff_weights))),2) ,
                   " (SD= ",round(sd(as.numeric(unlist(cutoff_weights))),2), "; range= ", round(min(as.numeric(unlist(cutoff_weights))),2), " - ", round(max(as.numeric(unlist(cutoff_weights))),2), ")"), "\n")
      # }
        cat("\n")

      cutoff_weights=data.frame(x=cutoff_weights)
      colnames(cutoff_weights)=name
      # browser()

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
      # data=cbind(data, cutoff_weights)

      # truncated_weights_exp=merge(truncated_weights_exp, new, by=c(colnames(truncated_weights_exp)[colnames(truncated_weights_exp) %in% colnames(new)]), all.x=T)
      # truncated_weights_d<-merge(data, truncated_weights_exp%>%dplyr::filter(.imp==k), by=c(ID), all.x=T)
      # data_for_model_with_weights_cutoff[[paste("fit_", k, "_", exposure, "-", outcome, "_", cutoff, sep="")]] <-truncated_weights_d

      # browser()
      #print histogram of new weights by cutoff value
      ggplot2::ggplot(data=as.data.frame(data), aes(x = as.numeric(data[,name]))) +
        geom_histogram(color = 'black', bins = 15)
      ggplot2::ggsave(paste("Hist_", exposure, "-", outcome, "_weights_cutoff_", cutoff, "imp_", k, ".png", sep=""), path=paste0(home_dir, "final weights/histograms/"), height=8, width=14)

      cat(paste0("A histogram with weights truncated at ", cutoff, " for ", "imp_", k, " have now been saved in the 'final weights/histograms/' folder."), "\n")
      write.csv(data,paste0(home_dir, exposure, "-", outcome, "_weights_cutoff_", cutoff, "imp_", k,"_", paste(all_cutoffs, sep=",", collapse=" "), ".csv") )

      data
    })
  })

  names(dat_w_t)=all_cutoffs
  # dat_w_t=lapply(dat_w_t,function(x){ names(x) <-all_cutoffs})




  #   } #ends k loop
  #
  #   all_trunc_weights=rbind(all_trunc_weights, truncated_weights_exp)
  #
  # } #ends cutoff loop

  #adds weights to mids object
  # colnames(all_weights)=c(ID, "weights")

  # a=MatchThem::complete(dat_w_t, "all")
  # a=mice::complete(imp_data_w, action="broad", include = TRUE)
  # a$outcome<-ifelse(a$WAVE==outcome_time_pt, a[,outcome], NA) #makes outcome variable for modeling w/ long dataset later
  # b=merge(a, all_trunc_weights, by=c(ID, ".imp"), all = T)
  # imp_data_w_t=mice::as.mids(b, .imp = ".imp")

  saveRDS(dat_w_t, paste0(home_dir, "final weights/values/",  exposure, "-", outcome,"_dat_w_t.rds"))
  # saveRDS(data_for_model_with_weights_cutoff, paste0(home_dir, "final weights/values/data_for_model_with_weights_cutoff.rds"))

  cat("USER ALERT: final cutoff weights using the user-specified cutoff values and 2 other values for subsequent sensiivity analyses have now each been saved as a dataset in 'final weights' folder","\n")
  cat(paste0("A histogram of weights for exposure ", exposure, " on ", outcome, ", all imputations has been saved to the 'final weights/histograms/' folder"),"\n")
  cat("\n")
  return(dat_w_t)


}
