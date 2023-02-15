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
truncateWeights <-function(object, wide_long_datasets, imp_data_w){

  home_dir=object$home_dir
  exposures=object$exposures
  exposure_epochs=object$exposure_epochs
  outcomes=object$outcomes
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))
  m=object$m

  if (weights_percentile_cutoff>1 | weights_percentile_cutoff<0){
    stop('Please select a percentile cutoff between 0 and 1 in the msmObject')}


  #creating list of all cutoff values to cycle through
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  data_for_model_with_weights_cutoff=list()
  # truncated_weights<-data.frame(ID=wide_long_datasets[[paste("imp", 1, "_widelong", sep="")]][,ID])
  # colnames(truncated_weights)<- ID
  ids=as.data.frame(mice::complete(imp_data_w,1)[,ID])
  ids=ids[!duplicated(ids),1]
  truncated_weights=data.frame(expand.grid(ids, 1:m))
  colnames(truncated_weights)<-c(ID,".imp")


  #creating truncated weights at 90th percentile for sensitivity analyses; changed to top coding to avoid exclusion
  #loops through all treatments and creates tx_weight_cutoff variable: replaces any weights above the 90th quantile with the 90th quantile value

  all_trunc_weights=data.frame()


  for (z in seq(length(outcomes))){ #cycling through all outcomes
    outcome=outcomes[z]

    for (t in 1:length(exposures)){ #cycling through all exposures
      exposure=exposures[t]

      data_all=readRDS(paste(paste0(home_dir, "original weights/data_for_model_with_weights.rds", sep="")))

      for (k in 1:m){ #cycling through imputations
        truncated_weights_exp=data.frame(expand.grid(ids,k))
        colnames(truncated_weights_exp)=c(ID, ".imp")


        for (c in 1:length(all_cutoffs)){ #cycling through all cutoff values
          cutoff=all_cutoffs[c]

          # data=complete(imp_data_w,k)
          data=as.data.frame(data_all[[paste("fit_", k, "_", exposure, "-", outcome, sep="")]])

          name=paste0(exposure, "-", outcome,"_", cutoff, "_weight_cutoff") #labeling new weight column
          # cutoff_weights=data_for_model_with_weights%>% #conducting perecentile cutoff and top coding
          #   dplyr::mutate(!!new_var := ifelse(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "-", outcome, "_weight")]<
          #                                       quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "-", outcome, "_weight")], probs=cutoff),
          #                                     data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure,"-", outcome,"_weight")],
          #                                     quantile(data_for_model_with_weights[,colnames(data_for_model_with_weights)==paste0(exposure, "-", outcome, "_weight")], probs=cutoff)
          #   ) )%>%dplyr::select(all_of(new_var))
          cutoff_weights=WeightIt::trim(data[,paste0(exposure, "-", outcome, "_weights")], cutoff, lower=F)

          if(cutoff==weights_percentile_cutoff){ #prints only values for user-specified cutoff
            cat(paste0("For ", exposure, "-", outcome, ", using the user-specified cutoff percentile of ", cutoff, ", the median weight is ", median(as.numeric(unlist(cutoff_weights))) ,
                       " (SD= ",sd(as.numeric(unlist(cutoff_weights))), "; range= ", min(as.numeric(unlist(cutoff_weights))), " - ", max(as.numeric(unlist(cutoff_weights))), ")"), "\n")
          }

          #binding on new cutoff weight

          cutoff_weights=data.frame(x=cutoff_weights)
          colnames(cutoff_weights)=name

          data_cutoff=cbind(data, cutoff_weights)


          new=data.frame(ID=data_cutoff[,ID])
          colnames(new)<-ID
          #calculates the mean value for each exposure for each exposure epoch
          for (e in 1:nrow(exposure_epochs)){
            epoch=exposure_epochs[e,1]
            temp=data.frame(row.names=1:nrow(data_cutoff))
            new_var=paste0(exposure, "_", epoch)

            #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
            for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
              level=as.numeric(unlist(exposure_epochs[e,2]))[l]
              z=data_cutoff[,which(grepl(paste0(exposure, ".", level), names(data_cutoff)))]
              temp=cbind(temp, z)
            }
            #adds a new variable of the exposure averaged within epoch
            new=new%>%dplyr::mutate(.imp=k,
                                    !!new_var :=rowMeans(temp, na.rm=T))
            new[,new_var]=as.numeric(new[,new_var])
          }

          new=cbind(new, cutoff_weights)

          # truncated_weights_exp=rbind(truncated_weights_exp, new)
          truncated_weights_exp=merge(truncated_weights_exp, new, by=c(colnames(truncated_weights_exp)[colnames(truncated_weights_exp) %in% colnames(new)]), all.x=T)

          # test=dplyr::rows_update(truncated_weights_exp, new, by=c(colnames(truncated_weights_exp)[colnames(truncated_weights_exp) %in% colnames(new)]), unmatched = "ignore")

          # test2=rquery::natural_join(truncated_weights_exp, new, by=c(colnames(truncated_weights_exp)[colnames(truncated_weights_exp) %in% colnames(new)]), jointype="FULL")

           # truncated_weights_exp[unlist(truncated_weights_exp$S_ID)==unlist(new$S_ID) &
          #                         truncated_weights_exp$.imp==new$.imp &
          #                         truncated_weights_exp$InRatioCor_Infancy==new$InRatioCor_Infancy, "InRatioCor-EF_avg_perc_0.92_weight_cutoff"]=
          #   new$`InRatioCor-EF_avg_perc_0.92_weight_cutoff`
          # truncated_weights_exp=merge(truncated_weights_exp, new, by=c(ID, ".imp", paste(exposure, exposure_epochs$epochs, sep="_")), all.x=T)


          truncated_weights_d<-merge(data, truncated_weights_exp%>%dplyr::filter(.imp==k), by=c(ID), all.x=T)


          data_for_model_with_weights_cutoff[[paste("fit_", k, "_", exposure, "-", outcome, "_", cutoff, sep="")]] <-truncated_weights_d

          #print histogram of new weights by cutoff value
          ggplot2::ggplot(data=as.data.frame(data_cutoff), aes(x = as.numeric(data_cutoff[,name]))) +
            geom_histogram(color = 'black', bins = 15)
          ggplot2::ggsave(paste("Hist_", exposure, "-", outcome, "_weights_cutoff_", cutoff, "imp_", k, ".png", sep=""), path=paste0(home_dir, "final weights/histograms/"), height=8, width=14)

          cat(paste0("A histogram of weights cut off at ", cutoff, " for exposure ", exposure, " on ", outcome, ", imputation ", k,  " has been saved to the 'final weights/histograms/' folder"),"\n")
          cat("\n")

          write.csv(data_cutoff,paste0(home_dir, exposure, "-", outcome, "_weights_cutoff_", cutoff, "imp_", k,"_", paste(all_cutoffs, sep=",", collapse=" "), ".csv") )
          cat("\n")
          cat("USER ALERT: final cutoff weights using the user-specified cutoff values and 2 other values for subsequent sensiivity analyses have now each been saved as a dataset in 'final weights' folder","\n")


        } #ends k loop

        all_trunc_weights=rbind(all_trunc_weights, truncated_weights_exp)

        # truncated_weights<-merge(truncated_weights, truncated_weights_exp, by=c(ID, ".imp"), all.x=T)

      } #ends exposure loop

    } #ends outcome

  } #ends cutoff loop

  #adds weights to mids object
  # colnames(all_weights)=c(ID, "weights")
  a=complete(imp_data_w, action="long", include = TRUE)
  # a=complete(imp_data_w,2)
  # a$weights=NA
  b=merge(a, all_trunc_weights, by=c(ID, ".imp"), all = T)
  imp_data_w_t=mice::as.mids(b, .imp = ".imp")


  return(data_for_model_with_weights_cutoff)


}
