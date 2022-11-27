#' Assess success of covariate balancing
#'
#' #Assesses how well balance was achieved for each of the covariates/potential confounds in relation to each of the exposures,
#' and returns a list of unbalanced covariates for each exposure to add to future models
#'
#' @param object msm object that contains all relevant user inputs
#' @param weights_models optional if createWeights has just run or will read in weights models from local storage
#' @param just_made_weights="no" no=weights are not assigned in global environment and are instead saved locally
#' @return list of unbalanced_covariates_for_models for each exposure
#' @seealso [CBPS::balance()] for more on the balance function
#' @seealso [msmObject()] for more on the weights_models param
#' @seealso [createWeights()] for more on the weights_models param
#' @export
#' @importFrom CBPS balance
#' @importFrom knitr kable
#' @examples assessBalance(object, weights_models=list(), just_made_weights="no")
assessBalance <- function (object, weights_models){

  home_dir=object$home_dir
  m=object$m
  exposures=object$exposures
  outcomes=object$outcomes
  exposure_time_pts=object$exposure_time_pts
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates

  #assessing balance and determining unbalanced covariates (i.e., those still correlated above 0.12)
  unbalanced_variables=list()
  all_post_balance_corrs=list()

  #iterating through saved weights models ("fits") to gather all post balance correlations
  for (f in 1:length(weights_models)){

    balance=CBPS::balance(weights_models[f][[1]])
    balance=as.data.frame(balance$balanced)
    balance=as.data.frame(balance)
    balance=data.table::setDT(balance, keep.rownames = TRUE)[]

    all_post_balance_corrs[[paste("unbalanced_vars_", names(weights_models)[f], sep="")]] <-balance

  }




  #averages across imputed datasets to determine final variables that remain correlated with exposure at r>0.12
  sig_unbalanced_averaged_k=list()
  unbalanced_covariates_for_models=list()
  all_post_balance_correlations=list()

  all_corrs_df={}


  for (z in seq(length(outcomes))){

    for (x in 1:length(exposures)){

      significant_corrs_remaining={}
      for (y in 1:length(exposure_time_pts)){

        #gathering all imputed data set lists for a given exposure, outcome, and time point
        exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposures[x], "-", outcomes[z], "_", exposure_time_pts[y]),
                                                                                 names(all_post_balance_corrs))]]
        exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[sapply(strsplit(names(all_post_balance_corrs), "_"),"[",6)==exposure_time_pts[y] &
                                                                             sapply(strsplit(names(all_post_balance_corrs), "_"),"[",5)==paste0(exposures[x], "-", outcomes[z])]]
        # exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposures[x], "_", time_pts[time_pts!=time_pts[y]]), names(all_post_balance_corrs))==F]] #making sure other time pts not included for nested cases (e.g., 15 and 154)

        exposure_time_corrs=do.call(cbind, exposure_list) #cbinds them
        names=exposure_time_corrs[,1] #extracts variable names
        exposure_time_corrs=as.data.frame((exposure_time_corrs))
        exposure_time_corrs=exposure_time_corrs[,unlist(lapply(exposure_time_corrs, is.numeric), use.names = FALSE)] #extracts numeric columns (i.e., correlations)
        exposure_time_corrs$mean_cor=rowMeans(exposure_time_corrs) #takes row means to find average correlation
        exposure_time_corrs=cbind(names, exposure_time_corrs) #adds names back


        #compiles all post-balance correlations
        all_post_balance_correlations[[paste0(exposures[x],"_", exposure_time_pts[y], "-", outcomes[z])]] <- exposure_time_corrs

        #gathers only remaining mean correlations over 0.12--the ones left unbalanced
        sig=exposure_time_corrs[abs(exposure_time_corrs$mean_cor)>balance_thresh,]
        colnames(sig)= c("covariate", apply(expand.grid("corr_imp", 1:m), 1, paste, sep="", collapse="_"), "mean_corr")
        sig_unbalanced_averaged_k[[paste("sig_unbalanced_vars_",  exposures[x], exposure_time_pts[y], "-", outcomes[z],"_", sep="")]] <- sig

        # write.csv(sig, file=paste(home_dir, "balance/post-balance correlation values/all_sig_post_balance_cors_", exposures[x],"_", exposure_time_pts[y],"-", outcomes[z],".csv", sep=""))
        # cat(paste0("USER ALERT: Check 'balance/post-balance correlation values/' folder for a csv file of all SIGNIFICANT correlations above ", as.character(balance_thresh), " for ",  exposures[x], "_", exposure_time_pts[y], "-", outcomes[z]," post-balancing and averages"),"\n")

        significant_corrs_remaining=rbind(significant_corrs_remaining, sig)

      }


      #compiles all remaining correlations over 0.12 for each treatment (including all time points) to print and save
      # browser()
      if (nrow(significant_corrs_remaining)>0){
      cat(paste0("USER ALERT: Inspect the following list of unbalanced covariates for exposure ", exposures[x], "-", outcomes[z], " :"),"\n")
        # browser()
      cat(knitr::kable(significant_corrs_remaining), sep="\n")
      cat("\n")

      # write.csv(significant_corrs_remaining, file=paste(home_dir, "balance/post-balance correlation values/all_sig_post_balance_cors_", exposures[x], "-", outcomes[z],".csv", sep=""))
      # cat(paste0("All covariates still significantly correlated with ", exposures[x],  "-", outcomes[z]," have been saved as a csv file in 'balance/post-balance correlation values/'"),"\n")
      # cat("\n")

      }else {
        cat(paste0("USER ALERT: There are no unbalanced covariates for exposure ", exposures[x], "-", outcomes[z]),"\n")
        # write.csv(significant_corrs_remaining, file=paste(home_dir, "balance/post-balance correlation values/all_sig_post_balance_cors_", exposures[x], "-", outcomes[z],".csv", sep=""))
        cat("\n")

      }



      covariates_for_model=unique(c(as.character(unlist(significant_corrs_remaining[,1]))))
      #renaming any factors (CBPS added a '1') to original names
      covariates_for_model[grepl(paste(factor_covariates, collapse="|"), paste0(covariates_for_model))]=
        substring(covariates_for_model[grepl(paste(factor_covariates, collapse="|"), paste0(covariates_for_model))], 1, nchar(covariates_for_model[grepl(paste(factor_covariates, collapse="|"), paste0(covariates_for_model))])-1)
      covariates_for_model=covariates_for_model[!grepl(exposures[x], covariates_for_model)] #removing exposure (as this will be modeled explicitly)
      covariates_for_model=covariates_for_model[! covariates_for_model %in% exposures[x]]
      covariates_for_model=paste0(covariates_for_model, sep="", collapse=" + ")
      covariates_for_model=gsub(".", "_", covariates_for_model, fixed=TRUE)

      unbalanced_covariates_for_models[[paste0(exposures[x], "-", outcomes[z])]] <-covariates_for_model

    }
  }

  #makes all corrs into data frame to sort and save out for easy viewing
  for (b in 1:length(all_post_balance_correlations)){
    f=all_post_balance_correlations[b]
    corrs=f[[1]]
    colnames(corrs)=c("Covariate", apply(expand.grid("Imp", seq(m)),1, paste0, sep="", collapse="."), "Mean Corr")
    corrs$exposure=sapply(strsplit(names(f), "_"), "[", 1)
    corrs$exposure_time=sapply(strsplit(sapply(strsplit(names(f), "_"), "[", 2), "-"), "[", 1)
    corrs$outcome=sapply(strsplit(names(f), "-"), "[", 2)
    all_corrs_df=rbind(all_corrs_df, corrs)
  }
  all_corrs_df=all_corrs_df[order(abs(all_corrs_df$`Mean Corr`), all_corrs_df$outcome, decreasing = T),]
  #save out correlations
  stargazer::stargazer(all_corrs_df,type="html", digits=2, column.labels = colnames(all_corrs_df),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/post-balance correlation values/all_post-balance_correlations.html"))

  cat(paste0("USER ALERT: Check 'balance/post-balance correlation values/' folder for a table of all post-balance correlations post-balancing and averaged across imputed datasets"),"\n")
  cat("\n")


  test=as.data.frame(unbalanced_covariates_for_models)
  #save out correlations
  stargazer::stargazer(test,type="html", digits=2, column.labels = colnames(test),summary=FALSE, rownames = FALSE, header=FALSE, no.space=F,column.sep.width = "5pt",

                       out=paste0(home_dir, "balance/unbalanced covariates/all_unbalanced_covariates.html"))
  cat(paste0("USER ALERT: Check 'balance/unbalanced covariates/' folder for table of unbalanced covariates for each exposure-outcome pair based on the cutoff correlation of ", as.character(balance_thresh), " averaged across imputed datasets"),"\n")



  return(unbalanced_covariates_for_models)

}


