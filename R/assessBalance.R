#' Assess success of covariate balancing
#'
#' #Assesses how well balance was achieved for each of the covariates/potential confounds in relation to each of the exposures,
#' and returns a list of unbalanced covariates for each exposure to add to future models
#'
#' @param home_dir path to home directory for the project
#' @param weights_models list of weights models from CBPS
#' @param m number of imputations from Amelia
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @param balance_thresh correlation value between covariate and exposure over which it is not considered balanced
#' @param just_made_weights="no" no=weights are not assigned in global environment and are instead saved locally
#' @return list of unbalanced_covariates_for_models for each exposure
#' @seealso [CBPS::balance()] for more on the balance function
#' @seealso [createWeights()] for more on the weights_models param
#' @export
#' @importFrom CBPS balance
#' @examples assessBalance(home_dir, weights_models=list(), m, exposures, time_pts, balance_thresh=0.12, just_made_weights="no")
assessBalance <- function (home_dir, weights_models=list(), m, exposures, outcomes, time_pts, balance_thresh=0.12, just_made_weights="no"){

  #read in weights saved locally if user has not just made them and they are not yet in global environment
  if (just_made_weights=="no"){

    weights_models=readRDS(paste0(home_dir, "weight fits/weights_models.rds"))

  }

  #assessing balance and determining unbalanced covariates (i.e., those still correlated above 0.12)
  unbalanced_variables=list()
  all_post_balance_corrs=list()

  #iterating through saved weights models ("fits")
  for (f in 1:length(weights_models)){

    balance=CBPS::balance(weights_models[f][[1]])
    balance=as.data.frame(balance$balanced)
    balance=as.data.frame(balance)
    balance=data.table::setDT(balance, keep.rownames = TRUE)[]

    all_post_balance_corrs[[paste("unbalanced_vars_", names(weights_models)[f], sep="")]] <-balance

    unbalanced_vars=balance[abs(balance$`Pearson Correlation`)>balance_thresh,] #0.12 corr default threshold

    write.csv(x=as.data.frame(unbalanced_vars), file=paste(home_dir, "balance/unbalanced_vars_", names(weights_models)[f], ".csv", sep=""))
    print(paste0("A list of any unbalanced variables for ", names(weights_models)[f] ," has now been added to the 'balance' folder as a csv file"))

    unbalanced_variables[[paste("unbalanced_vars_", names(weights_models)[f], sep="")]] <- unbalanced_vars
  }

  #averages across imputed datasets to determine final variables that remain correlated with exposure at r>0.12
  sig_unbalanced_averaged_k=list()
  unbalanced_covariates_for_models=list()

  for (z in seq(length(outcomes))){

    for (x in 1:length(exposures)){

      significant_corrs_remaining={}
      for (y in 1:length(time_pts)){

        #gathering all imputed data set lists for a given exposure, outcome, and time point
        exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposures[x], "-", outcomes[z], "_", time_pts[y]), names(all_post_balance_corrs))]]
        exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[sapply(strsplit(names(all_post_balance_corrs), "_"),"[",6)==time_pts[y] &
                                                                             sapply(strsplit(names(all_post_balance_corrs), "_"),"[",5)==paste0(exposures[x], "-", outcomes[z])]]
        # exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposures[x], "_", time_pts[time_pts!=time_pts[y]]), names(all_post_balance_corrs))==F]] #making sure other time pts not included for nested cases (e.g., 15 and 154)

        exposure_time_corrs=do.call(cbind, exposure_list) #cbinds them
        names=exposure_time_corrs[,1] #extracts variable names
        exposure_time_corrs=as.data.frame((exposure_time_corrs))
        exposure_time_corrs=exposure_time_corrs[,unlist(lapply(exposure_time_corrs, is.numeric), use.names = FALSE)] #extracts numeric columns (i.e., correlations)
        exposure_time_corrs$mean_cor=rowMeans(exposure_time_corrs) #takes row means to find average correlation
        exposure_time_corrs=cbind(names, exposure_time_corrs) #adds names back

        #saves out all correlations
        write.csv(exposure_time_corrs, file=paste(home_dir, "final weights/all_post_balance_cors_", exposures[x],"-", outcomes[z], "_", time_pts[y], ".csv", sep=""))
        print(paste0("Check 'weights' folder for a csv file for ",exposures[x],"-", outcomes[z], "_", time_pts[y] ," of all correlations post-balancing and averages"))

        #writes out only remaining mean correlations over 0.12--the ones left unbalanced
        sig=exposure_time_corrs[abs(exposure_time_corrs$mean_cor)>balance_thresh,]
        colnames(sig)= c("covariate", apply(expand.grid("corr_imp", 1:m), 1, paste, sep="", collapse="_"), "mean_corr")
        # assign(paste("sig_unbalanced_vars_",  exposures[x], "_", time_pts[y], sep=""), sig)
        sig_unbalanced_averaged_k[[paste("sig_unbalanced_vars_",  exposures[x], "-", outcomes[z],"_", time_pts[y], sep="")]] <- sig

        write.csv(sig, file=paste(home_dir, "final weights/all_sig_post_balance_cors_", exposures[x],"-", outcomes[z], "_", time_pts[y], ".csv", sep=""))
        print(paste0("USER ALERT: Check 'final weights' folder for a csv file of all SIGNIFICANT correlations above ", as.character(balance_thresh), " for ",  exposures[x], "-", outcomes[z],"_", time_pts[y]," post-balancing and averages"))

        significant_corrs_remaining=rbind(significant_corrs_remaining, sig)

      }

      #compiles all remaining correlations over 0.12 for each treatment (including all time points)
      print(paste0("USER ALERT: Inspect the following list of unbalanced covariates for exposure ", exposures[x], "-", outcomes[z], " :"))
      print(significant_corrs_remaining)

      write.csv(significant_corrs_remaining, file=paste(home_dir, "final weights/all_sig_post_balance_cors_", exposures[x], "-", outcomes[z],".csv", sep=""))
      print(paste0("All covariates still significantly correlated with ", exposures[x], " have been saved as a csv file in 'final weights"))

      covariates_for_model=unique(c(as.character(unlist(significant_corrs_remaining[,1]))))
      covariates_for_model=covariates_for_model[!grepl(exposures[x], covariates_for_model)] #removing exposure (as this will be modeled explicitly)
      covariates_for_model=covariates_for_model[! covariates_for_model %in% exposures[x]]
      covariates_for_model=paste0(covariates_for_model, sep="", collapse=" + ")
      covariates_for_model=gsub(".", "_", covariates_for_model, fixed=TRUE)

      unbalanced_covariates_for_models[[paste0(exposures[x], "-", outcomes[z])]] <-covariates_for_model

    }
  }

  return(unbalanced_covariates_for_models)

}

#' Get weights
#' Code to get and assign weights output if code has already been run to save weights locally
#' @param home_dir path to home directory for the project
#' @param just_made_weights="no" no=weights are not assigned in global environment and are instead saved locally
#' @return weights_models
#' @export
#' @examples getWeights(home_dir, just_made_weights="no")
#' @seealso [assessBalance()] for more on the weights_models param

getWeights <- function(home_dir, just_made_weights="no"){

  if  (just_made_weights=="no"){
    weights_models=readRDS(paste0(home_dir, "weight fits/weights_models.rds"))

    return(weights_models)
  }
}

