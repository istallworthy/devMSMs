#Code to asssess how well balance was achieved for each of the covariates/potential confounds

#' assessBalance
#' @param home_dir path to home directory for the project
#' @param weights_models list of weights models from CBPS
#' @param m number of imputations from Amelia
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @param balance_thresh correlation value between covariate and exposure over which it is not considered balanced
#' @param just_made_weights="no" no=weights are not assigned in global environment and are instead saved locally
#' @return
#' @export
#' @importFrom CBPS balance

#' @examples
assessBalance <- function (home_dir, weights_models=list(), m, exposures, time_pts, balance_thresh=0.12, just_made_weights="no"){

  #read in weights saved locally if user has not just made them and they are not yet in global environment
  if (just_made_weights=="no"){
    ReadRDSFiles <- function(fileDir, envir = .GlobalEnv) {
      p <- ".rds$"
      rds <- list.files(path = fileDir, pattern = p)
      out <- vapply(rds, FUN = function(.x) {
        nm <- sub(pattern = p, replacement = "", x = .x)
        # read data in
        # out <- readRDS(file = x)
        # load in global env
        assign(nm, value = readRDS(file = file.path(fileDir, .x)), envir = envir)
        if (!exists(nm, envir = envir)) return(FALSE)
        TRUE
      }, FUN.VALUE = logical(1L), USE.NAMES = FALSE)
      if (!all(out)) warning("Some `.rds` files not loaded.", call. = FALSE)
      spc <- paste0(rep('*', times = nchar(fileDir) + 1), collapse = "")
      cat("RDS Files loaded from:", fileDir, spc, rds[out], sep = "\n ", spc)
    }

    ReadRDSFiles(paste0(home_dir, "weight fits/"))
  }


  #added by IS 3.21.20 for more formal balance checking (prints out which variables are unbalanced)

  #assessing balance and determining unbalanced covariates (i.e., those still correlated above 0.12)
  unbalanced_variables=list()
  all_post_balance_corrs=list()

  for (f in 1:length(weights_models)){
    # for (f in 2:(length(fits)-1)){
    balance=CBPS::balance(weights_models[[f]])
    balance=as.data.frame(balance$balanced)
    balance=as.data.frame(balance)
    balance=data.table::setDT(balance, keep.rownames = TRUE)[]

    all_post_balance_corrs[[paste("unbalanced_vars_", names(weights_models)[f], sep="")]] <-balance

    unbalanced_vars=balance[abs(balance$`Pearson Correlation`)>balance_thresh,] #0.12 corr default threshold

    write.csv(x=as.data.frame(unbalanced_vars), file=paste(home_dir, "balance/unbalanced_vars_", names(weights_models)[f], ".csv", sep=""))
    message("A list of any unbalanced variables has now been added to the 'balance' folder as a csv file")

    unbalanced_variables[[paste("unbalanced_vars_", names(weights_models)[f], sep="")]] <- unbalanced_vars
  }

  #averages across imputed datasets to determine final variables that remain correlated with exposure at r>0.12
  sig_unbalanced_averaged_k=list()

  for (x in 1:length(exposures)){

    tx_all_t={}
    for (y in 1:length(time_pts)){

      #finds all unbalanced variable list for the given exposure and time point
      exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposure[x], "_", time_pts[y]), names(all_post_balance_corrs))]]
      exposure_time_corrs=do.call(cbind, exposure_list) #cbinds them
      names=exposure_time_corrs[,1] #extracts variable names
      exposure_time_corrs=as.data.frame((exposure_time_corrs))
      exposure_time_corrs=exposure_time_corrs[,unlist(lapply(exposure_time_corrs, is.numeric), use.names = FALSE)] #extracts numeric columns
      exposure_time_corrs$mean_cor=rowMeans(exposure_time_corrs) #takes row means to find average correlation
      exposure_time_corrs=cbind(names, exposure_time_corrs) #adds names back

      #saves out all correlations
      write.csv(exposure_time_corrs, file=paste(home_dir, "weights/all_post_balance_cors_", exposures[x], "_", time_pts[y], ".csv", sep=""))
      message("Check 'weights' folder for a csv file of all correlations with exposure post-balancing and averages")

      #writes out only remaining correlations over 0.12
      sig=exposure_time_corrs[abs(exposure_time_corrs$mean_cor)>balance_thresh,]
      # assign(paste("sig_unbalanced_vars_",  exposures[x], "_", time_pts[y], sep=""), sig)
      sig_unbalanced_averaged_k[[paste("sig_unbalanced_vars_",  exposures[x], "_", time_pts[y], sep="")]] <- sig

      write.csv(sig, file=paste(home_dir, "weights/all_sig_post_balance_cors_", exposures[x], "_", time_pts[y], ".csv", sep=""))
      message(paste0("Check 'weights' folder for a csv file of all SIGNIFICANT correlations above ", as.character(balance_thresh), " with exposure post-balancing and averages"))
      tx_all_t=rbind(tx_all_t, sig)

    }

    #compiles all remaining correlations over 0.12 for each treatment (including all time points)
    # assign(paste("sig_unbalanced_vars_",  exposures[x], sep=""), tx_all_t)

    print("Inspect the following list of unbalanced covariates for exposure ", exposures[x], " :")
    print(tx_all_t)

    write.csv(tx_all_t, file=paste(home_dir, "weights/all_sig_post_balance_cors_", exposures[x], ".csv", sep=""))

    # return(sig_unbalanced_averaged_k)

  }


  #### Fetching remaining correlated potential confounds to include in final models. Produces forms that can then easily be added to substantive models in future section.
  for (x in 1:length(exposures)){
    covs=get(paste("sig_unbalanced_vars_",  exposures[x], sep=""))$rn
    covs=sig_unbalanced_averaged_k[]
    covs=covs[!grepl(c(".24"), covs)] #these list out the older time pts that we want to avoid including bc of blocking
    covs=covs[!grepl(c(".35"), covs)]
    covs=covs[!grepl(c(".15"), covs)]
    covs=covs[!grepl(c("GrosPay"), covs)]
    covs=covs[!grepl(exposures[x], covs)] #removes lagged tx as these will already be in the substantive model

    b=covs
    covs=paste0(covs, sep="", collapse=" + ") #separates them with + for easy addition to models
    covs=gsub(".", "_", covs, fixed=TRUE)

    assign(paste("covs_",  treatments[x], sep=""), covs) #assigns to global environment

  }


}

