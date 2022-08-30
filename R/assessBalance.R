

#' Title
#'
#' @return
#' @export
#' @importFrom CBPS balance

#' @examples
assessBalance <- function (home_dir, exposures, time_pts, just_made_weights="no"){

  #read in weights saved locally if user has not just made them and they are not yet in global environment
  if (just_made_weights="no"){
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


  ####Creates variables in global environment of all post-balance correlations and variables that remain unbalanced after CBPS
  fits <- as.character(ls(pat="fit"))
  fits=fits[!grepl("balanced", fits)]
  for (f in 2:(length(fits)-1)){
    fit=get(fits[f])
    b=CBPS::balance(fit)
    b=b$balanced
    b=as.data.frame(b)
    b=setDT(b, keep.rownames = TRUE)[]
    all_post_bal_vars=b
    unbalanced_vars=b[abs(b$`Pearson Correlation`)>0.12,] #IS changed 10/13/20
    assign(paste("unbalanced_vars_", fits[f], sep=""), unbalanced_vars)
    assign(paste("all_post_balanced_vars_",fits[f], sep=""), all_post_bal_vars)
    unbalanced_vars=NULL
  }

  #averages across imputed datasets to determine final variables that remain correlated with exposure at r>0.12
  for (x in 1:length(exposures)){
    tx_all_t={}
    for (y in 1:length(time_pts)){
      #pulls in each imputed dataste for each treatment and time pt
      a=get(paste("all_post_balanced_vars_fit_1_",exposures[x], "_", time_pts[y], sep=""))
      colnames(a)=c("rn", "A")
      b=get(paste("all_post_balanced_vars_fit_2_",exposures[x], "_", time_pts[y], sep=""))
      colnames(b)=c("rn", "B")
      c=get(paste("all_post_balanced_vars_fit_3_",exposures[x], "_", time_pts[y], sep=""))
      colnames(c)=c("rn", "C")
      d=get(paste("all_post_balanced_vars_fit_4_",exposures[x], "_", time_pts[y], sep=""))
      colnames(d)=c("rn", "D")
      e=get(paste("all_post_balanced_vars_fit_5_",exposures[x], "_", time_pts[y], sep=""))
      colnames(e)=c("rn", "E")

      test=merge(a,b, by="rn")
      test=merge(test,c, by="rn")
      test=merge(test,d, by="rn")
      test=merge(test,e, by="rn")

      test$mean_cor=rowMeans(test[,2:6])

      #saves out all correlations
      write.csv(test, file=paste(home_dir, "weights/all_post_balance_cors_", exposures[x], "_", time_pts[y], ".csv", sep=""))

      #writes out only remaining correlations over 0.12
      sig=test[abs(test$mean_cor)>0.12,]
      assign(paste("sig_unbalanced_vars_",  exposures[x], "_", time_pts[y], sep=""), sig)
      write.csv(sig, file=paste(home_dir, "weights/all_sig_post_balance_cors_", exposures[x], "_", time_pts[y], ".csv", sep=""))

      tx_all_t=rbind(tx_all_t, sig)
      sig=NULL

    }

    #compiles all remaining correlations over 0.12 for each treatment (including all time points)
    assign(paste("sig_unbalanced_vars_",  exposures[x], sep=""), tx_all_t)

    print("Inspect the following list of unbalanced covariates for exposure ", exposures[x], " :")
    print(tx_all_t)

    write.csv(tx_all_t, file=paste(home_dir, "weights/all_sig_post_balance_cors_", exposures[x], ".csv", sep=""))

  }


  #### Fetching remaining correlated potential confounds to include in final models. Produces forms that can then easily be added to substantive models in future section.
  for (x in 1:length(exposures)){
    covs=get(paste("sig_unbalanced_vars_",  exposures[x], sep=""))$rn
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
