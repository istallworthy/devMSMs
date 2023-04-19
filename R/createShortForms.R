
#creating short forms that contain time invariant covariatse and time varying covariates at only t-1

createShortForms<- function(object, full_forms, keep=NULL){

  home_dir=object$home_dir
  exp_time_pts=object$exposure_time_pts
  short_form_lag=object$short_form_lag
  exposure=object$exposure
  outcome=object$outcome

  if (!is.null(keep) & length(unique(lengths(keep)))!=1){ #makes sure all keep fields are equal
    stop('Make sure the number of entries in each field of "keep" are equal')}

  short_forms=list()


  list=full_forms[names(full_forms)[grepl(exposure, names(full_forms))]]

  if (length(list)!=length(exp_time_pts)){ #makes sure all expossure time points are there
    stop('Make sure the forms list contains forms for each exposure time point for each exposure')}

  for (z in 1:length(list)){
    time=as.numeric(sapply(strsplit(names(list)[z], "-"), "[", length(unlist(strsplit(names(list)[z], "-")))))

    #find any user-specified covariates to keep
    if (!is.null(keep) & keep$exposure %in% exposure & keep$time %in% time){
      keep_cov=keep$tv_covar[which(keep$exposure==exposure & keep$time==time)]
    }else(keep_cov=NA)

    if (time ==exp_time_pts[1] | time ==exp_time_pts[2]){ #ignore first time point as there are no lagged values and second time pt bc only t-1 exists
      short_forms[[names(list)[z]]] <-list[[z]]

    }else{
      form=list[[z]]
      dv=form[[2]]
      covars=paste(deparse(form[[3]], width.cutoff = 500), collapse="")
      covar_time=sapply(strsplit(unlist(strsplit(as.character(covars), "\\+")), "\\."), "[", 2)
      covars=as.character(unlist(strsplit(covars, "\\+")))

      if(!is.na(keep_cov)){
        new_covars=c(covars[!as.numeric(covar_time)<exp_time_pts[z-1] | is.na(covar_time)], keep_cov)}
      else{ new_covars=covars[!as.numeric(covar_time)<exp_time_pts[z-short_form_lag] | is.na(covar_time)]}

      new_form=as.formula(paste0(dv, "~", paste(new_covars, sep="", collapse="+")))

      cat(paste0("The short form for ", names(list)[z], " including time-varying covariates at t-", short_form_lag, " only is:"), "\n")
      # print(new_form)

      # browser()
      print(new_form)
      cat("\n")

      short_forms[[names(list)[z]]] <- new_form
    }
  }


  saveRDS(short_forms, paste0(home_dir, "forms/short_forms.rds"))
  cat("Short formulas including time-varying covariates at t-", short_form_lag, " only have now been saved in the 'forms/' folder", "\n")

  return(short_forms)
}
