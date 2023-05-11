

updateForms <-function(object, forms, data_for_model_with_weights, balance_stats_full){

  exposure=object$exposure
  exposure_time_pts=object$exposure_time_pts
  outcome=object$outcome
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates
  home_dir=object$home_dir

  forms_csv=data.frame()


  # forms=short_forms
  bal_stats=balance_stats_full

  new_forms=forms

  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(bal_stats, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=bal_stats[[1]]$exp_time,
                               covar_time=bal_stats[[1]]$covar_time,
                               covariate=bal_stats[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0))%>%
    dplyr::filter(balanced_avg==0)


  # browser()
  #cycling thru time pts w/ lag >t-1
  new_forms=lapply(seq(1:length(exposure_time_pts)), function(i){
    exposure_time_pt=exposure_time_pts[i]

    f=forms[[names(forms)[grepl(paste0("form_", exposure, "-", outcome, "-", exposure_time_pt), names(forms))]]]

    if (i>2){

      temp=unbalanced_covars%>%dplyr::filter(exp_time==exposure_time_pt, as.numeric(covar_time)<exposure_time_pts[i-1],
                                             as.numeric(covar_time)>0)%>% #finds any lagged imbalanced covars
        dplyr::select(covariate)

      #renames factors (that were appended w/ level)
      # browser()
      if (nrow(temp)>0){
      temp$covariate[sapply(strsplit(sapply(strsplit(temp$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(temp$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(temp$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]

      temp=as.character(unlist(temp))

      cat(paste0("For ", exposure, " at t ", exposure_time_pt, ", the following covariate(s) will be added to the short form: "), temp, "\n")
      f=paste( paste(deparse(f, width.cutoff = 500), collapse=""), paste(temp, sep="", collapse=" + "), sep=" + ")
      # forms_csv_temp[1,2]=paste( paste(deparse(f, width.cutoff = 500), collapse=""), paste(temp, sep="", collapse=" + "), sep=" + ")
      }


    }

    f=as.formula(f)
    f


  })

  names(new_forms)  <-names(forms)

  # browser()
  forms_csv=data.frame(name=names(lapply(new_forms, function(f){paste(deparse(f, width.cutoff = 500), collapse="")})),
                        form=unlist(lapply(new_forms, function(f){paste(deparse(f, width.cutoff = 500), collapse="")})))

  write.csv(forms_csv, paste0(home_dir, "forms/",  exposure, "-", outcome, "_updated_short_balancing_formulas.csv", sep=""), row.names = F)

  return(new_forms)
}
