#' Creating updated short forms with any imbalanced covariates at time lags greater than the user-specified lag
#' @param object msm object that contains all relevant user inputs
#' @param forms short forms from createShortForms
#' @param data_for_model_with_weights all imputed datasets with weights
#' @param balance_stats_full balance assessment using full forms
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @return new_forms
#' @export
#' @examples updateForms(object, forms, data_for_model_with_weights, balance_stats_full)
updateForms <-function(object, forms, data_for_model_with_weights, balance_stats_full){

  exposure=object$exposure
  exposure_time_pts=object$exposure_time_pts
  outcome=object$outcome
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates
  home_dir=object$home_dir
  short_form_lag=object$short_form_lag

  forms_csv=data.frame()

  cat(paste0("Short balancing formulas at each exposure time point will be updated to include any time-varying covariates at lags greater than ",
             short_form_lag, " that remain imbalanced."), "\n")
  cat("\n")

  # forms=short_forms
  bal_stats=balance_stats_full
  new_forms=forms

  #averaging bal stats across all imputed datasets to determine residual imbalance
  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(bal_stats, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=bal_stats[[1]]$exp_time,
                               covar_time=bal_stats[[1]]$covar_time,
                               covariate=bal_stats[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0))%>%
    dplyr::filter(balanced_avg==0)

  #cycling thru time pts w/ lag >t-1
  new_forms=lapply(seq(1:length(exposure_time_pts)), function(i){
    exposure_time_pt=exposure_time_pts[i]

    f=forms[[names(forms)[grepl(paste0("form_", exposure, "-", outcome, "-", exposure_time_pt), names(forms))]]]

    if (i>2){ #ignore first 2 time points
      temp=unbalanced_covars%>%dplyr::filter(exp_time==exposure_time_pt, as.numeric(covar_time)<exposure_time_pts[i-1],
                                             as.numeric(covar_time)>0)%>% #finds any lagged imbalanced covars
        dplyr::select(covariate)

      #renames factors (that were appended w/ level)
      if (nrow(temp)>0){
      temp$covariate[sapply(strsplit(sapply(strsplit(temp$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(temp$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(temp$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]
      temp=as.character(unlist(temp))

      cat(paste0("For ", exposure, " at exposure time point ", exposure_time_pt, ", the following covariate(s) will be added to the short balancing formula: "), temp, "\n")
      f=paste( paste(deparse(f, width.cutoff = 500), collapse=""), paste(temp, sep="", collapse=" + "), sep=" + ")
      }
    }

    f=as.formula(f)
    f
  })

  names(new_forms)<-names(forms)

  forms_csv=data.frame(name=names(lapply(new_forms, function(f){paste(deparse(f, width.cutoff = 500), collapse="")})),
                        form=unlist(lapply(new_forms, function(f){paste(deparse(f, width.cutoff = 500), collapse="")})))
  cat("\n")
  cat(paste0("Across all updated balancing formulas at all exposure time points, there are a total of ",
             sum(unlist(lapply(1:nrow(forms_csv), function(x){length(unlist(strsplit(forms_csv[x,2], "\\+")))}))),
             " covariate confounders."))
  cat("\n")
  write.csv(forms_csv, paste0(home_dir, "forms/",  exposure, "-", outcome, "_updated_short_balancing_formulas.csv", sep=""), row.names = F)
  cat("Updated short forms have now been saved in the 'forms/' folder.")

  return(new_forms)
}
