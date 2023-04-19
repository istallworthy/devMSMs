
#examining initial imbalance for each exposure, each time point, each covariate for each history

preBalanceChecking <- function(object, wide_long_datasets, all_forms){

  exposure=object$exposure
  outcome=object$outcome

  .call <- match.call()
  form_type=.call[[4]]
  form_name=as.character(form_type)


  cat("The following statistics display covariate imbalance for each exposure at each exposure time point prior to weighting, using full formulas that reflect covariates at all lagged time points. Imbalance is first assessed for each exposure history up until the exposure time point and then average across histories.", "\n")

  # forms=all_forms

  #running balance stats function
  bal_summary=calcBalStats(object, wide_long_datasets, all_forms, form_name, exposure, outcome, k=1, weighted=0)

  return(bal_summary)
}

