#' Creates formulas for generating weights
#'
#' Creates full formulas for each exposure at each time point including all time-varying covariates at all prior lags
#'
#' @param object msm object that contains all relevant user inputs
#' @param wide_long_datasets from formatForWeights
#' @param all_potential_covariates from selectCoariates
#' @return forms
#' @export
#' @seealso [formatForWeights()], [identifyPotentialConfounds] for more on inputs
#' @examples createForms(object, wide_long_datasets, all_potential_covariates)
#'
createForms <- function(object, wide_long_datasets, all_potential_covariates){

  ID=object$ID
  home_dir=object$home_dir
  exposure=object$exposure
  outcome=object$outcome
  outcome_time_pt=object$outcome_time_pt
  time_pts=object$time_pts
  time_var_exclude=object$time_var_exclude
  exposure_time_pts=object$exposure_time_pts
  time_varying_covariates=object$time_varying_variables
  keep_concurrent_tv_vars=object$keep_concurrent_tv_vars


  cat("USER ALERT: Please inspect the full balancing formulas for each exposure time point below. They contain all confounders that will be used to assess balance at each exposure time point.", "\n",
      "By default, these formulas exclude any contemporaneous time-varying variables (given that they are difficult to disentangle from mediators).","\n",
      "If there aare contemporaneous time-varying variables that you wish to include in the balancing formula, please list them in the 'keep_concurrent_tv_vars' field of the msmObject and re-run this function.","\n")
  cat("\n")
  cat("All full balancing formulas have been saved out in the 'forms/' folder as .csv files", "\n")
  cat("\n")

  covariates_to_include=all_potential_covariates
  forms=list()
  forms_csv=data.frame()

  #Cycles through all exposure time points
  for (x in 1:length(exposure_time_pts)){
    #finds variables to include for that time point
    time=exposure_time_pts[x]
    vars_to_include=c(all_potential_covariates[as.numeric(sapply(strsplit(all_potential_covariates, "\\."), "[", 2))<time #all lagged tv covars only
                                               | is.na(sapply(strsplit(all_potential_covariates, "\\."), "[", 2))]) #inc baseline covars
    if(x>1){ #adds lagged exposure
      vars_to_include=c(vars_to_include,apply(expand.grid(exposure, as.character(exposure_time_pts[exposure_time_pts<time])), 1, paste, sep="", collapse="."))
    }

    vars_to_include=vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                            time_varying_covariates, #exclude time-varying covariates in long form (already in wide form)
                                                            paste(exposure, time, sep="."), #excludes exposure from concurrent time (it is the DV of the balancing form)
                                                            apply(expand.grid(time_varying_covariates[!time_varying_covariates %in% keep_concurrent_tv_vars], as.character(time)), 1, paste, sep="", collapse="."),#exclude any time-varying confounders current time point --cannot distinguish from mediators!
                                                            # apply(expand.grid(outcome, as.character(time_pts)), 1, paste, sep="", collapse="."), #exclude outcome at current time point
                                                            paste(outcome, time, sep="."), #exclude outcome at current time point
                                                            apply(expand.grid(outcome, as.character(time_pts[time_pts>time])), 1, paste, sep="", collapse="."),#excludes future outcome
                                                            time_var_exclude #exclude any time points that should not be there bc of planned missingness
    )]
    vars_to_include=vars_to_include[!duplicated(vars_to_include)]

    #Creates form for the given exposure time point
    f=as.formula(paste(paste0(exposure,".", time, " ~ "), paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + ")))

    #prints form for user inspection
    cat(paste0("The full formula for ", exposure, "-", outcome, " at ", exposure," time point ", as.character(time),
               " is:"), "\n")
    print(f)
    cat("\n")

    #writes to csv file
    forms_csv_temp=data.frame()
    forms_csv_temp[1,1]=paste0("Full formula for ", exposure, "-", outcome, " at ", exposure," time point ", as.character(time),":")
    forms_csv_temp[1,2]=paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + "))
    forms_csv=rbind(forms_csv, forms_csv_temp)

    #assigns to forms
    forms[[paste("form_", exposure,"-", outcome, "-", time, sep="")]] <- f
  }

  #tallies total number of covariates across all forms
  cat(paste0("Across all full balancing formulas at all exposure time points, there are a total of ",
             sum(unlist(lapply(1:nrow(forms_csv), function(x){length(unlist(strsplit(forms_csv[x,2], "\\+")))}))),
             " covariate confounders."))
  cat("\n")

  write.csv(forms_csv, paste0(home_dir, "forms/",  exposure, "-", outcome,"_full_balancing_formulas.csv", sep=""), row.names = F)
  cat("All full balancing formulas have been saved in the 'forms/' folder.")
  cat("\n")

  return(forms)
}


