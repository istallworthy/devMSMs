#' Creates formulas for generating weights
#'
#' Creates formulas for each exposure at each time point for the CBPS function that will make weights returning a list of forms
#'
#' @param object msm object that contains all relevant user inputs
#' @param covariates_to_include from identifyPotentialConfounds
#' @return forms
#' @export
#' @seealso [formatForWeights()], [identifyPotentialConfounds] for more on inputs
#' @examples createForms(object, wide_long_datasets,covariates_to_include, potential_colliders)
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
  potential_colliders=object$potential_colliders
  mandatory_keep_covariates=object$mandatory_keep_covariates
  exclude_covariates=object$exclude_covariates
  time_varying_covariates=object$time_varying_variables
  keep_concurrent_tv_vars=object$keep_concurrent_tv_vars

  covariates_to_include=all_potential_covariates


  # browser()
  ####Creates CBPS forms for each exposure at each time point (e.g., HOMEETA.6) that includes: all identified potential confounds for that treatment (at any time point) and lagged time points of the treatment, excludes other outcome a the given time point (as potential colliders).

  #determining column names
  imp_col_names=colnames(wide_long_datasets[[1]])

  forms=list()
  forms_csv=data.frame()

  cat("USER ALERT: Please inspect the balancing weights formulas below that contain all confounders that will be used to assess balance at each exposure time point.", "\n",
      "By default, these formulas exclude any contemporaneous time-varying variables (given that they are difficult to disentangle from mediators).","\n",
      "If there aare contemporaneous time-varying variables that you wish to include in the balancing formula, please list them in the 'keep_concurrent_tv_vars' field of the msmObject and re-run this function.","\n",
      "If there are any time-varying variables you wish to omit given their potential to be colliders, please list them in the 'potential_colliders' field of the msmObject and re-run,","\n",
      "You can also omit or include variables in all balancing forms using the 'exclude_covariates' and 'mandatory_keep_covariates' fields of the msmObject, respectively.", "\n")
  cat("All forms are saved out in the 'forms' folder as csv files", "\n")
  cat("\n")

  #cycles through outcome
  # for (z in seq(length(outcome))){
  #   outcome=outcome[z]

  # #Cycles through all exposure
  # for (y in 1:length(exposure)){
  #   exposure=exposure[y]

  #gathers all confounders for exposure-outcome pair at all time points
  # exposure_covariates=covariates_to_include[[paste0(exposure, "-", outcome, sep="")]]
  exposure_covariates=all_potential_covariates

  # browser()
  #assesses and reformat use input of exclude_covariates
  user_input_exclude_covariates=exclude_covariates
  exclude_covariates=NULL
  if (length(user_input_exclude_covariates)>0){
    for (c in 1:length(user_input_exclude_covariates)){
      if (grepl("\\.",user_input_exclude_covariates[c])){ #period indicates user specified something time-varying
        exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_exclude_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified
          temp=apply(expand.grid(user_input_exclude_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          exclude_covariates=c(exclude_covariates,temp)
        }else{ #otherwise it is time invariant
          exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c])
        }
      }
    }
  }


  # browser()
  # #gathers potential colliders
  # #checking to see if there are colliders listed for exposure-outcome pair
  # if (nrow(potential_colliders)>0){
  #   #gathers all colliders
  #   colliders=c(as.character(unlist(potential_colliders$colliders)))
  #   #finds if user wants to included lags of colliders
  #   colliders_lagged=potential_colliders$exclude_lags
  #   #finds out if user specified time pts for any colliders
  #   time_spec_colliders=colliders[grepl("\\.", colliders)]
  #   #finds any time-varying colliders with not time point specified
  #   all_time_colliders=colliders[colliders %in% time_varying_covariates]
  #   #finds any time invariant colliders
  #   time_invariant_colliders=colliders[!colliders %in% all_time_colliders & !colliders %in% time_spec_colliders]
  #
  # }else{
    colliders=NULL
    all_time_colliders=NULL
    colliders_lagged=F
    time_varying_colliders=NULL
    time_spec_colliders=NULL
    time_invariant_colliders=NULL
  # }


  #Cycles through all exposure time points
  for (x in 1:length(exposure_time_pts)){
    time=exposure_time_pts[x]

    vars_to_include=c(all_potential_covariates[as.numeric(sapply(strsplit(all_potential_covariates, "\\."), "[", 2))<time #all lagged tv covars
                                               | is.na(sapply(strsplit(all_potential_covariates, "\\."), "[", 2))]) #inc baseline covars
    if(x>1){ #adds lagged exposure
      vars_to_include=c(vars_to_include,apply(expand.grid(exposure, as.character(exposure_time_pts[exposure_time_pts<time])), 1, paste, sep="", collapse="."))
    }


    # #EXCLUSION STAGE
    #gathering all concurrent colliders
    concurrent_colliders={}
    lagged_colliders={}
    #finding all concurrent and lagged time-specified colliders
    if (length(time_spec_colliders)>0){
      if (sum(as.numeric(sapply(strsplit(time_spec_colliders, "\\."), "[", 2))==time)>0){ #add any time-speciic colliders to concurrent if time matches
        concurrent_colliders=c(concurrent_colliders, time_spec_colliders[as.numeric(sapply(strsplit(time_spec_colliders, "\\."), "[", 2))==time])
      }
      if(as.numeric(sapply(strsplit(time_spec_colliders, "\\."), "[", 2))<time){ #add time-specific colliders to lagged list if time is lagged
        lagged_colliders=c(lagged_colliders, time_spec[sapply(strsplit(time_spec_colliders, "\\."), "[", 2) %in% lags])
      }
    }

    if (length(time_invariant_colliders)>0){ #adds any time invariant colliders to concurrent colliders
      concurrent_colliders=c(concurrent_colliders, time_invariant_colliders)
    }

    if (length(all_time_colliders)>0){ #adds all_time colliders to concurrent colldiers
      concurrent_colliders=c(concurrent_colliders, paste(all_time_colliders, time, sep="."))
      lagged_colliders=c(lagged_colliders)
      if (colliders_lagged=="T"){ #applies lags to all_time_colliders if indicated by user
        lagged_colliders=c(lagged_colliders, apply(expand.grid(all_time_colliders, as.character(lags)), 1, paste, sep="", collapse="."))
      }
    }



    #USER: Place where you can manually delete variables that should not be included in forms such as index variables and possible colliders (e.g., outcome that could also cause a given tx )
    vars_to_include=vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                            time_varying_covariates, #exclude time-varying covariates in long form (already in wide form)
                                                            paste(exposure, time, sep="."), #excludes exposure from concurrent time (it is the DV of the balancing form)
                                                            apply(expand.grid(time_varying_covariates[!time_varying_covariates %in% keep_concurrent_tv_vars], as.character(time)), 1, paste, sep="", collapse="."),#exclude any time-varying confounders current time point --cannot distinguish from mediators!
                                                            # apply(expand.grid(outcome, as.character(time_pts)), 1, paste, sep="", collapse="."), #exclude outcome at current time point
                                                            paste(outcome, time, sep="."), #exclude outcome at current time point
                                                            apply(expand.grid(outcome, as.character(time_pts[time_pts>time])), 1, paste, sep="", collapse="."),#excludes future outcome
                                                            exclude_covariates[grepl(paste(exposure, time, sep="."), exclude_covariates)==F], #exclude any covariates specified by the user
                                                            exclude_covariates[!exclude_covariates %in% exposure], #exclude any covariates specified by the user (not counting exposure)
                                                            time_var_exclude, #exclude any time points that should not be there bc of planned missingness
                                                            concurrent_colliders, #exclude colliders specified at exposure time pt
                                                            # apply(expand.grid(colliders, as.character(outcome_time_pt)), 1, paste, sep="", collapse="."), #exclude colliders at outcome time points
                                                            lagged_colliders #any lagged colliders as indicated
    )]

    vars_to_include=vars_to_include[!duplicated(vars_to_include)]

    #Creates form for the given exposure time point
    # f=as.formula(paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + ")))
    f=as.formula(paste(paste0(exposure,".", time, " ~ "), paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + ")))

    #prints form for user inspection
    cat(paste0("The full formula for ", exposure, "-", outcome, " at ", exposure," time point ", as.character(time),
               ":"), "\n")
    print(f)
    cat("\n")

    # browser()

    forms_csv_temp=data.frame()
    forms_csv_temp[1,1]=paste0("Full formula for ", exposure, "-", outcome, " at ", exposure," time point ", as.character(time),":")
    forms_csv_temp[1,2]=paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + "))

    forms_csv=rbind(forms_csv, forms_csv_temp)

    #saves out form to global workspace and local directory
    # return(assign(paste("form_", exposure, "_", time, sep=""), f))
    forms[[paste("form_", exposure,"-", outcome, "-", time, sep="")]] <- f


  }
  #   }
  # # }
  write.csv(forms_csv, paste0(home_dir, "forms/",  exposure, "-", outcome,"_full_balancing_formulas.csv", sep=""), row.names = F)

  return(forms)
}


