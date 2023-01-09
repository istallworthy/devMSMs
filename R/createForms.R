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
createForms <- function(object, wide_long_datasets, covariates_to_include){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  outcome_time_pt=object$outcome_time_pt
  time_pts=object$time_pts
  time_var_exclude=object$time_var_exclude
  exposure_time_pts=object$exposure_time_pts
  potential_colliders=object$potential_colliders
  mandatory_keep_covariates=object$mandatory_keep_covariates
  exclude_covariates=object$exclude_covariates
  time_varying_covariates=object$time_varying_variables
  keep_concurrent_tv_vars=object$keep_concurrent_tv_vars

  #error checking
  if (length(potential_colliders$colliders)!=length(potential_colliders$exp_out_pair)){
    stop('Please provide at least one collider variable for each listed exposure-outcome pair in the msm object')}
  if (length(potential_colliders$exclude_lags)!=length(potential_colliders$exp_out_pair)){
    stop('Please indicate whether you want to exclude lagged values of colliders for each listed exposure-outcome pair in the msm object')}


  # browser()
  ####Creates CBPS forms for each exposure at each time point (e.g., HOMEETA.6) that includes: all identified potential confounds for that treatment (at any time point) and lagged time points of the treatment, excludes other outcomes a the given time point (as potential colliders).

  #determining column names
  imp_col_names=colnames(wide_long_datasets[[1]])

  forms=list()
  forms_csv=data.frame()

  cat("USER ALERT: Please inspect the balancing weights formulas below that contain all confounders that will be used to create weights at each exposure time point.", "\n",
  "These formula will reflect the associations shown in the potential_confound_correlations.html table excluding any contemporaneous time-varying variables (given that they are difficult to disentangle from mediators).","\n",
  "If there aare contemporaneous time-varying variables that you wish to include in the balancing formula, please list them in the 'keep_concurrent_tv_vars' field of the msmObject and re-run this function.","\n",
  "If there are any time-varying variables you wish to omit given their potential to be colliders, please list them in the 'potential_colliders' field of the msmObject and re-run,","\n",
  "You can also omit or include variables in all balancing forms using the 'exclude_covariates' and 'mandatory_keep_covariates' fields of the msmObject, respectively.", "\n")
  cat("All forms are saved out in the 'forms' folder as csv files", "\n")
  cat("\n")

  #cycles through outcomes
  for (z in seq(length(outcomes))){
    outcome=outcomes[z]

    #Cycles through all exposures
    for (y in 1:length(exposures)){
      exposure=exposures[y]

      #gathers all confounders for exposure-outcome pair at all time points
      exposure_covariates=covariates_to_include[[paste0(exposure, "-", outcome, sep="")]]

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
      #gathers potential colliders
      #checking to see if there are colliders listed for exposure-outcome pair
      if (nrow(potential_colliders)>0 && sum(potential_colliders$exp_out_pair== paste0(exposure, "-", outcome))>0){
        #gathers all colliders
        colliders=c(as.character(unlist(potential_colliders$colliders[potential_colliders$exp_out_pair== paste0(exposure, "-", outcome)])))
        #finds if user wants to included lags of colliders
        colliders_lagged=potential_colliders$exclude_lags[potential_colliders$exp_out_pair== paste0(exposure, "-", outcome)]
        #finds out if user specified time pts for any colliders
        time_spec_colliders=colliders[grepl("\\.", colliders)]
        #finds any time-varying colliders with not time point specified
        all_time_colliders=colliders[colliders %in% time_varying_covariates]
        #finds any time invariant colliders
        time_invariant_colliders=colliders[!colliders %in% all_time_colliders & !colliders %in% time_spec_colliders]

      }else{
        colliders=NULL
        all_time_colliders=NULL
        colliders_lagged=F
        time_varying_colliders=NULL
        time_spec_colliders=NULL
        time_invariant_colliders=NULL
      }


      #Cycles through all exposure time points
      for (x in 1:length(exposure_time_pts)){
        time=exposure_time_pts[x]

        # browser()

        #finds concurrent time-varying and time-invariant potential confounds
        concurrent_covariates=na.omit(exposure_covariates[exposure_covariates$exp_time==time,])
        concurrent_covariates=c(concurrent_covariates$row, concurrent_covariates$column)
        concurrent_covariates=unique(concurrent_covariates)

        #Finds relevant lagged time points
        lags=c(exposure_time_pts[exposure_time_pts<time])
        lagged_vars={} #variable for lagged variables
        if (length(lags)==0){ #no lags
          # vars_to_include=c(time_varying_concurrent_covariates, time_invariant_concurrent_covariates) #retains only concurrent variables
          vars_to_include=c(concurrent_covariates) #retains only concurrent variables

        }else{ #gets lagged exposures
          past_exposures=apply(expand.grid(exposure, as.character(lags)), 1, paste, sep="", collapse=".") #finds past exposure variables
          # keep_time_var_covars=apply(expand.grid(mandatory_keep_covariates[mandatory_keep_covariates %in% time_varying_covariates], as.character(lags)), 1, paste, sep="", collapse=".") #finds past time-varying mandatory covariates

          #populates data frames
          lagged_covariates=past_exposures
          # vars_to_include=c(time_varying_concurrent_covariates, time_invariant_concurrent_covariates, lagged_covariates)
          vars_to_include=c(concurrent_covariates, lagged_covariates)

        }



        # #EXCLUSION STAGE


        #gathering all concurrent colliders
        concurrent_colliders={}
        lagged_colliders={}
        #finding all concurrent and lagged time-specified colliders
        if (length(time_spec_colliders>0)){
          if (sum(as.numeric(sapply(strsplit(time_spec_colliders, "\\."), "[", 2))==time)>0){ #add any time-speciic colliders to concurrent if time matches
            concurrent_colliders=c(concurrent_colliders, time_spec_colliders[as.numeric(sapply(strsplit(time_spec_colliders, "\\."), "[", 2))==time])
          }
          if(sum(sapply(strsplit(time_spec_colliders, "\\."), "[", 2) %in% lags)>0){ #add time-specific colliders to lagged list if time is lagged
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



        #USER: Place where you can manually delete variables that should not be included in forms such as index variables and possible colliders (e.g., outcomes that could also cause a given tx )
        vars_to_include=vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                                time_varying_covariates, #exclude time-varying covariates in long form (already in wide form)
                                                                paste(exposure, time, sep="."), #excludes exposure from concurrent time (it is the DV of the balancing form)
                                                                apply(expand.grid(time_varying_covariates[!time_varying_covariates %in% keep_concurrent_tv_vars], as.character(time)), 1, paste, sep="", collapse="."),#exclude any time-varying confounders current time point --cannot distinguish from mediators!
                                                                # apply(expand.grid(outcome, as.character(time_pts)), 1, paste, sep="", collapse="."), #exclude outcome at current time point
                                                                paste(outcome, time, sep="."), #exclude outcome at current time point
                                                                apply(expand.grid(outcome, as.character(time_pts[time_pts>time])), 1, paste, sep="", collapse="."),#excludes future outcomes
                                                                exclude_covariates[grepl(paste(exposure, time, sep="."), exclude_covariates)==F], #exclude any covariates specified by the user
                                                                exclude_covariates[!exclude_covariates %in% exposure], #exclude any covariates specified by the user (not counting exposure)
                                                                time_var_exclude, #exclude any time points that should not be there bc of planned missingness
                                                                concurrent_colliders, #exclude colliders specified at exposure time pt
                                                                # apply(expand.grid(colliders, as.character(outcome_time_pt)), 1, paste, sep="", collapse="."), #exclude colliders at outcome time points
                                                                lagged_colliders #any lagged colliders as indicated
        )]

        vars_to_include=vars_to_include[!duplicated(vars_to_include)]

        #Creates form for the given exposure time point
        f=as.formula(paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + ")))

        #prints form for user inspection
        cat(paste0("Formula for ", exposure, "-", outcome, " at ", exposure," time point ", as.character(time),
                   ":"), "\n")
        print(f)
        cat("\n")

        forms_csv_temp=data.frame()
        forms_csv_temp[1,1]=paste0("Formula for ", exposure, "-", outcome, " at ", exposure," time point ", as.character(time),":")
        forms_csv_temp[1,2]=paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + "))

        forms_csv=rbind(forms_csv, forms_csv_temp)

        #saves out form to global workspace and local directory
        # return(assign(paste("form_", exposure, "_", time, sep=""), f))
        forms[[paste("form_", exposure,"-", outcome, "-", time, sep="")]] <- f


      }
    }
  }
  write.csv(forms_csv, paste0(home_dir, "forms/all_balancing_formulas.csv", sep=""), row.names = F)

  return(forms)
}


