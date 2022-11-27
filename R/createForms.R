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


  #error checking
  if (length(potential_colliders$colliders)!=length(potential_colliders$exp_out_pair)){
    stop('Please provide at least one collider variable for each listed exposure-outcome pair in the msm object')}
  if (length(potential_colliders$exclude_lags)!=length(potential_colliders$exp_out_pair)){
    stop('Please indicate whether you want to exclude lagged values of colliders for each listed exposure-outcome pair in the msm object')}


  ####Creates CBPS forms for each exposure at each time point (e.g., HOMEETA.6) that includes: all identified potential confounds for that treatment (at any time point) and lagged time points of the treatment, excludes other outcomes a the given time point (as potential colliders).

  #determining column names
  imp_col_names=colnames(wide_long_datasets[[1]])

  forms=list()

  cat("USER ALERT: Please inspect the weights formulas below. Balancing weights will attempt to balance on all of these potential confounding variables. If there are any time-varying variables you wish to omit at this time point, please list them in the 'potential_colliders' field of the msm object and re-run")
  cat("They are also saved out in the 'forms' folder as csv files", "\n")
  cat("\n")

  #cycles through outcomes
  for (z in seq(length(outcomes))){
    outcome=outcomes[z]

    #Cycles through all exposures
    for (y in 1:length(exposures)){
      exposure=exposures[y]

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


      #gathers potential colliders
      if (nrow(potential_colliders)>0 && sum(potential_colliders$exp_out_pair== paste0(exposure, "-", outcome))>0){
      colliders=c(as.character(unlist(potential_colliders$colliders[potential_colliders$exp_out_pair== paste0(exposure, "-", outcome)])))
      colliders_lagged=potential_colliders$exclude_lags[potential_colliders$exp_out_pair== paste0(exposure, "-", outcome)]

      # if (sum(grepl("\\.", colliders))>0){
      time_spec=colliders[grepl("\\.", colliders)]
      all_time_colliders=colliders[colliders %in% time_varying_covariates]
      time_invariant_colliders=colliders[!colliders %in% all_time_colliders & !colliders %in% time_spec]

      }else{
        colliders=NULL
        colliders_lagged=F
        time_varying_colliders=NULL
      }



      #Cycles through all exposure time points
      for (x in 1:length(exposure_time_pts)){
        time=exposure_time_pts[x]

        #finds concurrent time-varying and time-invariant potential confounds
        concurrent_covariates=na.omit(exposure_covariates[exposure_covariates$exp_time==time,])
        concurrent_covariates=c(concurrent_covariates$row, concurrent_covariates$column)
        concurrent_covariates=unique(concurrent_covariates)

        #splits covariates by time
        # time_varying_concurrent_covariates=paste0(concurrent_covariates[grepl(time_varying_covariates,concurrent_covariates)], paste0(".", time))
        time_varying_concurrent_covariates=concurrent_covariates[grepl(paste0(".",time), concurrent_covariates)]
        time_invariant_concurrent_covariates=concurrent_covariates[!concurrent_covariates %in% time_varying_concurrent_covariates]

        # #any mandatory covariates that are time-invariant
        # keep_time_invar_covars=mandatory_keep_covariates[!mandatory_keep_covariates %in% time_varying_covariates]


        #Finds relevant lagged time points
        lags=c(exposure_time_pts[exposure_time_pts<time])
        lagged_vars={} #variable for lagged variables
        if (length(lags)==0){ #no lags
          vars_to_include=c(time_varying_concurrent_covariates, time_invariant_concurrent_covariates) #retains only concurrent variables

        }else{ #lags present
          past_exposures=apply(expand.grid(exposure, as.character(lags)), 1, paste, sep="", collapse=".") #finds past exposure variables
          # keep_time_var_covars=apply(expand.grid(mandatory_keep_covariates[mandatory_keep_covariates %in% time_varying_covariates], as.character(lags)), 1, paste, sep="", collapse=".") #finds past time-varying mandatory covariates


          #populates data frames
          lagged_covariates=past_exposures
          vars_to_include=c(time_varying_concurrent_covariates, time_invariant_concurrent_covariates, lagged_covariates)
        }


        #EXCLUSION STAGE


        #gathering lagged colliders for exclusion if applicable
        if (colliders_lagged=="T" & length(all_time_colliders)>0){
          lagged_colliders=apply(expand.grid(all_time_colliders, as.character(lags)), 1, paste, sep="", collapse=".")
        }else{
          lagged_colliders=NULL
        }

        concurrent_colliders=c(time_spec[as.numeric(sapply(strsplit(time_spec, "\\."), "[", 2))==time], time_invariant_colliders)

        #USER: Place where you can manually delete variables that should not be included in forms such as index variables and possible colliders (e.g., outcomes that could also cause a given tx )
        vars_to_include=vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                                time_varying_covariates, #exclude time-varying covariates in long form (already in wide form)
                                                                paste(exposure, time, sep="."), #exclude exposure at current time point
                                                                apply(expand.grid(outcome, as.character(time_pts)), 1, paste, sep="", collapse="."), #exclude outcome at any time point
                                                                exclude_covariates[grepl(exposure, exclude_covariates)==F], #exclude any covariates specified by the user
                                                                time_var_exclude, #exclude any time points that should not be there bc of planned missingness
                                                                concurrent_colliders, #exclude colliders specified at exposure time pt
                                                                # apply(expand.grid(colliders, as.character(outcome_time_pt)), 1, paste, sep="", collapse="."), #exclude colliders at outcome time points
                                                                lagged_colliders #any lagged colliders as indicated
                                                                )]

        #Creates form for the given exposure time point
        f=paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + "))

        #prints form for user inspection
        cat(paste0("Formula for exposure ", exposure, "-", outcome, " at time point ", as.character(time),
                     ":"), "\n")
        print(f)
        cat("\n")

        #saves out form to global workspace and local directory
        # return(assign(paste("form_", exposure, "_", time, sep=""), f))
        forms[[paste("form_", exposure,"-", outcome, "_", time, sep="")]] <- f

        write.csv(f, paste0(home_dir, "forms/form_", exposure, "-", outcome, "_", time, ".csv", sep=""))

      }
    }
  }

  return(forms)
}


