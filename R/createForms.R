#' Creates formulae for generating weights
#'
#' Creates formulae for each exposure at each time point for the CBPS function that will make weights returning a list of forms
#'
#' @param wide_long_datasets from formatForWeights
#' @param covariates_to_include from identifyPotentialConfounds
#' @param exposures list of exposures of interest
#' @param outcomes list of outcomes of interest
#' @param time_pts list of time points
#' @param potential_colliders optional list of variables to be excluded from balancing at time point of exposure
#' @return forms
#' @export
#' @seealso [formatForWeights()], [identifyPotentialConfounds] for more on inputs
#' @examples createForms(wide_long_datasets,covariates_to_include, exposures, outcomes, time_pts, potential_colliders)
#'
createForms <- function(wide_long_datasets, covariates_to_include, exposures, outcomes, time_pts, potential_colliders=NULL, keep_covariates=NULL, time_var_exclude){

  ####Creates CBPS forms for each exposure at each time point (e.g., HOMEETA.6) that includes: all identified potential confounds for that treatment (at any time point) and lagged time points of the treatment, excludes other outcomes a the given time point (as potential colliders).

  #determining column names
  imp_col_names=colnames(wide_long_datasets[[1]])

  forms=list()

  print("USER ALERT: Please inspect the weights formula below. Balancing weights will attempt to balance on all of these potential confounding variables. If there are any time-varying variables you wish to omit at this time point, please list it in the 'potential_colliders' field and re-run")
  print("They are also saved out in 'forms' folder as csv files")

  #cycles through outcomes
  for (z in seq(length(outcomes))){
    outcome=outcomes[z]

    #Cycles through all exposures
    for (y in 1:length(exposures)){
      exposure=exposures[y]

      exposure_covariates=covariates_to_include[[paste(exposures[y],"-", outcomes[z], "_covars_to_include", sep="")]]


      #Cycles through all time points
      for (x in 1:length(time_pts)){
        time=time_pts[x]

        #finds concurrent time-varying and time-invariant potential confounds
        concurrent_covariates=exposure_covariates[exposure_covariates$time==time,]
        concurrent_covariates=c(concurrent_covariates$row, concurrent_covariates$column)
        concurrent_covariates=unique(concurrent_covariates)
        time_varying_concurrent_covariates=paste0(concurrent_covariates[concurrent_covariates %in% time_varying_covariates], paste0(".", time))
        time_invariant_concurrent_covariates=concurrent_covariates[!concurrent_covariates %in% time_varying_covariates]

        keep_time_invar_covars=keep_covariates[!keep_covariates %in% time_varying_covariates]

        #Finds relevant lagged time points
        lags=c(time_pts[time_pts<time])

        lagged_vars={}
        if (length(lags)==0){
          vars_to_include=c(time_varying_concurrent_covariates, time_invariant_concurrent_covariates, keep_time_invar_covars)

        }else{
          past_exposures=apply(expand.grid(exposure, as.character(lags)), 1, paste, sep="", collapse=".") #finds past exposures
          keep_time_var_covars=apply(expand.grid(keep_covariates[keep_covariates %in% time_varying_covariates], as.character(lags)), 1, paste, sep="", collapse=".")

          #finds past exposures
          lagged_covariates=past_exposures
          vars_to_include=c(time_varying_concurrent_covariates, time_invariant_concurrent_covariates, lagged_covariates, keep_time_var_covars)
        }

        # covariates=apply(expand.grid(time_varying_covariates, as.character(lags)), 1, paste, sep="", collapse=".")
        #exclude exposure at that time point
        # covariates=covariates[!covariates %in% paste0(exposure, ".", time)]

        # #Checks to see if lagged time-varying covariates exist in widelong dataset
        # lagged_vars=rbind(lagged_vars, covariates[covariates %in% imp_col_names])



        #Makes list of variables to include in form for a given tx and time point: USER CHECK THIS!
        # vars_to_include=c(covariates_to_include[[paste(exposure, "_covars_to_include", sep="")]], lagged_vars)


        #Removes time-varying vars that have now been coverted to wide;
        #USER: Place where you can manually delete variables that should not be included in forms such as index variables and possible colliders (e.g., outcomes that could also cause a given tx )
        vars_to_include=vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                                time_varying_covariates, #exclude time-varying covariates in long form (already in wide)
                                                                paste(exposure, time, sep="."), #exclude exposure at that time point
                                                                paste(outcome, time, sep="."), #exclude any outcomes at that time point
                                                                time_var_exclude, #exclude any time points that should not be there bc of planned missingness
                                                                apply(expand.grid(potential_colliders, as.character(time)), 1, paste, sep="", collapse="."))] #exclude any variables deemed colliders

        #Creates form for the given exposure time point
        f=paste(exposure, "~", paste0(vars_to_include[order(vars_to_include)], sep="", collapse=" + "))

        #prints form for user inspection
        print(paste0("Formula for exposure ", exposure, "-", outcome, " at time point ", as.character(time),
                     ": "))
        print(f)

        #saves out form to global workspace and local directory
        # return(assign(paste("form_", exposure, "_", time, sep=""), f))
        forms[[paste("form_", exposure,"-", outcome, "_", time, sep="")]] <- f

        write.csv(f, paste0(home_dir, "forms/form_", exposure, "-", outcome, "_", time, ".csv", sep=""))

      }
    }
  }

  return(forms)
}


