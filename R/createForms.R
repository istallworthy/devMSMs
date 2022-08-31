#code to create CBPS forms at each time point for each exposure
#' createForms
#'
#' @return forms
#' @export
#' @param wide_long_datasets from formatForWeights
#' @param covariates_to_include from identifyPotentialConfounds
#' @param exposures list of exposures of interest
#' @param outcomes list of outcomes of interest
#' @param time_pts list of time points
#' @param potential_colliders optional list of variables to be excluded from balancing at time point of exposure

#' @examples
createForms <- function(wide_long_datasets, covariates_to_include, exposures, outcomes, time_pts, potential_colliders=NULL){

  ####Creates CBPS forms for each treatment at each time point (e.g., HOMEETA.6) that includes: all identified potential confounds for that treatment (at any time point) and lagged time points of the treatment, excludes other outcomes a the given time point (as potential colliders).
  ####Call tx variables and lagged variables for each tx. Creates variables form_tx_time for each tx and time point. Could be helpful to spot check a few of these forms in the global environment.

  #determining column names
  imp_col_names=colnames(wide_long_datasets[[1]])

  forms=list()

  #Cycles through all exposures
  for (y in 1:length(exposures)){
    exposure=exposures[y]

    #Cycles through all time points
    for (x in 1:length(time_pts)){
      time=time_pts[x]

      #Finds relevant concurrent and lagged time points
      lags=as.data.frame(time_pts[time_pts<time | time_pts==time])

      lagged_vars={}

      covariates=apply(expand.grid(time_varying_covariates, as.character(lags)), 1, paste, sep="", collapse=".")
      #exclude exposure at that time point
      covariates=covariates[!covariates %in% paste0(exposure, ".", time)]

      #Checks to see if lagged time-varying covariates exist in widelong dataset
      if (sum(covariates %in% imp_col_names)<length(covariates)){
        message(paste0("The following lagged covariates for ", exposure, " at time ", time, " were not found in the dataset and thus skipped: ", c(covariates[!covariates %in% imp_col_names])))
      }

      lagged_vars=rbind(lagged_vars, covariates[covariates %in% imp_col_names])

      #Makes list of variables to include in form for a given tx and time point: USER CHECK THIS!
      # vars_to_include=c(as.character(unlist(get(paste(exposure, "_covars_to_include", sep="")))), lagged_vars)
      vars_to_include=c(covariates_to_include[[paste(exposure, "_covars_to_include", sep="")]], lagged_vars)

      #Removes time-varying vars that have now been coverted to wide;
      #USER: Place where you can manually delete variables that should not be included in forms such as index variables and possible colliders (e.g., outcomes that could also cause a given tx )
      vars_to_include=vars_to_include[!vars_to_include %in% c(ID, "WAVE",
                                                              time_varying_covariates, #exclude time-varying covariates in long form (already in wide)
                                                              paste(exposure, time_pt, sep=""), #exclude exposure at that time point
                                                              paste(outcomes, time, sep="."), #exclude any outcomes at that time point --collider bias
                                                              apply(expand.grid(potential_colliders, as.character(time)), 1, paste, sep="", collapse="."))] #exclude any variables deemed colliders


      #Creates form for the given exposure time point
      f=paste(exposure, "~", paste0(vars_to_include, sep="", collapse=" + "))

      #prints form for user inspection
      message(paste0("Please inspect the weights formula below for exposure ", exposure, " at time point ", as.character(time),
                   ". Balancing weights will attempt to balance on all of these potential confounding variables. If there are any time-varying variables you wish to omit at this time point, please list it in the 'potential_colliders' field and re-run"))
      print(f)

      #saves out form to global workspace and local directory
      # return(assign(paste("form_", exposure, "_", time, sep=""), f))
      forms[[paste("form_", exposure, "_", time, sep="")]] <- f

      write.csv(f, paste0(home_dir, "forms/form_", exposure, "_", time, ".csv", sep=""))
      message("Please see the 'forms' folder for csv files of each of the formulae to be used in creating the balancing weights")

    }
  }

  return(forms)
}


