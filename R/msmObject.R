
#' Create msm object
#' Gathers all information needed for the functions
#' @param data_path path to cleaned dataset
#' @param home_dir path to home directory for the project
#' @param ID person-level identifier in your dataset
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @param time_var variable in your long dataset that designates developmental time
#' @param missing missing data marker in your dataset
#' @param m number of imputations
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param exposure_time_pts times when exposures occurred
#' @param exposure_epochs data frame containing epochs and correponding time points for exposure histories
#' @param outcomes list of variables that represent your outcomes of interest
#' @param outcome_time_pts list of time points for outcomes
#' @param continuous_variables list of continuous variables in your dataset
#' @param time_var_exclude time-varying variables that should not be present because of planned missingness
#' @examples msmObject(data_path, home_dir, ID, time_pts, time_var, missing, m, exposures, exposure_time_pts, exposure_epochs, outcomes, outcome_time_pts, continuous_variables, time_var_exclude)

msmObject <- function(data_path, home_dir, ID, time_pts, time_var, missing, m, exposures, exposure_time_pts, exposure_epochs, outcomes, outcome_time_pts, continuous_variables=NULL, time_var_exclude=NULL){


  charOrNull <- function(x) {
    is.character(x) || is.null(x)
  }

  numOrNull <- function(x) {
    is.numeric(x) || is.null(x)
  }


  stopifnot(charOrNull(data_path))
  stopifnot(charOrNull(data_path))
  stopifnot(charOrNull(time_var))
  stopifnot(charOrNull(exposures))
  stopifnot(charOrNull(outcomes))
  stopifnot(charOrNull(continuous_variables))
  stopifnot(charOrNull(time_var_exclude))


  stopifnot(numOrNull(time_pts))
  stopifnot(numOrNull(m))
  stopifnot(numOrNull(exposure_time_pts))
  stopifnot(numOrNull(outcome_time_pts))

  object<-list(data_path=data_path, home_dir=home_dir, ID=ID, time_pts=time_pts, time_var=time_var, missing=missing, m=m, exposures=exposures,
               exposure_epochs=exposure_epochs, exposure_time_pts=exposure_time_pts, outcomes=outcomes, ID=ID, outcome_time_pts=outcome_time_pts,
               continuous_variables=continuous_variables, time_var_exclude=time_var_exclude)

  class(object)<-c("list","msmobject")

  return(object)

}
