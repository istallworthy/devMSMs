
#' Create msm object
#' Gathers all information needed for the functions
#' @param data_path path to cleaned dataset
#' @param home_dir path to home directory for the project
#' @param ID person-level identifier in your dataset
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @param time_var variable in your long dataset that designates developmental time
#' @param missing missing data marker in your dataset
#' @param time_varying_variables all variables in dataset that are time-varying
#' @param continuous_variables all continous variables
#' @param factor_covariates all variables that are factors
#' @param m number of imputations
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param exposure_time_pts times when exposures occurred
#' @param balance_thresh correlation value above which covariates are not considered balanced with respect to exposure
#' @param weights_percentile_cutoff percentile value for cutting off and replacing heavy tails of weights
#' @param exposure_epochs data frame containing epochs and correponding time points for exposure histories
#' @param outcomes list of variables that represent your outcomes of interest
#' @param outcome_time_pts list of time points for outcomes
#' @param continuous_variables list of continuous variables in your dataset
#' @param time_var_exclude time-varying variables that should not be present because of planned missingness
#' @param colors colors for plotting dose
#' @examples msmObject(data_path, home_dir, ID, time_pts, time_var, missing, m, exposures, exposure_time_pts, exposure_epochs, outcomes, outcome_time_pt, continuous_variables, time_var_exclude)

msmObject <- function(data_path, home_dir, ID, time_pts, time_var, missing, time_varying_variables=NULL, continuous_variables=NULL, factor_covariates=NULL, m=5, exposures, exposure_time_pts,  balance_thresh=0.12, weights_percentile_cutoff=0.90, exposure_epochs, reference="", comparisons="", hi_cutoff=.75,lo_cutoff=.25, mc_method="BH", outcomes, outcome_time_pt, mandatory_keep_covariates=NULL, exclude_covariates=NULL, time_var_exclude=NULL, potential_colliders=NULL, exposure_labels=NULL, outcome_labels=NULL, colors="Dark2"){

  charOrNull <- function(x) {
    is.character(x) || is.null(x)
  }
  numOrNull <- function(x) {
    is.numeric(x) || is.null(x)
  }

  #required
  stopifnot(is.character(data_path))
  stopifnot(is.character(home_dir))
  stopifnot(is.character(ID))
  stopifnot(is.character(time_var))
  stopifnot(is.character(exposures))
  stopifnot(is.character(outcomes))

  stopifnot(is.numeric(time_pts))
  stopifnot(is.numeric(exposure_time_pts))
  stopifnot(is.numeric(outcome_time_pt))

#optional
  stopifnot(charOrNull(time_varying_variables))
  stopifnot(charOrNull(continuous_variables))
  stopifnot(charOrNull(factor_covariates))
  stopifnot(charOrNull(time_var_exclude))
  stopifnot(charOrNull(exclude_covariates))
  stopifnot(charOrNull(mandatory_keep_covariates))
  stopifnot(charOrNull(exclude_covariates))
  stopifnot(charOrNull(mc_method))
  stopifnot(charOrNull(colors))


  stopifnot(numOrNull(m))
  stopifnot(numOrNull(balance_thresh))
  stopifnot(numOrNull(weights_percentile_cutoff))


  object<-list(data_path=data_path, home_dir=home_dir, ID=ID, time_pts=time_pts, time_var=time_var, missing=missing, time_varying_variables=time_varying_variables,
               continuous_variables=continuous_variables,factor_covariates=factor_covariates, m=m, exposures=exposures,exposure_time_pts=exposure_time_pts, balance_thresh=balance_thresh, weights_percentile_cutoff=weights_percentile_cutoff, exposure_epochs=exposure_epochs, reference=reference, comparisons=comparisons, hi_cutoff=hi_cutoff, lo_cutoff=lo_cutoff, mc_method=mc_method,
               outcomes=outcomes, outcome_time_pt=outcome_time_pt,mandatory_keep_covariates=mandatory_keep_covariates, exclude_covariates=exclude_covariates, time_var_exclude=time_var_exclude, potential_colliders=potential_colliders, exposure_labels=exposure_labels, outcome_labels=outcome_labels, colors=colors)

  class(object)<-c("list","msmobject")

  return(object)

}
