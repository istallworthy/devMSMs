
#' Create msm object
#' Gathers all information needed for the functions and checks user inputs
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
#' @param outcome list of variables that represent your outcome of interest
#' @param outcome_time_pts list of time points for outcome
#' @param continuous_variables list of continuous variables in your dataset
#' @param time_var_exclude time-varying variables that should not be present because of planned missingness
#' @param colors colors for plotting dose
#' @examples msmObject(data_file, home_dir, ID, time_pts, time_var, missing, time_varying_variables=NULL, continuous_variables=NULL, factor_covariates=NULL, m=5, imp_method="cart", exposure, exposure_time_pts, short_form_lag=1, weights_method="cbps", balance_thresh=0.1, weights_percentile_cutoff=0.95, exposure_epochs, reference=NA, comparisons=NA, hi_cutoff=.75,lo_cutoff=.25, mc_method="BH", outcome, outcome_time_pt, keep_concurrent_tv_vars=NULL, time_var_exclude=NULL,  exposure_labels=NULL, outcome_labels=NULL, dose_level="h", colors="Dark2")

msmObject <- function(data_file, home_dir, ID, time_pts, time_var, missing, time_varying_variables=NULL, continuous_variables=NULL, factor_covariates=NULL, m=5, imp_method="cart", exposure, exposure_time_pts, short_form_lag=1, weights_method="cbps", balance_thresh=0.1, weights_percentile_cutoff=0.95, exposure_epochs, reference=NA, comparisons=NA, hi_cutoff,lo_cutoff, mc_method="BH", outcome, outcome_time_pt, keep_concurrent_tv_vars=NULL, time_var_exclude=NULL,  exposure_labels=NULL, outcome_labels=NULL, dose_level="h", colors="Dark2"){


  if (!file.exists(home_dir)){
    stop('Please enter a valid home directory.')}

   data_path=paste0(home_dir, data_file)

  if (!file.exists(data_path)){
    stop('Please enter a valid data file that is in the home directory.')}


  charOrNull <- function(x) {
    is.character(x) || is.null(x)
  }
  numOrNull <- function(x) {
    is.numeric(x) || is.null(x)
  }
  logOrNull <- function(x) {
    is.logical(x) || is.null(x)
  }
  #required
  stopifnot(is.character(data_path))
  stopifnot(is.character(home_dir))
  stopifnot(is.character(ID))
  stopifnot(is.character(time_var))
  stopifnot(is.character(exposure))
  stopifnot(is.character(outcome))

  stopifnot(is.numeric(time_pts))
  stopifnot(is.numeric(exposure_time_pts))
  stopifnot(is.numeric(outcome_time_pt))
  stopifnot(is.numeric(short_form_lag))
  stopifnot(is.numeric(hi_cutoff))
  stopifnot(is.numeric(lo_cutoff))


#optional
  stopifnot(charOrNull(time_varying_variables))
  stopifnot(charOrNull(continuous_variables))
  stopifnot(charOrNull(factor_covariates))
  stopifnot(charOrNull(time_var_exclude))
  stopifnot(charOrNull(mc_method))
  stopifnot(charOrNull(colors))
  stopifnot(charOrNull(dose_level))
  stopifnot(charOrNull(keep_concurrent_tv_vars))
  stopifnot(charOrNull(imp_method))
  stopifnot(charOrNull(weights_method))


  stopifnot(numOrNull(m))
  stopifnot(numOrNull(balance_thresh))
  stopifnot(numOrNull(weights_percentile_cutoff))

  #creates 2 additional weights truncation values other than what the user supplied for sensitivity analyses
  weights_percentile_cutoffs_sensitivity=ifelse(weights_percentile_cutoff+0.03<1,
                                                     paste(weights_percentile_cutoff-0.03, weights_percentile_cutoff+0.03, sep=" ", collapse=" "),
                                                     paste(weights_percentile_cutoff-0.06, weights_percentile_cutoff-0.03, sep=" ", collapse=" "))


  object<-list(data_path=data_path, home_dir=home_dir, ID=ID, time_pts=time_pts, time_var=time_var,
               missing=missing, time_varying_variables=time_varying_variables,
               continuous_variables=continuous_variables,factor_covariates=factor_covariates, m=m, imp_method=imp_method,
               exposure=exposure,exposure_time_pts=exposure_time_pts, short_form_lag=short_form_lag, weights_method=weights_method,
               balance_thresh=balance_thresh, weights_percentile_cutoff=weights_percentile_cutoff,
               weights_percentile_cutoffs_sensitivity=weights_percentile_cutoffs_sensitivity,
               exposure_epochs=exposure_epochs, reference=reference, comparisons=comparisons, hi_cutoff=hi_cutoff,
               lo_cutoff=lo_cutoff, mc_method=mc_method,
               outcome=outcome, outcome_time_pt=outcome_time_pt,
               keep_concurrent_tv_vars=keep_concurrent_tv_vars, time_var_exclude=time_var_exclude,
               xposure_labels=exposure_labels, outcome_labels=outcome_labels,
               dose_level=dose_level, colors=colors)

  class(object)<-c("list","msmobject")

  return(object)

}
