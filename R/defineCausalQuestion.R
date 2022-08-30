#This code helps the individual delineate core aspects of their causal question

#' Title
#'
#' @param exposure string of the core exposure/treatment of interest variable
#' @param exposure_time_pts sequence of integer time points of exposure
#' @param outcome string of core outcome of interest variable
#' @param outcome_time_pts sequence of integer time points of outcome
#'
#' @return
#' @export
#'
#' @examples
#'
defineCausalQuestion<- function(exposures, exposure_time_pts, outcomes, outcome_time_pts){

  for (i in 1:length(exposures)){
    for (x in 1:length(outcomes)){
      print(paste0("Goal: examine causal effects of ", exposures[i], " at times ", paste0(exposure_time_pts,sep=",", collapse=""), " on ", outcomes[x], " at time", outcome_time_pts))
    }}
  return(exposures)
  return(exposure_time_pts)
  return(outcomes)
  return(outcome)
}
