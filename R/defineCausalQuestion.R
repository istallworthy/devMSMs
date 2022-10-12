#' Define causal questions at the heart of the subsequent marginal structural models
#'
#' This code helps the individual delineate core aspects of their causal questions
#'
#' @param object msm object that contains all relevant user inputs
#' @export
#' @examples defineCausalQuestion(object)
#'
defineCausalQuestion<- function(object){

  exposures=object$exposures
  exposure_time_pts=object$exposure_time_pts
  outcomes=object$outcomes
  outcome_time_pts=object$outcome_time_pts

  for (i in 1:length(exposures)){
    for (x in 1:length(outcomes)){
      print(paste0("Goal: examine causal effects of ", exposures[i], " at times ", paste0(exposure_time_pts,sep=",", collapse=""), " on ", outcomes[x], " at time(s) ", paste(c(outcome_time_pts), sep=",", collapse=",")))
    }}

}
