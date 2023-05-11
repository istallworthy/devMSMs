#' Format data for creating balancing weights
#'
#' This code formats the imputed datasets for calculating balancing weights using CBPS and reutrns a list of wide/long datasets and creates dataset for future modeling
#'
#' @param object msm object that contains all relevant user inputs
#' @param imputed_datasets output from imputeData
#' @param just_imputed "yes"= you have imputed datasets in global environment or "no" but they are saved locally from previous run
#' @return wide_long_datasets
#' @export
#' @importFrom readr read_csv
#' @importFrom stats reshape
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tidyr pivot_wider
#' @importFrom plyr join
#' @seealso [formatDataStruct()], [imputeData()]
#' @examples formatForWeights(object, data, imputed_datasets)
#'
formatForWeights <- function(object, data, imputed_datasets){

  ID=object$ID
  home_dir=object$home_dir
  m=object$m
  time_pts=object$time_pts
  time_var_exclude=object$time_var_exclude
  time_varying_covariates=object$time_varying_variables
  exposure=object$exposure
  outcome=object$outcome

  options(readr.num_columns = 0)
  cat("\n")
  cat("USER ALERT: Inspect the list of time-varying covariates (following the ID variable) and remove any that should not be there because of planned missingness design by adding them to 'time_var_exclude' in the msmObject and re-running","\n")
  cat("\n")


  #makes wide for each impute dataset: all time-varying covariates listed as long and wide

  #Cyles through imputed datasets and puts them in hybrid wide/long dataset
  wide_long_datasets<-lapply(1:m, function(k){
    imp=mice::complete(imputed_datasets,k)
    imp=as.data.frame(imp)

    time_varying_wide=apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep="", collapse=".")
    time_varying_wide=sort(time_varying_wide)
    time_varying_wide=c(ID, time_varying_wide)
    time_varying_wide=time_varying_wide[!time_varying_wide %in% time_var_exclude] #removes time-varying time pts that should not be there

    #Make wide so that all time-varying now in wide, #Then select all the new wide time varying vars --IS SHOULD AUTOMATE THIS STEP ....
    #MAKE SURE ALL YOUR TX AND OUTCOME VARIABLES AT EACH TIME PT ARE LISTED BELOW IN WIDE FORMAT
    # browser()
    imp_wide=suppressWarnings(stats::reshape(data=imp,
                                             idvar=ID,
                                             v.names= time_varying_covariates, #list ALL time-varying covariates
                                             timevar="WAVE",
                                             times=c(time_pts),
                                             direction="wide"))

    imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude] #only include what should be there
    msm_data=imp_wide
    msm_data=as.data.frame(msm_data)
  })
names(wide_long_datasets)<-1:m

  #create dataset for future modeling
  imp=as.data.frame(mice::complete(imputed_datasets,1)) #IS changed to use imputed dataset

  print(colnames(imp))

  imp_wide=suppressWarnings(stats::reshape(data=imp,
                                           idvar=ID,
                                           v.names= time_varying_covariates, #list ALL time-varying covariates
                                           timevar="WAVE",
                                           times=c(time_pts),
                                           direction="wide"))
  imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude]


  #
  write.csv(imp_wide, paste0(home_dir,  exposure, "-", outcome, "_data_for_final_model.csv"))
  # write.csv(wide_long_datasets[[1]], paste0(home_dir, "data_for_final_Mplus_model.csv"))

  cat("\n")
  cat("A dataset has been saved out as a csv file in home directory for later modeling ")


  return(wide_long_datasets)

}
