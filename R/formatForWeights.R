#' Format data for creating balancing weights
#'
#' This code formats the imputed datasets for calculating balancing weights
#' @param object msm object that contains all relevant user inputs
#' @param imputed_datasets output from imputeData
#' @return wide_long_datasets
#' @export
#' @importFrom readr read_csv
#' @importFrom stats reshape
#' @importFrom mice complete
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

  #Cyles through imputed datasets and puts them in wide dataset
  wide_long_datasets<-lapply(1:m, function(k){
    imp=mice::complete(imputed_datasets,k)
    imp=as.data.frame(imp)

    #finds all time-varying variables in wide format
    time_varying_wide=apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep="", collapse=".")
    time_varying_wide=sort(time_varying_wide)
    time_varying_wide=c(ID, time_varying_wide)
    time_varying_wide=time_varying_wide[!time_varying_wide %in% time_var_exclude] #removes time-varying time pts that should not be there

    #creates wide dataset
    imp_wide=suppressWarnings(stats::reshape(data=imp,
                                             idvar=ID,
                                             v.names= time_varying_covariates, #list ALL time-varying covariates
                                             timevar="WAVE",
                                             times=c(time_pts),
                                             direction="wide"))
    imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude] #confirm: only include what should be there
    msm_data=imp_wide
    msm_data=as.data.frame(msm_data)
  })
  names(wide_long_datasets)<-1:m

  #create dataset for future modeling --not actually sure if we need this anymore tbh
  imp=as.data.frame(mice::complete(imputed_datasets,1))
  imp_wide=suppressWarnings(stats::reshape(data=imp,
                                           idvar=ID,
                                           v.names= time_varying_covariates, #list ALL time-varying covariates
                                           timevar="WAVE",
                                           times=c(time_pts),
                                           direction="wide"))
  imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude]
  # write.csv(imp_wide, paste0(home_dir,  exposure, "-", outcome, "_data_for_final_model.csv"))
  #
  # cat("\n")
  # cat("A dataset has been saved out as a csv file in home directory for later modeling ")


  return(wide_long_datasets)

}
