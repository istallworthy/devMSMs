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

  options(readr.num_columns = 0)
  cat("\n")
  cat("USER ALERT: Inspect the list of time-varying covariates following the ID variable and remove any that should not be there because of planned missingness design by adding them to 'time_var_exclude' in the msmObject and re-running","\n")
  cat("\n")


  #makes hybrid "wide/long" dataset for each impute dataset: all time-varying covariates listed as long and wide

  #Cyles through imputed datasets and puts them in hybrid wide/long dataset
  wide_long_datasets=list()
  for (k in 1:m){

    imp=imputed_datasets[[paste0("imp", k)]]
    imp=as.data.frame(imp)

    time_varying_wide=apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep="", collapse=".")
    time_varying_wide=sort(time_varying_wide)
    time_varying_wide=c(ID, time_varying_wide)
    time_varying_wide=time_varying_wide[!time_varying_wide %in% time_var_exclude]

    #Make wide so that all time-varying now in wide, #Then select all the new wide time varying vars --IS SHOULD AUTOMATE THIS STEP ....
    #MAKE SURE ALL YOUR TX AND OUTCOME VARIABLES AT EACH TIME PT ARE LISTED BELOW IN WIDE FORMAT
    # browser()
    imp_wide=suppressWarnings(stats::reshape(data=imp,
                     idvar=ID,
                     v.names= time_varying_covariates, #list ALL time-varying covariates
                     timevar="WAVE",
                     times=c(time_pts),
                     direction="wide"))%>%
                    dplyr::select(dput(as.character(time_varying_wide)))

      # dplyr::select(cat(paste0("c(", unlist(paste0('"',time_varying_wide, '"', collapse=", ")), ")"))))




    imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude]

    #creates hybrid wide/long dataset
    msm_data=merge(imp, imp_wide, by=ID, all.x=T) #dont delete rows

    #remove time points that should not be there (but may have been added by imputation)
    msm_data=msm_data[msm_data$WAVE %in% time_pts,]

    #re-labels time points
    for (x in 1:length(time_pts)){
      msm_data$WAVE[msm_data$WAVE==time_pts[x]]=x}

    msm_data=as.data.frame(msm_data)

    # return(assign(paste("imp", k, "_widelong", sep=""), msm_data)) #save each out e.g., imp1_widelong
    wide_long_datasets[[paste("imp", k, "_widelong", sep="")]] <- msm_data

  }



  #create dataset for future modeling
  #run this code to make wide dataset for future use in actual msm model (i.e., non-imputed wide data)
  # imp=as.data.frame(data)
  imp=as.data.frame(imputed_datasets$imp1) #IS changed to ue imputed dataset

  imp_wide=suppressWarnings(stats::reshape(data=imp,
                                           idvar=ID,
                                           v.names= time_varying_covariates, #list ALL time-varying covariates
                                           timevar="WAVE",
                                           times=c(time_pts),
                                           direction="wide"))
  imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude]


  #
  # imp_wide=plyr::join(test,t, by=ID)
  #Create long/wide hybrid: merge this newly created wide dataset with long dataset
  write.csv(imp_wide, paste0(home_dir, "data_for_final_model.csv"))
  write.csv(wide_long_datasets[[1]], paste0(home_dir, "data_for_final_Mplus_model.csv"))

  cat("\n")
  cat("A dataset has been saved out as a csv file in home directory for later modeling ")


  return(wide_long_datasets)

}
