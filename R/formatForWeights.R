#this code formats the imputed datasets for calculating balancing weights
#' formatForWeights
#'
#' @param ID identifier
#' @param home_dir path to home directory for project
#' @param m number of imputed datasets from Amelia
#' @param data output from formatDataStruct
#' @param imputed_datasets
#' @param time_varying_covariates
#' @param time_pts identifier
#' @param time_var_exclude
#' @param just_imputed "yes"= you have imputed datasets in global environment or "no" but they are saved locally from previous run

#' @return wide_long_datasets
#' @export
#' @importFrom readr read_csv
#' @importFrom stats reshape
#' @importFrom dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom plyr join
#' @examples
#'
formatForWeights <- function(ID, home_dir, m, data, imputed_datasets=list(), time_varying_covariates, time_pts, time_var_exclude=NULL, just_imputed="yes"){
  options(readr.num_columns = 0)

  #if the user has not just imputed datasets (and imputations are instead saved locally from a prior run), read in imputed data
  if (just_imputed=="no"){
    imputed_datasets=list()
    for (x in 1:m){
      file_name=(paste("imp", x, '.csv', sep=""))
      name=paste("imp", x, sep="")
      imp=as.data.frame(readr::read_csv(paste(paste0(home_dir, "imputations/"), file_name, sep="")))
      imputed_datasets[[paste0("imp", x)]]<-imp
    }
  }

  #makes hybrid "wide/long" dataset for each impute dataset: all time-varying covariates listed as long and wide

  #Cyles through imputed datasets and puts them in hybrid wide/long dataset
  #USER INPUT: make sure that all time-varying covariates are listed below (long) and in the select() list (wide).
  #Also, make sure that all variables that should not be there (e.g., outcomes in final MSM models) are made to NULL below.

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
    imp_wide=reshape(data=imp,
                     idvar=ID,
                     v.names= time_varying_covariates, #list ALL time-varying covariates
                     timevar="WAVE",
                     times=c(time_pts),
                     direction="wide")%>%
      dplyr::select(dput(as.character(time_varying_wide)))


    message("Inspect the list below of time-varying covariates and remove any that should not be there because of planned missingness design by adding them to 'time_var_exclude' and re-running")
    print(time_varying_wide)

    imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude]

    #creates hybrid wide/long dataset
    msm_data=merge(imp, imp_wide, by=ID, all.x=T) #dont delete rows

    # #Delete the time-varying long covariates since we focus on lags in the forms below
    # msm_data=msm_data[!colnames(msm_data) %in% time_varying_covariates]

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
  imp=as.data.frame(data)
  test=tidyr::pivot_wider(imp,
                   id_cols=ID,
                   names_from=WAVE,
                   values_from=c(time_varying_covariates))

  invar=colnames(imp)[colnames(imp) %in% time_varying_covariates==FALSE]
  invar=invar[!grepl("WAVE", invar)]
  invars=imp[invar] #pulls only invariant vars

  t=invars%>%
    dplyr::group_by(get(ID)) %>%
    dplyr::summarise(across(everything(), mean, na.rm=T))

  imp_wide=plyr::join(test,t, by=ID)
  #Create long/wide hybrid: merge this newly created wide dataset with long dataset
  # data_for_model=merge(imp, imp_wide, by="s_id", all.x=T) #dont delete rows
  write.csv(imp_wide, paste0(home_dir, "data_for_final_model.csv"))
  message("Saved out dataset for final modeling as a csv file in home directory")


  return(wide_long_datasets)

}
