#this code formats the imputed datasets for calculating balancing weights
#' Title
#'
#' @param ID identifier
#' @param home_dir path to home directory for project
#' @param imputed_datasets
#' @param time_varying_covariates
#' @param time_pts identifier
#' @param time_var_exclude
#' @param just_imputed "yes"= you have imputed datasets in global environment or "no" but they are saved locally from previous run

#' @return
#' @export
#' @importFrom readr read_csv
#' @examples
#'
formatForWeights <- function(ID, home_dir, imputed_datasets, time_varying_covariates, time_pts, time_var_exclude=NULL, just_imputed="yes"){

  #if the user has not just imputed datasets (and imputations are instead saved locally from a prior run), read in imputed data
  if (just_imputed=="no"){
    for (x in 1:5){
      file_name=(paste("imp", x, '.csv', sep=""))
      name=paste("imp", x, sep="")
      imp=readr::read_csv(paste(paste0(home_dir, "imputations/"), file_name, sep=""))
      assign(paste("Imp", x, sep=""),imp)
      imp=NULL

      return(Imp1)
      return(Imp2)
      return(Imp3)
      return(Imp4)
      return(Imp5)
    }
  }

  #makes hybrid "wide/long" dataset for each impute dataset: all time-varying covariates listed as long and wide

  #Cyles through imputed datasets and puts them in hybrid wide/long dataset
  #USER INPUT: make sure that all time-varying covariates are listed below (long) and in the select() list (wide).
  #Also, make sure that all variables that should not be there (e.g., outcomes in final MSM models) are made to NULL below.

  for (k in 1:5){

    imp=get(paste("Imp", k, sep="")) #retrieves imputed data
    imp=as.data.frame(imp)

    time_varying_wide=apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep="", collapse=".")
    time_varying_wide=sort(time_varying_wide)
    time_varying_wide=c(ID, time_varying_wide)

    #Make wide so that all time-varying now in wide, #Then select all the new wide time varying vars --IS SHOULD AUTOMATE THIS STEP ....
    #MAKE SURE ALL YOUR TX AND OUTCOME VARIABLES AT EACH TIME PT ARE LISTED BELOW IN WIDE FORMAT
    imp_wide=reshape(data=imp,
                     idvar=ID,
                     v.names= time_varying_covariates, #list ALL time-varying covariates
                     timevar="WAVE",
                     times=c(time_pts),
                     direction="wide")%>%
      dplyr::select(.dots=dput(as.character(time_varying_wide)))


    print("Inspect the list below of time-varying covariates and remove any that should not be there because of planned missingness design")
    print(time_varying_wide)

    imp_wide=imp_wide[,!colnames(imp_wide) %in% time_var_exclude]

    #creates hybrid wide/long dataset
    msm_data=merge(imp, imp_wide, by="s_id", all.x=T) #dont delete rows

    # #Delete the timevarying long covariates since we focus on lags in the forms below
    # # MSMDATA$AGE=NULL
    # # MSMDATA$INR=NULL
    # MSMDATA$WIND=NULL
    # MSMDATA$zid=NULL
    # MSMDATA$pcx_sensitive=NULL
    # MSMDATA$pcx_neg=NULL
    # MSMDATA$pcx_pos=NULL
    # MSMDATA$pcx_pomd=NULL
    # MSMDATA$pcx_nomd=NULL
    # MSMDATA$pcx_engaged=NULL
    # MSMDATA$pcx_qrel=NULL
    # MSMDATA$pcx_CompTwo=NULL
    # MSMDATA$GrosPay1=NULL
    # MSMDATA$sAABASE=NULL

    #remove time points that should not be there (but may have been added by imputation)
    msm_data=msm_data[msm_data$WAVE %in% time_pts,]

    #re-labels time points
    for (x in 1:length(time_pts)){
    msm_data$WAVE[msm_data$WAVE==time_pts[x]]=x}

    msm_data=as.data.frame(msm_data)

    return(assign(paste("imp", k, "_widelong", sep=""), msm_data)) #save each out e.g., imp1_widelong

  }

}
