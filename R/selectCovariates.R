#' Identifies the available covariates
#'
#' Identifies the covariates in in the dataset that will be candidates for confounding
#'
#' identifyCovariates
#' @param data output from formatDataStruct
#' @param ID person-level identifier in your dataset
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param outcomes list of variables that represent your outcomes of interest
#' @param exclude_covariates list of variables to exclude based on theoretical or practical reasons
#' @return potential_covariates
#' @export
#' @seealso [formatDataStruct()]
#' @examples identifyCovariates(data, ID, exposures, outcomes, exclude_covariates)
#'
identifyCovariates <- function(data, ID, exposures, outcomes, exclude_covariates=NULL){
  #identifying potential covariates (time-varying and time invariant)
  potential_covariates=colnames(data)[colnames(data) %in% (c(ID, "WAVE", exposures, outcomes, exclude_covariates))==FALSE]

  print(paste0("Below are the ", as.character(length(potential_covariates)), " covariates that will be considered as potential confounding variables. If you would like to exclude any variables, please add them to 'exclude_covariates' and re-run. "))
  print(potential_covariates)

  return(potential_covariates)
}





#' Identify potential confounds
#'
#' Identifies potential confounds from the covariates based on those that are correlated with either exposures or outcomes above 0.1
#'
#' @param home_dir path to home directory for the project
#' @param data output from formatDataStruct
#' @param time_pt_datasets output from makeTimePtDatasets
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @param exclude_covariates list of variables to exclude based on theoretical or practical reasons
#' @return covariates_to_include
#' @export
#' @importFrom Hmisc rcorr
#' @seealso [formatDataStruct()], [makeTimePtDatasets()]
#' @examples identifyPotentialConfounds(home_dir, data, time_pt_datasets, time_pts, exclude_covariates=NULL)
#'
identifyPotentialConfounds <- function(home_dir, data, time_pt_datasets, time_pts, exclude_covariates=NULL){
  #Formats rcorr output for easy manipulation
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }

  #Loop through different time point datasets and determine all correlations at each time point
  #generates all correlations above 0.1 for all data variables in dataset
  covariate_correlations={}
  for (x in 1:length(time_pts)){
    time_pt=time_pts[x]
    # d=get(paste("Data_", time_pt, sep=""))
    d=time_pt_datasets[[x]] #gathers data at given time point
    d=d[,3:ncol(d)]
    d=as.data.frame(lapply(d, as.numeric))
    c <- suppressWarnings(Hmisc::rcorr(as.matrix(d))) #makes corr table of all vars; cannot have dates or non-factored characters
    c=flattenCorrMatrix(c$r, c$P)
    filtered=c%>%
      filter(abs(c$cor)>0.1) #accounts for d=0.2 see Stuart paper for rationale
    filtered$time=time_pt
    covariate_correlations=rbind(covariate_correlations, filtered)
    filtered=NULL
    c=NULL
  }

  #save out correlations
  write.csv(covariate_correlations, paste0(home_dir, "balance/covariate_correlations.csv"))
  print("Check the 'balance' folder to view a csv file of all correlations between covariates and exposures/outcomes > 0.1")

  #for each exposure, list correlated covariates that include only that exposure and/or outcomes of interest
  variables_to_include=c(ID, "WAVE", exposures, outcomes)
  covariates_to_include=list()
  for (y in 1:length(exposures)){
    variables=covariate_correlations[covariate_correlations$row %in% c(exposures[y], outcomes) | covariate_correlations$column %in% c(exposures[y], outcomes),]
    unique_variables=unique(c(variables$row, variables$column))
    unique_variables=unique_variables[!unique_variables %in% c(exclude_covariates, "WAVE", ID)]
    variables_to_include=c(variables_to_include, unique_variables)
    # assign(paste(exposures[y], "_covars_to_include", sep=""), unique_variables) #makes data frames for each tx with only relevant covariates
    covariates_to_include[[paste(exposures[y], "_covars_to_include", sep="")]] <- unique_variables #saves to list
  }

  return(covariates_to_include)
}



#' Creates data structure for imputation
#'
#' Creates data structure for imputation with only the variables needed (i.e., exposures, outcomes, potential confounds)
#'
#' @param ID person-level identifier in your dataset
#' @param data output from formatDataStruct
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param outcomes list of variables that represent your outcomes of interest
#' @param covariates_to_include output from identifyPotentialConfounds
#' @param exclude_covariates list of variables to exclude based on theoretical or practical reasons
#' @importFrom Hmisc rcorr
#' @importFrom knitr kable
#' @return data_to_impute
#' @export
#' @seealso [formatDataStruct()], [identifyPotenialConfounds()]
#' @examples dataToImpute(ID, data, exposures, outcomes, covariates_to_include, exclude_covariates)
#'
dataToImpute <-function(ID, data, exposures, outcomes, covariates_to_include, exclude_covariates){
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }

  #creates final dataset with only relevant variables
  variables_to_include=unique(c(ID, "WAVE", exposures, outcomes, unique(unlist(covariates_to_include))))
  data2=as.data.frame(data[names(data)[names(data) %in% variables_to_include] ])
  data2=data2[,!colnames(data2) %in% exclude_covariates]

  #inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  hi_corr_covars <- suppressWarnings(Hmisc::rcorr(as.matrix(data2[,3:ncol(data2)])))

  hi_corr_covars=flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)
  View_hi_corr_covars=hi_corr_covars%>%
    filter(abs(hi_corr_covars$cor)>0.6)

  print("USER ALERT: To simplify the balancing models consider removing any highly correlated, redundant covariates by listing them in the 'exclude_covariates' field above and re-running this function:")
  print(knitr::kable(View_hi_corr_covars))

  data_to_impute=data2
  write.csv(data_to_impute, paste0(home_dir, "imputations/data_to_impute.csv"))
  print("See the 'imputations' folder for a csv file of the data to be imputed")

  return(data_to_impute)
}
