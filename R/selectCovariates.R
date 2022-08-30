#The functions below allow the user to select the appropriate covariates, or potential confounds in their
#causal investigation


#' Title
#'
#' @param data
#' @param potential_covariates
#'
#' @return
#' @export
#'
#' @examples
identifyCovariates <- function(home_dir, data, time_pts, potential_covariates, factor_covariates, time_varying_covariates, exclude_covariates=NULL){

  #defines functions
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

  #creates subfolders in homedir for each imputed dataste
  # subfolder_names <- c("Imp1", "Imp2", "Imp3", "Imp4", "Imp5")
  # for (name in subfolder_names){
  #   path <- paste0(home_dir, name)
  #   if (dir.exists(path)==F){
  #     dir.create(path)
  #   }
  # }


  #Loop through different time point datasets and determine all correlations at each time point
  # and populate them all into covariate_correlations

  #generates all correlations above 0.1 for all data variables in dataset
  covariate_correlations={}
  for (x in 1:length(time_pts)){
    time_pt=time_pts[x]
    d=get(paste("Data_", time_pt, sep=""))
    c <- Hmisc::rcorr(as.matrix(d)) #makes corr table of all vars; cannot have dates or non-factored characters
    c=flattenCorrMatrix(c$r, c$P)
    filtered=c%>%
      filter(abs(c$cor)>0.1) #accounts for d=0.2 see Stuart paper for rationale
    filtered$time=time_pt
    covariate_correlations=rbind(covariate_correlations, filtered)
    filtered=NULL
    c=NULL
  }

  #for each exposure, list correlated covariates that include only that exposure and/or outcomes of interest
  variables_to_include=c(ID, "WAVE", exposures, outcomes)
  for (y in 1:length(exposures)){
    variables=covariate_correlations[covariate_correlations$row %in% c(exposures[y], outcomes) | covariate_correlations$column %in% c(exposures[y], outcomes),]
    unique_variables=unique(c(variables$row, variables$column))
    unique_variables=unique_variables[!unique_variables %in% c(exclude_covariates, "WAVE", ID)]
    variables_to_include=c(variables_to_include, unique_variables)
    assign(paste(exposures[y], "_covars_to_include", sep=""), unique_variables) #makes data frames for each tx with only relevant covariates
  }

  #creates final dataset with only relevant variables
  variables_to_include=unique(variables_to_include)
  data2=as.data.frame(data[names(data)[names(data) %in% variables_to_include] ])

  #inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  #View corr in final data for EACH TX to see if anything is highly correlated
  #creates correlation matrix of final dataset for inspection of highly correlated variables
  hi_corr_covars <- Hmisc::rcorr(as.matrix(data2[,3:ncol(data2)]))
  hi_corr_covars=flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)
  View_hi_corr_covars=hi_corr_covars%>%
    filter(abs(hi_corr_covars$cor)>0.6)

  print("Consider removing redundancies in the following covariates by listing in 'exclude_covariates' field above and re-running this function:")
  print(View_hi_corr_covars)

  data_to_impute=data2

  return(data_to_impute)
}
