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
#' @importFrom stargazer stargazer
#' @seealso [formatDataStruct()], [makeTimePtDatasets()]
#' @examples identifyPotentialConfounds(home_dir, data, time_pt_datasets, time_pts, exclude_covariates=NULL)
#'
identifyPotentialConfounds <- function(home_dir, data, exposures, outcomes, time_pt_datasets, time_pts, time_varying_covariates,exclude_covariates=NULL){

  user_input_exclude_covariates=exclude_covariates
  exclude_covariates=NULL
  if (length(user_input_exclude_covariates)>0){
    for (c in 1:length(user_input_exclude_covariates)){
      if (grepl("\\.",user_input_exclude_covariates[c])){ #period indicates user specified something time-varying
        exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_exclude_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified

          temp=apply(expand.grid(user_input_exclude_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          exclude_covariates=c(exclude_covariates,temp)
        }else{ #otherwise it is time invariant
          exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c])
        }
      }
    }
  }



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

  #keeps only those that involve either an exposure or outcome
  # test2=subset(covariate_correlations, subset=grepl(exposures, row) | grepl(exposures, column) | grepl(outcomes, row) | grepl(outcomes, column))
  covariate_correlations=covariate_correlations[covariate_correlations$row %in% c(exposures, outcomes) | covariate_correlations$column %in% c(exposures,outcomes),]
  covariate_correlations=covariate_correlations[!covariate_correlations$row %in% c(exclude_covariates) & !covariate_correlations$column %in% c(exclude_covariates),]

  #save out correlations
  stargazer::stargazer(covariate_correlations,type="html", digits=2, column.labels = colnames(covariate_correlations),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/all_potential_counfound_correlations.html"))
  # write.csv(covariate_correlations, paste0(home_dir, "balance/covariate_correlations.csv"))
  print("Check the 'balance' folder to view an html file of a table of all correlations between potential confounds and exposures/outcomes > 0.1 at each time point")

  #for each exposure-outcome pair, list correlated covariates that include only that exposure and/or outcomes of interest
  variables_to_include=c(ID, "WAVE", exposures, outcomes)
  covariates_to_include=list()

  for (z in seq(length(outcomes))){

    for (y in 1:length(exposures)){
      variables=covariate_correlations[covariate_correlations$row %in% c(exposures[y], outcomes[z]) | covariate_correlations$column %in% c(exposures[y], outcomes[z]),]
      unique_variables=unique(c(variables$row, variables$column)) #finds variables correlated at any time point
      unique_variables=unique_variables[!unique_variables %in% c(user_input_exclude_covariates, "WAVE", ID)]
      variables_to_include=c(variables_to_include, unique_variables)
      # assign(paste(exposures[y], "_covars_to_include", sep=""), unique_variables) #makes data frames for each tx with only relevant covariates
      # covariates_to_include[[paste(exposures[y], "_covars_to_include", sep="")]] <- unique_variables #saves to list
      covariates_to_include[[paste(exposures[y],"-", outcomes[z], "_covars_to_include", sep="")]] <- variables #saves relevant covariates at different time points

      #save out correlations by exposure-outcome pairing
      stargazer::stargazer(covariate_correlations,type="html", digits=2, column.labels = colnames(covariate_correlations),summary=FALSE, rownames = FALSE, header=FALSE,
                           out=paste0(home_dir, "balance/", exposures[y],"-", outcomes[z], "potential_counfound_correlations.html"))
      print(paste0("Check the 'balance' folder to view an html file of a table of all correlations between potential confounds for ", exposures[y],"-", outcomes[z], " at > 0.1 for each time point"))

    }
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
dataToImpute <-function(ID, data, exposures, outcomes, covariates_to_include, exclude_covariates=NULL, keep_covariates=NULL){
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }

  #assesses and reformat use input of exclude_covariates
  user_input_exclude_covariates=exclude_covariates
  exclude_covariates=NULL
  if (length(user_input_exclude_covariates)>0){
    for (c in 1:length(user_input_exclude_covariates)){
      if (grepl("\\.",user_input_exclude_covariates[c])){ #period indicates user specified something time-varying
        exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_exclude_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified

          temp=apply(expand.grid(user_input_exclude_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          exclude_covariates=c(exclude_covariates,temp)
        }else{ #otherwise it is time invariant
          exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c])
        }
      }
    }
  }

  #assesses and reformat use input of keep covariates
  user_input_keep_covariates=keep_covariates
  keep_covariates=NULL
  if (length(user_input_keep_covariates)>0){
    for (c in 1:length(user_input_keep_covariates)){
      if (grepl("\\.",user_input_keep_covariates[c])){ #period indicates user specified something time-varying
        keep_covariates=c(keep_covariates, user_input_keep_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_keep_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified

          temp=apply(expand.grid(user_input_keep_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          keep_covariates=c(keep_covariates,temp)
        }else{ #otherwise it is time invariant
          keep_covariates=c(keep_covariates, user_input_keep_covariates[c])
        }
      }
    }
  }




  #creates final dataset with only relevant variables
  covariates_to_include=unique(unlist(lapply(covariates_to_include, function(x) c(x$row, x$column))))

  variables_to_include=unique(c(ID, "WAVE", exposures, outcomes, covariates_to_include, keep_covariates))
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
