#' Identifies the available covariates
#'
#' Identifies the covariates in in the dataset that will be candidates for confounding
#'
#' identifyCovariates
#' @param data data from formatDataStruct
#' @param object msm object that contains all relevant user inputs
#' @param exclude_covariates list of variables to exclude based on theoretical or practical reasons
#' @return potential_covariates
#' @export
#' @seealso [formatDataStruct()]
#' @examples identifyCovariates(object,data, exclude_covariates)
#'
identifyCovariates <- function(object, data, exclude_covariates=NULL){

  ID=object$ID
  exposures=object$exposures
  outcomes=object$outcomes

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
#' @param object msm object that contains all relevant user inputs
#' @param time_pt_datasets output from makeTimePtDatasets
#' @param exclude_covariates list of variables to exclude based on theoretical or practical reasons
#' @param time_varying_covariates list of covariates that are time-varying
#' @return covariates_to_include
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom stargazer stargazer
#' @seealso [formatDataStruct()], [makeTimePtDatasets()]
#' @examples identifyPotentialConfounds(home_dir, data, time_pt_datasets, time_pts, exclude_covariates=NULL)
#'
identifyPotentialConfounds <- function(object, time_pt_datasets, time_varying_covariates=NULL, exclude_covariates=NULL){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  time_pts=object$time_pts
  exposure_time_pts=object$exposure_time_pts
  outcome_time_pts=object$outcome_time_pts

  #formatting any covariates to exclude
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


  covariates_to_include=list()

  #cycles through outcomes
  for (z in seq(length(outcomes))){
    outcome=outcomes[z]

    #cycles through exposures
    for (y in 1:length(exposures)){
      exposure=exposures[y]

      #Gathers all confounds relevant to an exposure at a given time point (note: lagged values of exposure are added later in the forms function)

      #gathers correlations with exposure at given exposure time point
      covariate_correlations={}
      for (x in 1:length(exposure_time_pts)){
        time_pt=exposure_time_pts[x]
        d=time_pt_datasets[[x]]
        d=d[,3:ncol(d)]
        d=as.data.frame(lapply(d, as.numeric))
        c <- suppressWarnings(Hmisc::rcorr(as.matrix(d))) #makes corr table of all vars; cannot have dates or non-factored characters
        c=flattenCorrMatrix(c$r, c$P)
        filtered=c%>%
          filter(abs(c$cor)>0.1) #accounts for d=0.2 see Stuart paper for rationale
        filtered$exp_time=time_pt
        filtered=filtered[filtered$row %in% exposure | filtered$column %in% exposure,] #gets only those associated with exposure at given time pt
        covariate_correlations=rbind(covariate_correlations, filtered)
        filtered=NULL
        c=NULL

        #gathers correlations with outcome at any outcome time point
        o=data[data$WAVE %in% outcome_time_pts,]
        o=o[,3:ncol(o)]
        o=as.data.frame(lapply(o, as.numeric))
        o <- suppressWarnings(Hmisc::rcorr(as.matrix(o))) #makes corr table of all vars; cannot have dates or non-factored characters
        o=flattenCorrMatrix(o$r, o$P)
        filtered=o%>%
          filter(abs(o$cor)>0.1) #accounts for d=0.2 see Stuart paper for rationale
        filtered$exp_time=time_pt
        filtered=filtered[filtered$row %in% outcome |filtered$column %in% outcome,] #gets only those associated with exposure at given time pt
        covariate_correlations=rbind(covariate_correlations, filtered)
        filtered=NULL
        c=NULL
      }



      # #keeps only those that involve either an exposure or outcome
      # # test2=subset(covariate_correlations, subset=grepl(exposures, row) | grepl(exposures, column) | grepl(outcomes, row) | grepl(outcomes, column))
      # covariate_correlations=covariate_correlations[covariate_correlations$row %in% c(exposures, outcomes) | covariate_correlations$column %in% c(exposures,outcomes),]
      # covariate_correlations=covariate_correlations[!covariate_correlations$row %in% c(exclude_covariates) & !covariate_correlations$column %in% c(exclude_covariates),]

      #save out correlations
      stargazer::stargazer(covariate_correlations,type="html", digits=2, column.labels = colnames(covariate_correlations),summary=FALSE, rownames = FALSE, header=FALSE,
                           out=paste0(home_dir, "balance/potential confounds/", exposure, "-", outcome, "_potential_counfound_correlations.html"))
      # write.csv(covariate_correlations, paste0(home_dir, "balance/covariate_correlations.csv"))
      print(paste0("Check the 'balance/potential confounds/' folder to view an html file of a table of the correlations for the potential confounds for ", exposure, "-", outcome))

      covariates_to_include[[paste0(exposure, "-", outcome, sep="")]] <-covariate_correlations


      #
      #   variables_to_include=c(variables_to_include, )
      #
      #   #for each exposure-outcome pair, list correlated covariates that include only that exposure and/or outcomes of interest
      #   variables_to_include=c(ID, "WAVE", exposures, outcomes)
      #
      #       variables=covariate_correlations
      #       unique_variables=unique(c(variables$row, variables$column)) #finds variables correlated at any time point
      #       unique_variables=unique_variables[!unique_variables %in% c(user_input_exclude_covariates, "WAVE", ID)]
      #       variables_to_include=c(variables_to_include, unique_variables)
      #       # assign(paste(exposures[y], "_covars_to_include", sep=""), unique_variables) #makes data frames for each tx with only relevant covariates
      #       # covariates_to_include[[paste(exposures[y], "_covars_to_include", sep="")]] <- unique_variables #saves to list
      #       covariates_to_include[[paste(exposures[y],"-", outcomes[z], "_covars_to_include", sep="")]] <- variables #saves relevant covariates at different time points
      #
      #       #save out correlations by exposure-outcome pairing
      #       stargazer::stargazer(covariate_correlations,type="html", digits=2, column.labels = colnames(covariate_correlations),summary=FALSE, rownames = FALSE, header=FALSE,
      #                            out=paste0(home_dir, "balance/", exposures[y],"-", outcomes[z], "potential_counfound_correlations.html"))
      #       print(paste0("Check the 'balance' folder to view an html file of a table of all correlations between potential confounds for ", exposures[y],"-", outcomes[z], " at > 0.1 for each time point"))

    }
  }

  print(paste0("A total of ", as.character(length(unique(unlist(lapply(covariates_to_include, function(x) c(x$row, x$column)))))), " potential confounds will be considered for balancing."))

  return(covariates_to_include)
}



#' Creates data structure for imputation
#'
#' Creates data structure for imputation with only the variables needed (i.e., exposures, outcomes, potential confounds)
#'
#' @param object msm object that contains all relevant user inputs
#' @param covariates_to_include output from identifyPotentialConfounds
#' @param exclude_covariates list of variables to exclude based on theoretical or practical reasons
#' @param keep_covariates list of covariates to keep
#' @importFrom Hmisc rcorr
#' @importFrom knitr kable
#' @return data_to_impute
#' @export
#' @seealso [formatDataStruct()], [identifyPotenialConfounds()]
#' @examples dataToImpute(ID, data, exposures, outcomes, covariates_to_include, exclude_covariates)
#'
dataToImpute <-function(object, covariates_to_include, exclude_covariates=NULL, keep_covariates=NULL){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  time_pts=object$time_pts

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
  covariates_to_include=covariates_to_include[order(covariates_to_include)]

  variables_to_include=unique(c(ID, "WAVE", exposures, outcomes, covariates_to_include, keep_covariates))
  data2=as.data.frame(data[names(data)[names(data) %in% variables_to_include] ])
  data2=data2[,!colnames(data2) %in% c(exclude_covariates)]

  data_to_impute=data2

  data2=data2[,!colnames(data2) %in% c(exposures, outcomes)]

  #inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  hi_corr_covars <- suppressWarnings(Hmisc::rcorr(as.matrix(data2[,3:ncol(data2)])))

  hi_corr_covars=flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)
  View_hi_corr_covars=hi_corr_covars%>%
    filter(abs(hi_corr_covars$cor)>0.6)

  print("USER ALERT: To simplify the balancing models consider removing any highly correlated, redundant covariates by listing them in the 'exclude_covariates' field above and re-running this function:")
  print(knitr::kable(View_hi_corr_covars))

  write.csv(data_to_impute, paste0(home_dir, "imputations/data_to_impute.csv"))
  print("See the 'imputations' folder for a csv file of the data to be imputed")

  return(data_to_impute)
}
