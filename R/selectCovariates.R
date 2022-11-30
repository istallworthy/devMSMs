#' Identifies the available covariates
#'
#' Identifies the covariates in in the dataset that will be candidates for confounding
#'
#' identifyCovariates
#' @param data data from formatDataStruct
#' @param object msm object that contains all relevant user inputs
#' @return potential_covariates
#' @export
#' @seealso [formatDataStruct()]
#' @examples identifyCovariates(object,data, exclude_covariates)
#'
identifyCovariates <- function(object, data){

  ID=object$ID
  exposures=object$exposures
  outcomes=object$outcomes
  exclude_covariates=object$exclude_covariates

  #identifying potential covariates (time-varying and time invariant)
  potential_covariates=colnames(data)[colnames(data) %in% (c(ID, "WAVE", exposures, outcomes, exclude_covariates))==FALSE]

  cat(paste0("Below are the ", as.character(length(potential_covariates)), " covariates that will be considered as potential confounding variables. If you would like to exclude any variables, please add them to 'exclude_covariates' when creating the msms object and re-run."),"\n")
  print(potential_covariates)

  return(potential_covariates)
}





#' Identify potential confounds
#'
#' Identifies potential confounds from the covariates based on those that are correlated with either exposures or outcomes above 0.1
#'
#' @param object msm object that contains all relevant user inputs
#' @return covariates_to_include
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom dplyr filter
#' @importFrom stargazer stargazer
#' @seealso [formatDataStruct()], [makeTimePtDatasets()]
#' @examples identifyPotentialConfounds(home_dir, data, time_pt_datasets, time_pts, exclude_covariates=NULL)
#'
identifyPotentialConfounds <- function(object){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  time_pts=object$time_pts
  exposure_time_pts=object$exposure_time_pts
  outcome_time_pt=object$outcome_time_pt
  time_varying_covariates=object$time_varying_variables
  exclude_covariates=object$exclude_covariates
  balance_thresh=object$balance_thresh
  mandatory_keep_covariates=object$mandatory_keep_covariates



  #formatting any covariates to EXCLUDE
  user_input_exclude_covariates=exclude_covariates
  exclude_covariates=NULL
  if (length(user_input_exclude_covariates)>0){
    for (c in 1:length(user_input_exclude_covariates)){
      if (grepl("\\.",user_input_exclude_covariates[c])){ #period indicates user specified something time-varying at specific time point
        exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_exclude_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified, add times

          temp=apply(expand.grid(user_input_exclude_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          exclude_covariates=c(exclude_covariates,temp)
        }else{ #otherwise it is time invariant
          exclude_covariates=c(exclude_covariates, user_input_exclude_covariates[c])
        }
      }
    }
  }


  #formatting any covariates to INCLUDE
  user_input_include_covariates=mandatory_keep_covariates
  include_covariates=NULL
  if (length(user_input_include_covariates)>0){
    for (c in 1:length(user_input_include_covariates)){
      if (grepl("\\.",user_input_include_covariates[c])){ #period indicates user specified something time-varying at specific time point
        include_covariates=c(include_covariates, user_input_include_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_include_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified, add times

          temp=apply(expand.grid(user_input_include_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          include_covariates=c(include_covariates,temp)
        }else{ #otherwise it is time invariant
          include_covariates=c(include_covariates, user_input_include_covariates[c])
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


  #makes wide dataset
  data_wide=suppressWarnings(stats::reshape(data=data,
                                           idvar=ID,
                                           v.names= time_varying_covariates, #list ALL time-varying covariates
                                           timevar="WAVE",
                                           times=c(time_pts),
                                           direction="wide"))

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

        temp_corr={}

        #gathers correlations with exposure from variables at exposure time point or below
        time_pt=exposure_time_pts[x]
        data_long=data%>%dplyr::filter(WAVE<time_pt+1) #filters to current time point or prior
        d=suppressWarnings(stats::reshape(data=data_long,
                                          idvar=ID,
                                          v.names= time_varying_covariates, #list ALL time-varying covariates
                                          timevar="WAVE",
                                          times=c(time_pts),
                                          direction="wide"))
        d=d[,2:ncol(d)]
        d=as.data.frame(lapply(d, as.numeric))
        c <- suppressWarnings(Hmisc::rcorr(as.matrix(d))) #makes corr table of all vars; cannot have dates or non-factored characters
        c=flattenCorrMatrix(c$r, c$P)
        #denotes time varying variables
        # c$row[c$row %in% time_varying_covariates]=paste0(c$row[c$row %in% time_varying_covariates], ".", time_pt)
        # c$column[c$column %in% time_varying_covariates]=paste0(c$column[c$column %in% time_varying_covariates], ".", time_pt)
        keep=c[c$row %in% include_covariates | c$column %in% include_covariates, ]
        keep=keep[grepl(paste0(exposure, ".", time_pt),keep$row) | grepl(paste0(exposure, ".", time_pt),keep$column),]
        keep$exp_time=time_pt
        filtered=c%>%
          filter(abs(c$cor)>balance_thresh) #accounts for d=0.2 see Stuart paper for rationale
        filtered$exp_time=time_pt
        filtered=filtered[filtered$row %in% paste0(exposure, ".", time_pt) | filtered$column %in% paste0(exposure, ".", time_pt),] #finds those that include exposure
                          # | filtered$column %in% paste0(outcome, ".", outcome_time_pt) | filtered$row %in% paste0(outcome, ".", outcome_time_pt),] #gets only those associated with exposure at given time pt OR outcome at outcome time pt
        filtered=filtered[!filtered$row %in% exclude_covariates[grepl(exposure, exclude_covariates)==F] & !filtered$column %in% exclude_covariates[grepl(exposure, exclude_covariates)==F],] #rejects those the user wants to exclude (as long as it's not the exposure itself)
        temp_corr=rbind(temp_corr, filtered)
        filtered=NULL
        c=NULL
        data_long=NULL

        #makes sure keep_covars are included in table so that table contents match forms
        if (sum(temp_corr$row %in% include_covariates | temp_corr$row %in% include_covariates)< nrow(keep)){
          temp_corr=rbind(temp_corr, keep)
        }

        #gathers correlations between variables at time pt or prior with outcome at outcome time point
        data_long=data%>%dplyr::filter(WAVE<time_pt+1 |WAVE==outcome_time_pt) #filters to current time point or prior
        o=suppressWarnings(stats::reshape(data=data_long,
                                         idvar=ID,
                                         v.names= time_varying_covariates, #list ALL time-varying covariates
                                         timevar="WAVE",
                                         times=c(time_pts),
                                         direction="wide"))
        o=o[,2:ncol(o)]
        o=as.data.frame(lapply(o, as.numeric))
        o <- suppressWarnings(Hmisc::rcorr(as.matrix(o))) #makes corr table of all vars; cannot have dates or non-factored characters
        o=flattenCorrMatrix(o$r, o$P)
        filtered=o%>%
          filter(abs(o$cor)>balance_thresh) #accounts for d=0.2 see Stuart paper for rationale
        filtered$exp_time=time_pt
        filtered=filtered[filtered$row %in% paste0(outcome, ".", outcome_time_pt) | filtered$column %in% paste0(outcome, ".", outcome_time_pt),] #gets only those associated with outcome at outcome time point
        filtered=filtered[as.numeric(sapply(strsplit(filtered$row, "\\."), "[", 2))<time_pt+1 | as.numeric(sapply(strsplit(filtered$column, "\\."), "[", 2))<time_pt+1,] #making sure it includes a variable from exposure time pt or prior
        filtered=filtered[!filtered$row %in% exclude_covariates[grepl(exposure, exclude_covariates)==F] & !filtered$column %in% exclude_covariates[grepl(exposure, exclude_covariates)==F],] #rejects those the user wants to exclude
        temp_corr=rbind(temp_corr, filtered)
        filtered=NULL
        c=NULL
        data_long=NULL



        temp_corr=temp_corr[order(temp_corr$row, temp_corr$column),] #orders by covariate

        covariate_correlations=rbind(covariate_correlations, temp_corr)


      }



      #save out correlations
      suppressMessages(stargazer::stargazer(covariate_correlations,type="html", digits=2, column.labels = colnames(covariate_correlations),summary=FALSE, rownames = FALSE, header=FALSE,
                           out=paste0(home_dir, "balance/potential confounds/", exposure, "-", outcome, "_potential_counfound_correlations.html"), show))


      cat("\n")
      # write.csv(covariate_correlations, paste0(home_dir, "balance/covariate_correlations.csv"))
      cat(paste0("Check the 'balance/potential confounds/' folder to view an html file of a table of the correlations for the potential confounds for ", exposure, "-", outcome),"\n")

      covariates_to_include[[paste0(exposure, "-", outcome, sep="")]] <-covariate_correlations

    }
  }

  cat("\n")
  cat(paste0("A total of ", as.character(length(unique(unlist(lapply(covariates_to_include, function(x) c(x$row, x$column)))))), " potential confounds will be considered for balancing."),"\n")

  return(covariates_to_include)
}



#' Creates data structure for imputation
#'
#' Creates data structure for imputation with only the variables needed (i.e., exposures, outcomes, potential confounds)
#'
#' @param object msm object that contains all relevant user inputs
#' @param covariates_to_include output from identifyPotentialConfounds
#' @importFrom Hmisc rcorr
#' @importFrom knitr kable
#' @importFrom dplyr filter
#' @return data_to_impute
#' @export
#' @seealso [formatDataStruct()], [identifyPotenialConfounds()]
#' @examples dataToImpute(ID, data, exposures, outcomes, covariates_to_include, exclude_covariates)
#'
dataToImpute <-function(object, covariates_to_include){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  time_pts=object$time_pts
  exclude_covariates=object$exclude_covariates
  mandatory_keep_covariates=object$mandatory_keep_covariates
  time_varying_covariates=object$time_varying_variables

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
  user_input_mandatory_keep_covariates=mandatory_keep_covariates
  mandatory_keep_covariates=NULL
  if (length(user_input_mandatory_keep_covariates)>0){
    for (c in 1:length(user_input_mandatory_keep_covariates)){
      if (grepl("\\.",user_input_mandatory_keep_covariates[c])){ #period indicates user specified something time-varying
        mandatory_keep_covariates=c(mandatory_keep_covariates, user_input_mandatory_keep_covariates[c]) #can just be added as is
      }else{
        if (sum(grepl(user_input_mandatory_keep_covariates[c], time_varying_covariates))>0){ #if it is time-varying but no time specified

          temp=apply(expand.grid(user_input_mandatory_keep_covariates[c], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".") #appends all time points
          mandatory_keep_covariates=c(mandatory_keep_covariates,temp)
        }else{ #otherwise it is time invariant
          mandatory_keep_covariates=c(mandatory_keep_covariates, user_input_mandatory_keep_covariates[c])
        }
      }
    }
  }




  #creates final dataset with only relevant variables
  covariates_to_include=unique(unlist(lapply(covariates_to_include, function(x) c(x$row, x$column))))
  covariates_to_include=covariates_to_include[order(covariates_to_include)]

  variables_to_include=unique(c(ID, "WAVE", exposures, outcomes, covariates_to_include, mandatory_keep_covariates, time_varying_covariates))
  data2=as.data.frame(data[names(data)[names(data) %in% variables_to_include] ])
  # data2=data2[,!colnames(data2) %in% c(exclude_covariates)]

  data_to_impute=data2

  data2=data2[,!colnames(data2) %in% c(exposures, outcomes)]

  #inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  hi_corr_covars <- suppressWarnings(Hmisc::rcorr(as.matrix(data2[,3:ncol(data2)])))

  hi_corr_covars=flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)
  View_hi_corr_covars=hi_corr_covars%>%
    dplyr::filter(abs(hi_corr_covars$cor)>0.6)

  cat("USER ALERT: To simplify the balancing models, consider removing any highly correlated, redundant covariates by listing them in the 'exclude_covariates' field in the msm object and re-running:")
  print(knitr::kable(View_hi_corr_covars))

  write.csv(data_to_impute, paste0(home_dir, "imputations/data_to_impute.csv"))
  cat("See the 'imputations' folder for a csv file of the data to be imputed","\n")

  return(data_to_impute)
}
