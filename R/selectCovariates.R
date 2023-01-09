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
  time_varying_covariates=object$time_varying_variables
  time_pts=object$time_pts
  time_var_exclude=object$time_var_exclude

  #gathering and delineating exclude covariates
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

  # browser()


  #identifying potential covariates (time-varying and time invariant)
  potential_covariates=colnames(data)[colnames(data) %in% (c(ID, "WAVE"))==FALSE]

  if (length(exposures)==1){ #remove exposures from list if only one
    potential_covariates=colnames(data)[colnames(data) %in% exposures==FALSE]

  }

  # if(all(is.na(data%>%dplyr::filter(WAVE<outcome_time_pt)%>%dplyr::select(outcomes)))){
  #   potential_covariates=colnames(data)[colnames(data) %in% outcome)==FALSE]
  # }

  time_invar_covars=potential_covariates[!potential_covariates %in% time_varying_covariates]

  time_var_covars=apply(expand.grid(potential_covariates[potential_covariates %in% time_varying_covariates], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".")

  all_potential_covariates=c(time_invar_covars, time_var_covars)

  #exclusions
  all_potential_covariates=all_potential_covariates[!all_potential_covariates %in% c(exclude_covariates, paste(outcomes, outcome_time_pt, sep="."), time_var_exclude)]

  all_potential_covariates=all_potential_covariates[order(all_potential_covariates)]

  #formatting for table output to visualize available covariates by time point
  covar_table=data.frame(variable=sapply(strsplit(all_potential_covariates, "\\."), "[", 1),
                         time_pt=sapply(strsplit(all_potential_covariates, "\\."), "[", 2)
  )
  covar_table=covar_table[order(covar_table$time_pt,  covar_table$variable),]
  covar_table=aggregate(variable ~ time_pt, covar_table, toString)
  covar_table=covar_table[order(as.numeric(covar_table$time_pt)),]
  write.csv(covar_table, paste0(home_dir, "balance/covariates_considered_by_time_pt.csv"), row.names = F)

  unique_vars=length(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))

  test=as.data.frame(matrix(nrow=length(time_pts), ncol=unique_vars))
  colnames(test)=unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1)))[order(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))]
  rownames(test)=time_pts
  for (l in 1:nrow(test)){
    test[l, c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]), all_potential_covariates)], "\\."), "[",1), time_invar_covars)]= 1
  }
  NumTimePts=data.frame(NumTimePts=colSums(test, na.rm=T))
  test=rbind(test, t(NumTimePts))
  NumVars=data.frame(NumVars=rowSums(test, na.rm=T))
  test[1:nrow(test),ncol(test)+1]=NumVars
  write.csv(test, paste0(home_dir, "balance/matrix_of_covariates_considered_by_time_pt.csv"), row.names = T)

  cat("See the balance folder for a table and matrix displaying all covariates considered for each time point", "\n")



  cat(paste0("Below are the ", as.character(length(all_potential_covariates)), " variables, across ", unique_vars, " unique constructs, that will be evalated empirically as potential confounding variables for all exposure-outcome pairs."), "\n",
      "Please inspect this list carefully. It should include all time-varying covariates (excluding time points when they were not collected), time invariant covariates, exposures (if you have listed more than one), and outcome variables if they were collected at time points earlier than the outcome time point.", "\n",
      "If you would like to exclude any variables, please add them to 'exclude_covariates' field in the msmObject and re-run.","\n")
  cat("\n")
  print(all_potential_covariates)

  return(all_potential_covariates)
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
identifyPotentialConfounds <- function(object, all_potential_covariates){


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
  bal_only_exp=object$bal_only_exp


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

  if (sum(time_varying_covariates %in% colnames(data))!=length(time_varying_covariates)){
    stop("Please make sure all time_varying_covariates entered in the msmObject match column names in your dataset")
  }

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
      # keep=NA
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
        d=d[,colnames(d)!=ID]
        d=d[,colnames(d) %in% c(paste(exposure, exposure_time_pts[exposure_time_pts<time_pt+1], sep="."),
                                paste(outcome, outcome_time_pt, sep="."),
                                all_potential_covariates)] #making sure only the right variables are being considered (what user inspected)
        d=as.data.frame(lapply(d, as.numeric))
        c <- suppressWarnings(Hmisc::rcorr(as.matrix(d))) #makes corr table of all vars; cannot have dates or non-factored characters
        c=flattenCorrMatrix(c$r, c$P)
        keep=c[c$row %in% include_covariates | c$column %in% include_covariates, ]
        keep=keep[grepl(paste0(exposure, ".", time_pt),keep$row) | grepl(paste0(exposure, ".", time_pt),keep$column),]
        if (nrow(keep)>0){keep$exp_time=time_pt}
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

        # browser()

        if (bal_only_exp==T){ #if the user wants to define confounders only based on exposure
          temp_corr=temp_corr[!temp_corr %in% paste0(outcome, ".", outcome_time_pt)]
          cat("Per user specification, only the exposure(s) will be considered in empirically identifying confounders.", "\n")

        }else{ #otherwise, consider correlations w/ outcome

          #gathers correlations between variables at time pt or prior with outcome at outcome time point
          data_long=data%>%dplyr::filter(WAVE<time_pt+1 |WAVE==outcome_time_pt) #filters to current time point or prior
          o=suppressWarnings(stats::reshape(data=data_long,
                                            idvar=ID,
                                            v.names= time_varying_covariates, #list ALL time-varying covariates
                                            timevar="WAVE",
                                            times=c(time_pts),
                                            direction="wide"))

          o=o[,colnames(o)!=ID]
          o=o[,colnames(o) %in% c(paste(exposure, exposure_time_pts[exposure_time_pts<time_pt+1], sep="."),
                                  paste(outcome, outcome_time_pt, sep="."),
                                  all_potential_covariates)] #make sure only using correct variables as specified by user
          o=as.data.frame(lapply(o, as.numeric))
          o[,colnames(o) %in% factor_covariates]=as.data.frame(lapply(o[colnames(o) %in% factor_covariates], as.factor))

          o <- suppressWarnings(Hmisc::rcorr(as.matrix(o))) #makes corr table of all vars; cannot have dates or non-factored characters
          o=flattenCorrMatrix(o$r, o$P)
          filtered=o%>%
            filter(abs(o$cor)>balance_thresh) #accounts for d=0.2 see Stuart paper for rationale
          filtered$exp_time=time_pt
          filtered=filtered[filtered$row %in% paste0(outcome, ".", outcome_time_pt) | filtered$column %in% paste0(outcome, ".", outcome_time_pt),] #gets only those associated with outcome at outcome time point
          # filtered=filtered[filtered$row %in% outcome | filtered$column %in% outcome,] #gets only those associated with outcome at outcome time point
          filtered=filtered[as.numeric(sapply(strsplit(filtered$row, "\\."), "[", 2))<time_pt+1 | as.numeric(sapply(strsplit(filtered$column, "\\."), "[", 2))<time_pt+1 |
                              filtered$row %in% time_invar_covars | filtered$column %in% time_invar_covars
                            ,] #making sure it includes a variable from exposure time pt or prior or is time invariant
          filtered=filtered[!filtered$row %in% exclude_covariates[grepl(exposure, exclude_covariates)==F] & !filtered$column %in% exclude_covariates[grepl(exposure, exclude_covariates)==F],] #rejects those the user wants to exclude
          filtered=na.omit(filtered)
          temp_corr=rbind(temp_corr, filtered)
          filtered=NULL
          c=NULL
          data_long=NULL
        }



        temp_corr=temp_corr[order(temp_corr$row, temp_corr$column),] #orders by covariate
        # temp_corr=temp_corr[temp_corr$row %in% c(paste(exposure, exposure_time_pts, sep=".")
        #                                          temp_corr$row[grepl(exposure, temp_corr$row)==F][temp_corr$row[grepl(exposure, temp_corr$row)==F] %in% all_potential_covariates]
        #                                          )]

        covariate_correlations=rbind(covariate_correlations, temp_corr)


      }



      #save out correlations
      sink(paste0(home_dir, "balance/potential confounds/", exposure, "-", outcome, "_potential_counfound_correlations.html"))
      stargazer::stargazer(covariate_correlations,type="html", digits=2, column.labels = colnames(covariate_correlations),summary=FALSE, rownames = FALSE, header=FALSE,
                           out=paste0(home_dir, "balance/potential confounds/", exposure, "-", outcome, "_potential_counfound_correlations.html"))
      sink()

      cat("\n")
      # write.csv(covariate_correlations, paste0(home_dir, "balance/covariate_correlations.csv"))
      cat(paste0("USER ALERT: Check the 'balance/potential confounds/' folder to view an html file of a table of the correlations for the potential confounders for ", exposure, "-", outcome),"\n")

      covariates_to_include[[paste0(exposure, "-", outcome, sep="")]] <-covariate_correlations

    }
  }

  cat("\n")
  cat(paste0("A total of ", as.character(length(unique(unlist(lapply(covariates_to_include, function(x) c(x$row, x$column)))))), " covariate confounders will be considered for balancing across all exposure-outcome pairs and all time points."),"\n")

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
#' @importFrom corrplot corrplot
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
  data2=data2[,!colnames(data2) %in% c(exclude_covariates[!exclude_covariates %in% c(exposures, outcomes)])] #makes sure not to exclude exposures or outcomes even if listed to exclude

  data_to_impute=data2

  #makes correlation table
  pdf(file = paste0(home_dir, "all_vars_corr_plot.pdf"))
  suppressWarnings(corrplot::corrplot(cor(
    as.data.frame(lapply(data_to_impute[,colnames(data_to_impute)[colnames(data_to_impute)!=ID]], as.numeric)), use="pairwise.complete.obs"
  ), method="color", order='alphabet', diag=FALSE, type="lower", tl.cex = 0.5, tl.col="black"))
  dev.off()

  cat("A correlation plot of all variables has been saved in the home directory", "\n")


  data2=data2[,!colnames(data2) %in% c(exposures, outcomes)]

  #inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  hi_corr_covars <- suppressWarnings(Hmisc::rcorr(as.matrix(data2[,3:ncol(data2)])))

  # browser()

  hi_corr_covars=flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)
  View_hi_corr_covars=hi_corr_covars%>%
    dplyr::filter(abs(hi_corr_covars$cor)>0.7)

  View_hi_corr_covars=View_hi_corr_covars[!View_hi_corr_covars$row %in% c("WAVE") & !View_hi_corr_covars$column %in% c("WAVE"),]

  cat("USER ALERT: To simplify the balancing models, consider removing any highly correlated, redundant covariates by listing them in the 'exclude_covariates' field in the msm object and re-running:", "\n")
  print(knitr::kable(View_hi_corr_covars))

  write.csv(data_to_impute, paste0(home_dir, "imputations/data_to_impute.csv"))
  cat("See the 'imputations' folder for a csv file of the data to be imputed","\n")

  return(data_to_impute)
}
