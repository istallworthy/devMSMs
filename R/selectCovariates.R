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
  exposure=object$exposure
  outcome=object$outcome
  exclude_covariates=object$exclude_covariates
  time_varying_covariates=object$time_varying_variables
  time_pts=object$time_pts
  time_var_exclude=object$time_var_exclude
  outcome_time_pt=object$outcome_time_pt
  home_dir=object$home_dir

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

  if (length(exposure)==1){ #remove exposure from list if only one
    potential_covariates=colnames(data)[colnames(data) %in% exposure==FALSE]

  }


  time_invar_covars=potential_covariates[!potential_covariates %in% time_varying_covariates]

  time_var_covars=apply(expand.grid(potential_covariates[potential_covariates %in% time_varying_covariates], as.character(as.numeric(time_pts))), 1, paste0, sep="", collapse=".")

  all_potential_covariates=c(time_invar_covars, time_var_covars)

  #exclusions
  all_potential_covariates=all_potential_covariates[!all_potential_covariates %in% c(exclude_covariates, paste(outcome, outcome_time_pt, sep="."), time_var_exclude)]

  all_potential_covariates=all_potential_covariates[order(all_potential_covariates)]

  #formatting for table output to visualize available covariates by time point
  covar_table=data.frame(variable=sapply(strsplit(all_potential_covariates, "\\."), "[", 1),
                         time_pt=sapply(strsplit(all_potential_covariates, "\\."), "[", 2)
  )
  covar_table=covar_table[order(covar_table$time_pt,  covar_table$variable),]
  covar_table=aggregate(variable ~ time_pt, covar_table, toString)
  covar_table=covar_table[order(as.numeric(covar_table$time_pt)),]
  write.csv(covar_table, paste0(home_dir, "balance/", exposure, "-", outcome, "_covariates_considered_by_time_pt.csv"), row.names = F)

  unique_vars=length(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))

  test=as.data.frame(matrix(nrow=length(time_pts), ncol=unique_vars))
  colnames(test)=unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1)))[order(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))]
  rownames(test)=time_pts
  for (l in 1:nrow(test)){
    test[l, c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]), all_potential_covariates)], "\\."), "[",1), time_invar_covars)]= 1
  }
  test=test[,colnames(test)[!colnames(test) %in% c(ID, "WAVE")]]
  NumTimePts=data.frame(NumTimePts=colSums(test, na.rm=T))
  test=rbind(test, t(NumTimePts))
  NumVars=data.frame(NumVars=rowSums(test, na.rm=T))
  test[1:nrow(test),ncol(test)+1]=NumVars
  write.csv(test, paste0(home_dir, "balance/",  exposure, "-", outcome, "_matrix_of_covariates_considered_by_time_pt.csv"), row.names = T)

  cat("See the balance folder for a table and matrix displaying all covariates considered for each time point.", "\n")



  #-2 to exclude ID and WAVE
  cat(paste0("Below are the ", as.character(length(all_potential_covariates)-2), " variables, across ", unique_vars-2, " unique constructs, that will be evalated empirically as potential confounding variables for all exposure-outcome pairs."), "\n",
      "Please inspect this list carefully. It should include all time-varying covariates (excluding time points when they were not collected), time invariant covariates, exposure (if you have listed more than one), and outcome variables if they were collected at time points earlier than the outcome time point.", "\n",
      "If you would like to exclude any variables from consideration as covariate confounders, please add them to 'exclude_covariates' field in the msmObject and re-run.","\n")
  cat("\n")
  print(all_potential_covariates[!all_potential_covariates %in% c(ID, "WAVE")])

  return(all_potential_covariates)
}




#' Creates data structure for imputation
#'
#' Creates data structure for imputation with only the variables needed (i.e., exposure, outcome, potential confounds)
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
#' @examples dataToImpute(ID, data, exposure, outcome, covariates_to_include, exclude_covariates)
#'
dataToImpute <-function(object, all_potential_covariates){

  ID=object$ID
  home_dir=object$home_dir
  exposure=object$exposure
  outcome=object$outcome
  time_pts=object$time_pts
  exclude_covariates=object$exclude_covariates
  mandatory_keep_covariates=object$mandatory_keep_covariates
  time_varying_covariates=object$time_varying_variables
  exposure_epochs=object$exposure_epochs

  covariates_to_include=all_potential_covariates

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
  covariates_to_include=covariates_to_include[order(covariates_to_include)]

  variables_to_include=unique(c(ID, "WAVE", exposure, outcome, covariates_to_include, mandatory_keep_covariates, time_varying_covariates))
  data2=as.data.frame(data[names(data)[names(data) %in% variables_to_include] ])
  data2=data2[,!colnames(data2) %in% c(exclude_covariates[!exclude_covariates %in% c(exposure, outcome)])] #makes sure not to exclude exposure or outcome even if listed to exclude

  data_to_impute=data2

  #makes correlation table
  pdf(file = paste0(home_dir,  exposure, "-", outcome,"_all_vars_corr_plot.pdf"))
  suppressWarnings(corrplot::corrplot(cor(
    as.data.frame(lapply(data_to_impute[,colnames(data_to_impute)[colnames(data_to_impute)!=ID]], as.numeric)), use="pairwise.complete.obs"
  ), method="color", order='alphabet', diag=FALSE, type="lower", tl.cex = 0.5, tl.col="black"))
  dev.off()

  cat("A correlation plot of all variables has been saved in the home directory", "\n")


  data2=data2[,!colnames(data2) %in% c(exposure, outcome)]

  #inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  hi_corr_covars <- suppressWarnings(Hmisc::rcorr(as.matrix(data2[,3:ncol(data2)])))


  hi_corr_covars=flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)
  View_hi_corr_covars=hi_corr_covars%>%
    dplyr::filter(abs(hi_corr_covars$cor)>0.7)

  View_hi_corr_covars=View_hi_corr_covars[!View_hi_corr_covars$row %in% c("WAVE") & !View_hi_corr_covars$column %in% c("WAVE"),]

  cat("USER ALERT: To simplify the balancing models, consider removing any highly correlated, redundant covariates by listing them in the 'exclude_covariates' field in the msm object and re-running:", "\n")
  # print(knitr::kable(View_hi_corr_covars))
  cat(knitr::kable(View_hi_corr_covars, caption="Correlated Covariates", format='pipe'),  sep="\n")


 #  #adding in exposure epoch interactions pre-imputation
 #  time_varying_wide=apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep="", collapse=".")
 #  time_varying_wide=sort(time_varying_wide)
 #  time_varying_wide=c(ID, time_varying_wide)
 #  dat_wide=suppressWarnings(stats::reshape(data=data_to_impute,idvar=ID,v.names= time_varying_covariates, timevar="WAVE", times=c(time_pts), direction="wide"))
 #  exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")
 #  interactions<-paste0(unlist(lapply(2:length(exp_epochs), function(z){
 #    apply(combn(exp_epochs,z), 2, paste, sep="", collapse=":")
 #  })))
 #  new=data.frame(ID=dat_wide[,ID])
 #  names(new)=ID
 #  #calculates the mean value for each exposure for each exposure epoch
 #  for (e in 1:nrow(exposure_epochs)){
 #    epoch=exposure_epochs[e,1]
 #    temp=data.frame(row.names=1:nrow(dat_wide))
 #    new_var=paste0(exposure, "_", epoch)
 #    #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
 #    for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
 #      level=as.numeric(unlist(exposure_epochs[e,2]))[l]
 #      z=dat_wide[,which(grepl(paste0(exposure, ".", level), names(dat_wide)))]
 #      temp=cbind(temp, z)
 #    }
 #    #adds a new variable of the exposure averaged within epoch
 #    new=new%>%dplyr::mutate(!!new_var :=rowMeans(temp, na.rm=T))
 #    new[,new_var]=as.numeric(new[,new_var])
 #  }
 #  new[new=="NaN"]=NA
 #
 #  for (x in interactions){
 #    temp=new%>%dplyr::select(c(unlist(strsplit(x, ":"))))
 #    new=new%>%dplyr::mutate(!!x :=apply(temp,1,prod))
 #  }
 # data_to_impute=merge(data_to_impute, new, by=ID, all.x=T)



  write.csv(data_to_impute, paste0(home_dir, "imputations/",  exposure, "-", outcome,"_data_to_impute.csv"))
  cat("See the 'imputations' folder for a csv file of the data to be imputed","\n")

  return(data_to_impute)
}
