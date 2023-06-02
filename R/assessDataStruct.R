#'Formats dataset and creates required directory
#'
#'Creates the required directories, assesses the data structure and makes it uniform for future functions and returns data
#'This function requires a clean dataset in long format that contain columns for ID, time, exposure, outcome, and potential covariate confounds
#'
#' @param object msm object that contains all relevant user inputs
#' @return data formatted dataset
#' @export
#' @importFrom gtools permutations
#' @importFrom readr read_csv
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @examples formatDataStruct(object)
#'
formatDataStruct <-function(object) {

  data_path=object$data_path
  home_dir=object$home_dir
  missing=object$missing
  time_var=object$time_var
  time_varying_covariates=object$time_varying_variables
  factor_covariates=object$factor_covariates
  id=object$ID
  exposure=object$exposure
  outcome=object$outcome
  exposure_time_pts=object$exposure_time_pts
  exposure_epochs=object$exposure_epochs
  outcome_time_pt=object$outcome_time_pt
  time_varying_covariates=object$time_varying_variables
  time_pts=object$time_pts
  time_var_exclude=object$time_var_exclude
  hi_cutoff=object$hi_cutoff
  lo_cutoff=object$lo_cutoff

  cat("Formatting data")

  options(readr.num_columns = 0)

  #error checking
  if (!file.exists(data_path)){
    stop('Please provide a valid directory for your data in data_path when creating the msm object')
  }
  if (!dir.exists(home_dir)){
    stop('Please provide a valid home directory in home_dir when creating the msm object')
  }


  #creating all necessary directories within the home directory
  if(dir.exists(paste0(home_dir, "imputations/"))==F){dir.create(paste0(home_dir, "imputations/"))}
  if(dir.exists(paste0(home_dir, "original weights/"))==F){dir.create(paste0(home_dir, "original weights/"))}
  if(dir.exists(paste0(home_dir, "original weights/values/"))==F){dir.create(paste0(home_dir, "original weights/values/"))}
  if(dir.exists(paste0(home_dir, "original weights/histograms/"))==F){dir.create(paste0(home_dir, "original weights/histograms/"))}

  if(dir.exists(paste0(home_dir, "final weights/"))==F){dir.create(paste0(home_dir, "final weights/"))}
  if(dir.exists(paste0(home_dir, "final weights/values/"))==F){dir.create(paste0(home_dir, "final weights/values/"))}
  if(dir.exists(paste0(home_dir, "final weights/histgrams/"))==F){dir.create(paste0(home_dir, "final weights/histograms/"))}

  if(dir.exists(paste0(home_dir, "forms/"))==F){dir.create(paste0(home_dir, "forms/"))}

  if(dir.exists(paste0(home_dir, "pre balance/"))==F){dir.create(paste0(home_dir, "pre balance/"))}
  if(dir.exists(paste0(home_dir, "pre balance/plots/"))==F){dir.create(paste0(home_dir, "pre balance/plots/"))}

  if(dir.exists(paste0(home_dir, "balance/"))==F){dir.create(paste0(home_dir, "balance/"))}
  if(dir.exists(paste0(home_dir, "balance/plots/"))==F){dir.create(paste0(home_dir, "balance/plots/"))}
  if(dir.exists(paste0(home_dir, "balance/comparison values/"))==F){dir.create(paste0(home_dir, "balance/comparison values/"))}

  if(dir.exists(paste0(home_dir, "msms/"))==F){dir.create(paste0(home_dir, "msms/"))}
  if(dir.exists(paste0(home_dir, "msms/estimated means/"))==F){dir.create(paste0(home_dir, "msms/estimated means/"))}
  if(dir.exists(paste0(home_dir, "msms/contrasts/"))==F){dir.create(paste0(home_dir, "msms/contrasts/"))}

  if(dir.exists(paste0(home_dir, "results figures/"))==F){dir.create(paste0(home_dir, "results figures/"))}



  cat("\n")
  #reading and formatting LONG dataset
  data=suppressWarnings(as.data.frame(readr::read_csv(data_path, show_col_types=FALSE)))

  colnames(data)[colnames(data)==time_var] <- "WAVE" #assigning time variable
  data[data == missing] <- NA #makes NA the missingness indicator



  #exposure summary
  exp=as.data.frame(data[,c("WAVE", colnames(data)[colnames(data) %in% exposure])])
  exp=exp%>%dplyr::filter(WAVE %in% exposure_time_pts)%>%dplyr::group_by(WAVE)%>%
    dplyr::summarise_all(list(mean=mean, sd=sd, min=min, max=max), na.rm=TRUE)
  exp=exp[,order(colnames(exp), decreasing=TRUE)]
  cat(knitr::kable(exp, caption=paste0("Summary of ", exposure, " Exposure Information"), format='pipe'),  sep="\n")
  invisible(knitr::kable(exp, caption=paste0("Summary of ", exposure, " Exposure Information"), format='html')%>%
              kableExtra::kable_styling()%>%
              kableExtra::save_kable(file=paste0(home_dir, exposure, "_exposure_info.html")))
  cat(paste0(exposure, " exposure descriptive statistics have now been saved in the home directory"), "\n")
  cat("\n")


  #exposure history summary
  time_varying_wide=apply(expand.grid(time_varying_covariates, as.character(time_pts)), 1, paste, sep="", collapse=".")
  time_varying_wide=sort(time_varying_wide)
  time_varying_wide=c(id, time_varying_wide)
  time_varying_wide=time_varying_wide[!time_varying_wide %in% time_var_exclude] #rem
  dat_wide=suppressWarnings(stats::reshape(data=data,idvar=id,v.names= time_varying_covariates, timevar="WAVE", times=c(time_pts), direction="wide"))
  new=data.frame(ID=dat_wide[,id])
  colnames(new)<-id
  #averages exposure across time points that constitute the exposure epochs (e.g., infancy = 6 &15)
  for (e in 1:nrow(exposure_epochs)){
    epoch=exposure_epochs[e,1]
    temp=data.frame(row.names=1:nrow(dat_wide))
    new_var=paste0(exposure, "_", epoch)
    #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
    for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
      level=as.numeric(unlist(exposure_epochs[e,2]))[l]
      z=dat_wide[,which(grepl(paste0(exposure, ".", level), names(dat_wide)))]
      temp=cbind(temp, z) }
    new=new%>%dplyr::mutate(!!new_var :=rowMeans(temp, na.rm=T))
    new[,new_var]=as.numeric(new[,new_var])
  }
  #assigning history (e.g., h-h-h) based on user-specified hi/lo cutoffs
  tot_hist=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")
  new$history<- lapply(1:nrow(new), function(x){
    paste(lapply(1:nrow(exposure_epochs), function(y){
      if(is.nan(new[x,y+1])){
        return(NA)}
      if(new[x,y+1]>=as.numeric(quantile(new[,y+1],probs=hi_cutoff, na.rm=T))){
        return("h")}
      if(new[x,y+1]<=as.numeric(quantile(new[,y+1],probs=lo_cutoff, na.rm=T))){
        return("l")}
    }), collapse="-")
  })
  #summarizing n's by history
  his_summ=as.data.frame(as.data.frame(new)%>%dplyr::group_by(history)%>%dplyr::summarize(n=dplyr::n()))
  his_summ=his_summ[!grepl("NULL", unlist(his_summ$history)),]
  his_summ=his_summ[!grepl("NA", unlist(his_summ$history)),]
  cat(paste0("USER ALERT: Out of the total of ", nrow(dat_wide), " individuals in the sample, below is the distribution of the ", sum(his_summ$n), " (",
             round((sum(his_summ$n)/nrow(dat_wide))*100,2), "%) that fall into ", nrow(his_summ), " out of the ", length(tot_hist) ,
             " the total user-defined exposure histories created from ",
             lo_cutoff*100, "th and ", hi_cutoff*100, "th percentile values for low and high levels of exposure ", exposure,
             ", respectively, across ", paste(exposure_epochs$epochs, collapse=", "),
             ". Please inspect the distribution of the sample across the following exposure histories and ensure there is sufficient spread to avoid extrapolation and low precision :"), "\n")
  if(nrow(his_summ)!=length(tot_hist)){
    cat(paste0("USER ALERT: There are no individuals in your sample that fall into ",tot_hist[!tot_hist %in% his_summ$history] ,
               " exposure history/histories. You may wish to consider different high/low cutoffs or choose a different measure to avoid extrapolation."), "\n")
  }
  cat(knitr::kable(as.data.frame(his_summ), caption=paste0("Summary of User-Specified Exposure ", exposure, " Histories Based on Exposure Epochs"), format='pipe', row.names = F ), sep="\n")
  cat("\n")


  #outcome summary
  exp=as.data.frame(data[,c("WAVE", colnames(data)[colnames(data) %in% outcome])])
  exp=suppressWarnings(exp%>%dplyr::filter(WAVE %in% outcome_time_pt)%>%dplyr::group_by(WAVE)%>%dplyr::summarise_all(list(mean=mean, sd=sd, min=min, max=max), na.rm=TRUE))
  exp=exp[,order(colnames(exp), decreasing=TRUE)]
  cat(knitr::kable(as.data.frame(exp), caption=paste0("Summary of Outcome ", outcome, " Information"), format='pipe'), sep="\n")
  knitr::kable(exp, caption=paste0("Summary of Outcome ", outcome, " Information"), format='html')%>%
    kableExtra::kable_styling()%>%
    kableExtra::save_kable(file=paste0(home_dir, outcome, "_outcome_info.html"))
  cat(paste0(outcome, " outcome descriptive statistics have now been saved in the home directory"), "\n")

  if (sum(factor_covariates %in% colnames(data))<length(factor_covariates)){
    stop('Please provide factor covariates that correspond to columns in your data when creating the msm object')
  }

  #formatting factor covariates
  data[,factor_covariates] <- lapply(data[,factor_covariates] , factor)
  data[,id]=factor(data[,id])

  #formatting numeric covariates
  numeric_vars=colnames(data)[!colnames(data) %in% c(factor_covariates, id)]
  data[,numeric_vars] <- lapply(data[,numeric_vars] , as.numeric)


  return(data)
}



