#'Formats dataset and creates required directory
#'
#'Creates the required directories, assesses the data structure and makes it uniform for future functions and returns data
#'This function requires a clean dataset in long format that contain columns for ID, time, exposure, outcome, and potential covariate confounds
#'
#' @param object msm object that contains all relevant user inputs
#' @return data formatted dataset
#' @export
#' @importFrom readr read_csv
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @examples formatDataStruct(object, factor_covariates)
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
  outcome_time_pt=object$outcome_time_pt


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

  if(dir.exists(paste0(home_dir, "combined weights/"))==F){dir.create(paste0(home_dir, "combined weights/"))}
  if(dir.exists(paste0(home_dir, "combined weights/values/"))==F){dir.create(paste0(home_dir, "combined weights/values/"))}
  if(dir.exists(paste0(home_dir, "combined weights/histograms/"))==F){dir.create(paste0(home_dir, "combined weights/histograms/"))}

  if(dir.exists(paste0(home_dir, "final weights/"))==F){dir.create(paste0(home_dir, "final weights/"))}
  if(dir.exists(paste0(home_dir, "final weights/values/"))==F){dir.create(paste0(home_dir, "final weights/values/"))}
  if(dir.exists(paste0(home_dir, "final weights/histgrams/"))==F){dir.create(paste0(home_dir, "final weights/histograms/"))}

  if(dir.exists(paste0(home_dir, "forms/"))==F){dir.create(paste0(home_dir, "forms/"))}

  if(dir.exists(paste0(home_dir, "pre-balance/"))==F){dir.create(paste0(home_dir, "pre-balance/"))}
  if(dir.exists(paste0(home_dir, "pre-balance/plots/"))==F){dir.create(paste0(home_dir, "pre-balance/plots/"))}


  if(dir.exists(paste0(home_dir, "balance/"))==F){dir.create(paste0(home_dir, "balance/"))}
  if(dir.exists(paste0(home_dir, "balance/plots/"))==F){dir.create(paste0(home_dir, "balance/plots/"))}
  if(dir.exists(paste0(home_dir, "balance/comparison values/"))==F){dir.create(paste0(home_dir, "balance/comparison values/"))}
  if(dir.exists(paste0(home_dir, "balance/unbalanced covariates/"))==F){dir.create(paste0(home_dir, "balance/unbalanced covariates/"))}
  if(dir.exists(paste0(home_dir, "balance/unbalanced covariates/"))==F){dir.create(paste0(home_dir, "balance/unbalanced covariates/"))}

  if(dir.exists(paste0(home_dir, "balance/potential confounds/"))==F){dir.create(paste0(home_dir, "balance/potential confounds"))}

  if(dir.exists(paste0(home_dir, "msms/"))==F){dir.create(paste0(home_dir, "msms/"))}
  if(dir.exists(paste0(home_dir, "msms/original/"))==F){dir.create(paste0(home_dir, "msms/original/"))}
  if(dir.exists(paste0(home_dir, "msms/sensitivity checks/"))==F){dir.create(paste0(home_dir, "msms/sensitivity checks/"))}
  if(dir.exists(paste0(home_dir, "msms/linear hypothesis testing/"))==F){dir.create(paste0(home_dir, "msms/linear hypothesis testing/"))}
  if(dir.exists(paste0(home_dir, "msms/linear hypothesis testing/original/"))==F){dir.create(paste0(home_dir, "msms/linear hypothesis testing/original/"))}
  if(dir.exists(paste0(home_dir, "msms/linear hypothesis testing/sensitivity checks/"))==F){dir.create(paste0(home_dir, "msms/linear hypothesis testing/sensitivity checks/"))}

  if(dir.exists(paste0(home_dir, "results figures/"))==F){dir.create(paste0(home_dir, "results figures/"))}
  if(dir.exists(paste0(home_dir, "results figures/original/"))==F){dir.create(paste0(home_dir, "results figures/original/"))}
  if(dir.exists(paste0(home_dir, "results figures/sensitivity checks/"))==F){dir.create(paste0(home_dir, "results figures/sensitivity checks/"))}

  if(dir.exists(paste0(home_dir, "for Mplus/"))==F){dir.create(paste0(home_dir, "for Mplus/"))}


  cat("\n")
  #reading and formatting LONG dataset
  data=suppressWarnings(as.data.frame(readr::read_csv(data_path, show_col_types=FALSE)))

  colnames(data)[colnames(data)==time_var] <- "WAVE"
  data[data == missing] <- NA #makes NA the missingness indicator



  #exposure summary
  exp=as.data.frame(data[,c("WAVE", colnames(data)[colnames(data) %in% exposure])])
  exp=exp%>%dplyr::filter(WAVE %in% exposure_time_pts)%>%dplyr::group_by(WAVE)%>%
    dplyr::summarise_all(list(mean=mean, sd=sd, min=min, max=max), na.rm=TRUE)
  exp=exp[,order(colnames(exp), decreasing=TRUE)]
  cat(knitr::kable(exp, caption="Summary of Exposure Information", format='pipe'),  sep="\n")

  invisible(knitr::kable(exp, caption="Summary of Exposure Information", format='html')%>%
        kableExtra::kable_styling()%>%
        kableExtra::save_kable(file=paste0(home_dir, "exposure_info.html")))
  cat("Exposure descriptive statistics have now been saved in the home directory", "\n")
  cat("\n")

  #outcome summary
  exp=as.data.frame(data[,c("WAVE", colnames(data)[colnames(data) %in% outcome])])
  exp=suppressWarnings(exp%>%dplyr::filter(WAVE %in% outcome_time_pt)%>%dplyr::group_by(WAVE)%>%dplyr::summarise_all(list(mean=mean, sd=sd, min=min, max=max), na.rm=TRUE))
  exp=exp[,order(colnames(exp), decreasing=TRUE)]

  cat(knitr::kable(as.data.frame(exp), caption="Summary of Outcome Information", format='pipe'), sep="\n")
  knitr::kable(exp, caption="Summary of Outcome Information", format='html')%>%
    kableExtra::kable_styling()%>%
    kableExtra::save_kable(file=paste0(home_dir, "outcome_info.html"))
  cat("Outcome descriptive statistics have now been saved in the home directory", "\n")


  if (sum(factor_covariates %in% colnames(data))<length(factor_covariates)){
    stop('Please provide factor covariates that correspond to columns in your data when creating the msm object')
  }


  data[,factor_covariates] <- lapply(data[,factor_covariates] , factor)

  data[,id]=factor(data[,id])
  numeric_vars=colnames(data)[!colnames(data) %in% c(factor_covariates, id)]

  data[,numeric_vars] <- lapply(data[,numeric_vars] , as.numeric)


  return(data)
}



