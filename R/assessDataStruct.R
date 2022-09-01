
#This function creates the required directories, assesses the data structure including time points, identifiers, and potential covariates.
#' assessDataStruct
#'
#' @param data_path path to cleaned dataset
#' @param home_dir path to home directory for the project
#' @param factor_covariates list of covariates that are factors


#' @return data formatted dataset

#' @export
#' @importFrom readr read_csv
#' @examples
#'
#This function requires a clean dataset in long format that contain columns for ID, time, exposure, outcome, and potential covariate confounds
formatDataStruct <-function(data_path, home_dir, missing, factor_covariates) {
# require(readr)
  options(readr.num_columns = 0)

  #creating all necessary directories within the home directory
  if(dir.exists(paste0(home_dir, "imputations/"))==F){dir.create(paste0(home_dir, "imputations/"))}
  if(dir.exists(paste0(home_dir, "weight fits/"))==F){dir.create(paste0(home_dir, "weight fits/"))}
  if(dir.exists(paste0(home_dir, "weights/"))==F){dir.create(paste0(home_dir, "weights/"))}
  if(dir.exists(paste0(home_dir, "combined weights/"))==F){dir.create(paste0(home_dir, "combined weights/"))}
  if(dir.exists(paste0(home_dir, "final weights/"))==F){dir.create(paste0(home_dir, "final weights/"))}
  if(dir.exists(paste0(home_dir, "plots/"))==F){dir.create(paste0(home_dir, "plots/"))}
  if(dir.exists(paste0(home_dir, "forms/"))==F){dir.create(paste0(home_dir, "forms/"))}
  if(dir.exists(paste0(home_dir, "balance/"))==F){dir.create(paste0(home_dir, "balance/"))}

  #reading and formatting LONG dataset
  data=as.data.frame(readr::read_csv(data_path), col_types=cols(), show_col_types=FALSE)
  colnames(data)[colnames(data)==time_var] <- "WAVE"
  data[data == missing] <- NA
  data[,factor_covariates] <- lapply(data[,factor_covariates] , factor)

  return(data)
}


#' makeTimePtDatasets
#' @param data output from formatDataStruct
#' @param ID person-level identifier in your dataset
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @param time_var variable in your long dataset that designates developmental time
#' @param missing missing data marker in your dataset
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param outcomes list of variables that represent your outcomes of interest
#'
#' @return time_pt_datasets datasets split by time point
#' @export
#' @importFrom dplyr filter
#' @examples
makeTimePtDatasets <-function(data, ID, time_pts, time_var, missing, exposures, outcomes){
  require(dplyr)
  #creating dataset for each time point and store in a list
  data_list={}
  time_pt_datasets=list()
  for (t in 1:length(time_pts)){
    # assign(paste0("Data_", time_pts[t]), data%>% dplyr::filter(WAVE==time_pts[t]))
    time_pt_datasets[[paste0("Data_", time_pts[t])]] <-data%>%
      dplyr::filter(WAVE==time_pts[t])
    data_list=c(data_list, paste0("Data_", time_pts[t]))
  }

  return(time_pt_datasets)

}
