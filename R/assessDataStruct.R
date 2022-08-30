
#This function assess the data structure including time points, identifiers, and potential covariates.
#' Title
#'
#' @param data
#'
#' @return
#' @export
#' @importFrom readr read_csv
#' @importFrom dplyr

#' @examples
#'

#This function requires a clean dataset in long format that contain columns for ID, time, exposure, outcome, and potential covariate confounds
assessDataStruct <-function(data_path, home_dir, ID, time_pts, time_var, missing, exposures, outcomes){

  #creating all necessary directories within the home directory
  if(dir.exists(paste0(home_dir, "imputations/"))==F){dir.create(paste0(home_dir, "imputations/"))}
  if(dir.exists(paste0(home_dir, "weight fits/"))==F){dir.create(paste0(home_dir, "weight fits/"))}
  if(dir.exists(paste0(home_dir, "weights/"))==F){dir.create(paste0(home_dir, "weights/"))}
  if(dir.exists(paste0(home_dir, "combined weights/"))==F){dir.create(paste0(home_dir, "combined weights/"))}
  if(dir.exists(paste0(home_dir, "final weights/"))==F){dir.create(paste0(home_dir, "final weights/"))}
  if(dir.exists(paste0(home_dir, "plots/"))==F){dir.create(paste0(home_dir, "plots/"))}
  if(dir.exists(paste0(home_dir, "forms/"))==F){dir.create(paste0(home_dir, "forms/"))}
  if(dir.exists(paste0(home_dir, "balance/"))==F){dir.create(paste0(home_dir, "balance/"))}


  #reading and formatting dataset
  data=as.data.frame(readr::read_csv(data_path))
  colnames(data)[colnames(data)==time_var] <- "WAVE"
  data[data == missing] <- NA
  data[,factor_covariates] <- lapply(data[,factor_covariates] , factor)


  #identifying potential covariates (time-varying and time invariant)
  potential_covariates=colnames(data)[colnames(data) %in% (c(ID, "WAVE", exposures, outcomes))==FALSE]

  #creating dataset for each time point
  data_list={}
  for (t in 1:length(time_pts)){
    assign(paste0("Data_", time_pts[t]), data%>% dplyr::filter(WAVE==time_pts[t]))
    data_list=c(data_list, paste0("Data_", time_pts[t]))
  }

  return(data, potential_covariates)

}
