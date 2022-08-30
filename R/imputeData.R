#Function to create 5 imputed datasets from original dataset
#' Title
#'
#' @param data_to_impute
#'
#' @return
#' @export
#' @importFrom amelia amelia

#' @examples
imputeData <- function(home_dir, ID, data_to_impute, continuous_variables){

  data_to_impute=as.data.frame(data_to_impute)
  to_remove=c(ID, "WAVE")

  #finds ordinal variables --all others assumed continuou
  ordinal_vars=colnames(data_to_impute)[!colnames(data_to_impute) %in% c(to_remove, continuous_variables)]

  #creates 5 imputed datasets --add more detail here
  imputed=amelia::amelia(x = data_to_impute, m = 5, idvars = ID, #m=5 for final
                 ts = "WAVE", cs = NULL, priors = NULL, lags = NULL, intercs = FALSE, leads = NULL, splinetime = NULL,
                 logs = NULL, sqrts = NULL, lgstc = NULL, noms = NULL,
                 ords=ordinal_vars,
                 empri=0.01*nrow(ALLdata2),
                 # autopri=1,
                 bounds=NULL,  max.resample = 100, #100 seems to be the default
                 tolerance = 1e-04
  )

  #this contains each of the 5 imputed datasets
  imputed_datasets<- imputed$imputations

  #save out imputed datasets
  write.csv(imputed_datasets$imp1, file=paste0(home_dir,"Imputations/imp1.csv"))
  write.csv(imputed_datasets$imp2, file=paste0(home_dir,"Imputations/imp2.csv"))
  write.csv(imputed_datasets$imp3, file=paste0(home_dir,"Imputations/imp3.csv"))
  write.csv(imputed_datasets$imp4, file=paste0(home_dir,"Imputations/imp4.csv"))
  write.csv(imputed_datasets$imp5, file=paste0(home_dir,"Imputations/imp5.csv"))

  return(imputed_datasets)

}
