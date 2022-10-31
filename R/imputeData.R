#' Imputes dataset so there is no missing at each time point
#'
#' Creates m imputed datasets from original dataset using the Aeelia package. See Amelia documentation for more detail: https://cran.r-project.org/web/packages/Amelia/Amelia.pdf

#' @param object msm object that contains all relevant user inputs
#' @param data_to_impute output from dataToImpute
#' @param max.resample Amelia input for how many times amelia should redraw imputed values when trying to meet logical constraints of obounds
#' @param cs Amelia input for col num or var indicating cross section variable
#' @param priors Amelia input col matrix indicateing priors for indiv missing obs or variable-wide missing
#' @param lags Amelia input for cols that should have lags included in imputation model
#' @param intercs Amelia input logical var indicating if time effects of polytime should vary across cross-section
#' @param leads Amelia input for cols that should have leads included in imputation model
#' @param splinetime Amelia input integer of 0+ to control cubic smoothing splines of time
#' @param logs Amelia input for cols that require log-linear transformation
#' @param sqrts Amelia input for cols that require sqrt transformation
#' @param lgstc Amelia input for cols that require lgstc transformation
#' @param noms Amelia input for cols that are nominal
#' @param bounds Amelia input for 3-col matrix for logical bounds on imputations
#' @return imputed_datasets imputation results
#' @export
#' @importFrom Amelia amelia
#' @seealso [dataToImput()], [Amelia::amelia()]
#' @examples imputeData(object, data_to_impute, read_imps_from_file="no", max.resample = 100, cs=NULL, priors=NULL, lags=NULL, intercs=FALSE, leads=NULL, splinetime=NULL, logs=NULL, sqrts=NULL, lgstc=NULL, noms=NULL, bounds=NULL)
#'
imputeData <- function(object, data_to_impute, read_imps_from_file="no", max.resample = 100, cs=NULL, priors=NULL, lags=NULL, intercs=FALSE, leads=NULL, splinetime=NULL, logs=NULL, sqrts=NULL, lgstc=NULL, noms=NULL, bounds=NULL){

  home_dir=object$home_dir
  ID=object$ID
  continuous_variables=object$continuous_variables
  m=object$m

  if (read_imps_from_file=="yes"){
    imputed_datasets=list()
    for (x in 1:m){
      file_name=(paste("imp", x, '.csv', sep=""))
      name=paste("imp", x, sep="")
      imp=suppressWarnings(as.data.frame(readr::read_csv(paste(paste0(home_dir, "imputations/"), file_name, sep=""))))
      imputed_datasets[[paste0("imp", x)]]<-imp
    }
    return(imputed_datasets)

  }else{


    data_to_impute=as.data.frame(data_to_impute)
    to_remove=c(ID, "WAVE")

    #finds ordinal variables --all others assumed continuou
    ordinal_vars=colnames(data_to_impute)[!colnames(data_to_impute) %in% c(to_remove, continuous_variables)]

    #creates 5 imputed datasets --add more detail here
    imputed=Amelia::amelia(x = data_to_impute, m = m, idvars = ID, #m=5 for final
                           ts = "WAVE", cs =cs, priors = priors, lags = lags, intercs = intercs, leads = leads, splinetime = splinetime,
                           logs = logs, sqrts = sqrts, lgstc = lgstc, noms = noms,
                           ords=ordinal_vars,
                           empri=0.01*nrow(data_to_impute),
                           bounds=bounds,  max.resample = max.resample, #100 seems to be the default
                           tolerance = 1e-04
    )

    #this contains each of the imputed datasets
    imputed_datasets<- imputed$imputations

    #save out imputed datasets
    for (k in 1:m){
      write.csv(imputed_datasets[[paste0("imp", k)]], file=paste0(home_dir,"imputations/imp", k,".csv"))
    }
    cat("See the 'imputations' folder for a csv file of each imputed dataset","\n")

    return(imputed_datasets)
  }

}
