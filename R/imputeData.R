#' Imputes dataset so there is no missing at each time point
#'
#' Creates m imputed datasets from original dataset using the Aeelia package. See Amelia documentation for more detail: https://cran.r-project.org/web/packages/Amelia/Amelia.pdf

#' @param object msm object that contains all relevant user inputs
#' @param data_to_impute output from dataToImpute
#' @return imputed_datasets imputation results
#' @export
#' @importFrom Amelia amelia
#' @seealso [dataToImpute()], [mice::mice()]
#' @examples imputeData(object, data_to_impute, read_imps_from_file="no", max.resample = 100, cs=NULL, priors=NULL, lags=NULL, intercs=FALSE, leads=NULL, splinetime=NULL, logs=NULL, sqrts=NULL, lgstc=NULL, noms=NULL, bounds=NULL)
#'
imputeData <- function(object, data_to_impute, read_imps_from_file="no"){

  home_dir=object$home_dir
  ID=object$ID
  continuous_variables=object$continuous_variables
  m=object$m
  imp_method=object$imp_method
  time_varying_covariates=object$time_varying_variables
  factor_covariates=object$factor_covariates
  exposure=object$exposure
  outcome=object$outcome


  if (read_imps_from_file=="yes"){

    imputed_datasets=list()
    imp=readRDS(paste0(home_dir,"imputations/",  exposure, "-", outcome,"_all_imp.rds"))
    imputed_datasets<-imp
    cat(paste0("Reading in ", imputed_datasets$m, " imputations from local folder"))
    cat("\n")
    return(imputed_datasets)

  }else{

    library(mice)
    set.seed(123)

    cat(paste0("Creating ", m, " imputed datasets using the ", imp_method, " imputation method in mice. This may take some time to run."))
    cat("\n")

    ## Configure parallelization (code from Kazuki Yoshida; https://rpubs.com/kaz_yos/mice-exclude)
    ## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
    library(doParallel)
    ## Reproducible parallelization
    library(doRNG)
    ## Detect core count
    nCores <- min(parallel::detectCores(), 8)
    # nCores<-1
    ## Used by parallel::mclapply() as default
    options(mc.cores = nCores)
    ## Used by doParallel as default
    options(cores = nCores)
    ## Register doParallel as the parallel backend with foreach
    ## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
    doParallel::registerDoParallel(cores = nCores)
    ## Report multicore use
    cat("### Using", foreach::getDoParWorkers(), "cores\n")
    cat("### Using", foreach::getDoParName(), "as backend\n")


    data_to_impute=as.data.frame(data_to_impute)
    data_to_impute[,ID]=as.factor(data_to_impute[,ID])
    # data_to_impute$S_ID=as.factor(data_to_impute$S_ID)

    #adding in missing longitudinal time points in long format
    data_to_impute_full=data_to_impute%>%
      tidyr::complete(.data[[ID]], WAVE)

    #fills in all person-level variables
    time_invar=colnames(data_to_impute)[!colnames(data_to_impute) %in% time_varying_covariates]
    time_invar=time_invar[!time_invar %in% c(ID, "WAVE")]
    person_info=data_to_impute[!duplicated(data_to_impute$S_ID),]%>%dplyr::select(c(all_of(ID), all_of(time_invar)))

    for (x in 1:nrow(person_info)){
      id=as.character(person_info[x, ID])
      data_to_impute_full[data_to_impute_full[colnames(data_to_impute)==ID]==id, colnames(data_to_impute_full) %in% time_invar]=person_info[x,time_invar]
    }

    #making variable types
    data_to_impute_full[,colnames(data_to_impute_full) %in% factor_covariates]=lapply(data_to_impute_full[,colnames(data_to_impute_full) %in% factor_covariates], factor)


    ## Set seed for reproducibility
    # set.seed(123)

    # M=m
    # Parallelized execution
    miceout <- foreach(i = seq_len(m), .combine = mice::ibind) %dorng% {
    cat("### Started iteration", i, "\n")
    miceout <- mice::mice(data_to_impute_full, m=1, method=imp_method, maxit = 5,
                          print = F)
    # miceout <- mice(data = df_before, m = 1, print = TRUE,
    #                 predictorMatrix = predictorMatrix, method = dryMice$method,
    #                 MaxNWts = 2000)
    cat("### Completed iteration", i, "\n")
    ## Make sure to return the output
      miceout
    }
    imputed=miceout

    # library(miceadds)
    saveRDS(imputed, paste0(home_dir,"imputations/",  exposure, "-", outcome,"_all_imp.rds"))

    #print warnings
    cat("USER ALERT: Please view any logged events from the imputation below:", "\n")
    # print(imputed$loggedEvents)
    cat(knitr::kable(imputed$loggedEvents, caption="Logged Events from mice", format='pipe'),  sep="\n")
    cat("\n")

    #this contains each of the imputed datasets
    imputed_datasets<- imputed

    #save out imputed datasets
    for (k in 1:m){
      write.csv(mice::complete(imputed,k), file=paste0(home_dir,"imputations/",  exposure, "-", outcome,"_imp", k,".csv"))

    }
    cat("See the 'imputations' folder for a .csv file of each imputed dataset","\n")

    return(imputed_datasets)
  }

}
