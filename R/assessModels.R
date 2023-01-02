
#' Assess marginal structural models
#' Code to assess the various models and decide which fits the best
#'
#' @param object msm object that contains all relevant user inputs
#' @param all_models output from fitModel
#' @importFrom sjPlot tab_model
#' @importFrom jtools export_summs
#' @return best_models list of best-fitting models for each exposure=outcome pairing
#' @seealso [fitModel()]
#' @examples assessModel(object, all_models)
assessModel <-function(object, all_models, best_model=NULL){

  require(huxtable) #not sure if this is the best way to load these? not associated w/ particular function but needed for jtools
  require(officer)
  require(flextable)

  home_dir=object$home_dir
  weights_percentile_cutoff=object$weights_percentile_cutoff

  best_models=list()
  cat("USER ALERT: please inspect the following best-fitting models for each exposure-outcome pair for original and sensitivity check weight cutoff values:", "\n")
  cat("\n")

  for (x in 1:length(all_models)){

    exp_out=names(all_models)[x]

    model_list=c(names(all_models[[exp_out]]))

    aics=lapply(all_models[[exp_out]], function(x) x$aic)

    aics=aics[!names(aics) %in% c("m1", "m3")] #we only want to compare baseline (m0) to model with sig covars (m2) to model with sig interaction (m4)

    if (is.null(best_model)){ #if user did not specify which model they want
    best_fit=names(which.min(unlist(aics)))
    }else{ #otherwise use user-specified bests model
      best_fit=best_model
    }

    best_fit_model=all_models[[exp_out]][names(all_models[[exp_out]])==best_fit]

    #looks to see if best model is in the model list (may not be as not all exposure-outcome pairs had additional covariates to try for covariate models)
    # t=try(best_models[[exp_out]]<-best_fit_model[[1]], silent=T)
    # if("try-error" %in% class(t)){
    if(length(best_fit_model)==0){
      cat(paste0("USER ALERT: the user-specified model, ", best_fit, ", does not exist for ", exp_out, " so the existing model with the lowest AIC will be used."), "\n")
      best_fit=names(which.min(unlist(aics))) #if model does not exist, choose lowest AIC model
      best_fit_model=all_models[[exp_out]][names(all_models[[exp_out]])==best_fit]
    }

    best_models[[exp_out]]<-best_fit_model[[1]]

    cat(paste0("The identified best-fitting model for ", exp_out, " is ", best_fit, ":"),"\n")
    print(summary(best_fit_model[[1]]))

    # browser()
    folder_label=ifelse(as.numeric(sapply(strsplit(exp_out, "_"), "[", 3))==weights_percentile_cutoff, "original/", "sensitivity checks/")

    #print table of model evidence
    suppressWarnings(jtools::export_summs(all_models[[exp_out]], to.file="docx", statistics= c(N="nobs", AIC="AIC", R2="r.squared"),
                                          file.name =paste0(home_dir, "msms/",folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_table_mod_ev.docx", sep="")))


    cat(paste0("See the 'msms/original/' folder for tables of model evidence for ", exp_out, " with original weight cutoffs"), "\n")
    cat("\n")
    cat(paste0("See the 'msms/sensitivity checks/' folder for tables of model evidence for ", exp_out, " with sensitivity check weight cutoffs"), "\n")
    cat("\n")

    cat("\n")


  }
  return(best_models)
}
