
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
assessModel <-function(object, all_models){

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

    best_fit=names(which.min(unlist(aics)))
    best_fit_model=all_models[[exp_out]][names(all_models[[exp_out]])==best_fit]



    best_models[[exp_out]]<-best_fit_model[[1]]


    cat(paste0("The best-fitting model for ", exp_out, " is ", best_fit),"\n")
    print(summary(best_fit_model[[1]]))

    # browser()
    folder_label=ifelse(as.numeric(sapply(strsplit(exp_out, "_"), "[", 3))==weights_percentile_cutoff, "original/", "sensitivity checks/")

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
