
#' Assess marginal structural models
#' Code to assess the various models and decide which fits the best
#'
#' @param home_dir path to home directory for the project
#' @param all_models output from fitModel
#' @importFrom sjPlot tab_model
#' @importFrom jtools export_summs
#' @return best_models list of best-fitting models for each exposure=outcome pairing
#' @seealso [fitModel()]
#' @examples assessModel(home_dir, all_models)
assessModel <-function(home_dir, all_models){

  best_models=list()

  for (x in 1:length(all_models)){

    exp_out=names(all_models)[x]

    model_list=c(names(all_models[[exp_out]]))

    # survey::regTermTest(all_models[[exp_out]])
    # extractAIC(all_models[[exp_out]])

    aics=lapply(all_models[[exp_out]], function(x) x$aic)

    aics=aics[!names(aics) %in% c("m1", "m3")] #we only want to compare baseline (m0) to model with sig covars (m2) to model with sig interaction (m4)

    best_fit=names(which.min(unlist(aics)))
    best_fit_model=all_models[[exp_out]][names(all_models[[exp_out]])==best_fit]

    # summary(best_fit_model[[1]])

    best_models[[exp_out]]<-best_fit_model[[1]]


    print(paste0("The best-fitting model for ", exp_out, " is ", best_fit))

    jtools::export_summs(all_models[[exp_out]], to.file="docx", file.name =paste0(home_dir, "msms/", sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_table_mod_ev.docx", sep=""))

    #saves out table of model evidence
    # # library(sjPlot)
    # suppressWarnings(sjPlot::tab_model(all_models[[exp_out]], show.p=T, show.aicc=T, show.loglik=T,
    #                   dv.labels = model_list,
    #                   p.style = "numeric_stars",
    #                   file=paste0(home_dir, "msms/", sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_table_mod_ev.html", sep="")
    # ))
    print(paste0("See the 'msms' folder for tables of model evidence for ", exp_out))


  }
  return(best_models)
}
