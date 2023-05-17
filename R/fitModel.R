
#' Fit weighted model
#' This code fits a weighted marginal structural model to examine the effects of different exposure histories on outcome
#' @param object msm object that contains all relevant user inputs
#' @param data_for_model_with_weights_cutoff dataset with truncated weights see truncateWeights
#' @param unbalanced_covariates_for_models unbalanced covariates see assessBalance
#' @importFrom survey svydesign
#' @importFrom survey svyglm
#' @return marginal structural models
#' @export
#' @seealso [truncateWeights()], [asesssBalance()]
#' @examples fitModel(object, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models)

fitModel <- function(object, data_for_model_with_weights_cutoff, balance_stats_final, model="m3"){

  home_dir=object$home_dir
  ID=object$ID
  exposure=object$exposure
  exposure_epochs=object$exposure_epochs
  outcome=object$outcome
  outcome_time_pt=object$outcome_time_pt
  factor_covariates=object$factor_covariates
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))
  m=object$m
  balance_thresh=object$balance_thresh
  time_varying_covariates=object$time_varying_variables
  time_pts=object$time_pts


  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(balance_stats_final, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=balance_stats_final[[1]]$exp_time,
                               covar_time=balance_stats_final[[1]]$covar_time,
                               covariate=balance_stats_final[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0))%>%
    dplyr::filter(balanced_avg==0)%>%
    dplyr::select(covariate)
  #getting only time invariant or t1 imbalanced covariates
  unbalanced_covars$keep=ifelse(grepl(paste0(".", time_pts[1]), unlist(unbalanced_covars$covariate)), 1, 0)
  unbalanced_covars$keep=ifelse(unbalanced_covars$keep==0 &
                                  !gsub(" ", "",sapply(strsplit(unlist(unbalanced_covars$covariate), "\\."), "[",1)) %in% time_varying_covariates, 1, unbalanced_covars$keep)
  baseline_covars=unlist(unbalanced_covars%>%dplyr::filter(keep==1)%>%dplyr::select(covariate))
  nonb_covars=unlist(unbalanced_covars%>%dplyr::filter(keep==0)%>%dplyr::select(covariate))

  #listing any imbalanced covariates
  #renames factors (that were appended w/ level)
  baseline_covars[sapply(strsplit(sapply(strsplit(baseline_covars, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(baseline_covars, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(baseline_covars, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]
  baseline_covars=baseline_covars[!grepl(exposure, baseline_covars)] #excludes exposure
  covariate_list= paste(as.character(baseline_covars), sep="", collapse=" + ")


  if (!model %in% c("m0", "m1", "m2", "m3")){
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")')
  }

  if (length(outcome_time_pt)>1){
    stop('This function is designed only for single time point outcome')}

  if(model=="m1" | model =="m3"){
    #if there are no imbalanced covariates
    if(covariate_list[1]==""){
      stop("You have selected a covariate model but there are no imbalanced baseline covariates. Please choose another model.")
    }else{ #there are imbalanced covariates
      cat("USER ALERT: The following imbalanced baseline covariates are included in the model: ", "\n")
      cat(covariate_list, "\n")
      cat("\n")


      cat("\n")
      cat("USER ALERT: The following imbalanced covariates will NOT be included in the model: ", "\n")

      cat(nonb_covars, "\n")
      cat("\n")

    }
  }

  #lists out exposure-epoch combos
  exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")

  #getting interactions
  interactions<-paste(lapply(2:length(exp_epochs), function(z){
    paste(apply(combn(exp_epochs,z), 2, paste, sep="", collapse=":"), sep="", collapse=" + ")
  }), collapse= " + ")

  fits= lapply(1:length(data_for_model_with_weights_cutoff), function(y){ #cutoffs
    d=data_for_model_with_weights_cutoff[[y]]
    cutoff=all_cutoffs[y]
    cutoff_label=ifelse(cutoff==weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))

    lapply(d, function(x){ #imputed datasets
      data=x
      data$weights=NULL
      data$weights<-data[,colnames(data)[grepl("weight", colnames(data))]]

      s=survey::svydesign(id=~1, #
                          data=data, #adds list of imputation data
                          weights=~weights)

      f0=paste(paste0(outcome, ".", outcome_time_pt), "~", paste0(exp_epochs, sep="", collapse=" + "))
      m0=survey::svyglm(as.formula(f0), design=s) #list of model fitted to all imputed datasets

      #baseline model (main effects) only
      if(model=="m0"){
        return(m0) #save model
      }else{
        #baseline + sig covar model OR baseline + sig covar + int model
        if(model=="m1" | model =="m3"){
          f1.temp=paste(f0, "+", covariate_list) #baseline + covariate model
          f1=f1.temp
          m1=survey::svyglm(as.formula(f1), design=s)

          #baseline + imbalanced covars
          if(model=="m1"){
            return(m1)
          }
        }

        #baseline + interaction OR baseline + covars + interactions
        if (model=="m2" | model=="m3"){
          f2=paste(f0, "+", paste(interactions, sep="", collapse=" + "))

          #baseline + interactions
          if (model=="m2"){
            m2=survey::svyglm(as.formula(f2), design=s)
            #baseline + interactions
            return(m2)
          }

          #baseline + covars + interactions
          if (model=="m3"){
            f3=paste(f1, "+", paste(interactions, sep="", collapse=" + "))
            m3=survey::svyglm(as.formula(f3), design=s)
            #baseline + covars+ interactions
            return(m3)
          }
        }
      }
    })
  })
  names(fits)=all_cutoffs

  cat(paste0("The marginal model, ", model, ", run for each imputed dataset using the user-specified weights truncation percentile value of ",
             paste(weights_percentile_cutoff), " is summarized here:"), "\n")
  # browser()

  print(lapply(fits[paste0(weights_percentile_cutoff)],function(x){lapply(x, summary)}))

  # browser()

  #print table of model evidence comparing imputations per cutoff level
  lapply(1:length(fits), function(y){
    i=fits[[y]]
    suppressWarnings(jtools::export_summs(i, to.file="docx", statistics= c(N="nobs", AIC="AIC", R2="r.squared"),
                                          model.names= c(paste0("Imp.", 1:length(i))),
                                          file.name =paste0(home_dir, "msms/", names(fits)[y],"_", model, "_table_mod_ev.docx", sep="")))
  })

  cat("Tables of model evidence for all truncation percentile values have now been saved in the 'msm' folder.")

  m0=NULL
  m1=NULL
  m2=NULL
  m3=NULL


  cat("\n")

  # file_label=ifelse(cutoff==weights_percentile_cutoff, "original", "sensitivity checks")
  saveRDS(fits, file = paste(paste0(home_dir, "msms/exposure-outcome_", model, "_model.rds", sep="")))


  cat("\n")
  # cat("A new data file has been saved as a .csv file in the in the 'msms' folder","\n")

  return(fits)

}
