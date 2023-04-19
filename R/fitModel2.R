
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
  #makes wide for modeling
  unbalanced_covars=sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1)
  #listing any imbalanced covariates
  baseline_covars=unlist(unbalanced_covars)[!gsub(" ", "",sapply(strsplit(unlist(unbalanced_covars), "\\."), "[",1)) %in% time_varying_covariates]
  #renames factors (that were appended w/ level)
  baseline_covars[sapply(strsplit(sapply(strsplit(baseline_covars, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(baseline_covars, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(baseline_covars, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]
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
      cat("The following imbalanced baseline covariates are included in the covariate model: ", "\n")
      cat(covariate_list)
    }
  }

  library(MatchThem)


  #lists out exposure-epoch combos
  exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")
  #getting interactions
  interactions=paste(apply(combn(exp_epochs,2), 2, paste, sep="", collapse=":"), sep="", collapse=" + ")


  fits= lapply(1:length(dat_w_t), function(y){ #cutoffs
    d=dat_w_t[[y]]
    cutoff=all_cutoffs[y]
    cutoff_label=ifelse(cutoff==weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))

    lapply(d, function(x){ #imputed datasets
      d=x
      s=survey::svydesign(~ID,#~1, #
                          data=d, #adds list of imputation data
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

          #baseline + interactions
          if (model=="m2"){
            f2=paste(f0, "+", paste(interactions, sep="", collapse=" + "))
            m2=survey::svyglm(as.formula(f2), design=s)
            #baseline + interactions
            return(m2)
          }

          #baseline + covars + interactions
          if (model=="m3"){
            f3=paste(f2, "+", paste(interactions, sep="", collapse=" + "))
            m3=survey::svyglm(as.formula(f3), design=s)
            #baseline + covars+ interactions
            return(m3)
          }
        }
      }

    })
  })
  names(fits)=all_cutoffs

  cat(paste0("The marginal model, ", model, ", run for each imputed dataset using the weights truncation cutoffs ", paste(all_cutoffs, collapse=","), " is summarized here:"), "\n")
  # print(summary(models[[model]]))
  print(lapply(fits,function(x){lapply(x, summary)}))

  # all_models[[paste0("fit", k, "_", exposure, "-", outcome, "_cutoff_", cutoff)]]<-models

  # models=NULL
  m0=NULL
  m1=NULL
  m2=NULL
  m3=NULL


  cat("\n")

  file_label=ifelse(cutoff==weights_percentile_cutoff, "original", "sensitivity checks")
  saveRDS(fits, file = paste(paste0(home_dir, "msms/", file_label, "/all_exposure-outcome_models_weights.rds", sep="")))




  # })

  # browser()
  # names(all_models)<-paste0(exposure, "-", outcome, "_cutoff_", all_cutoffs)

  # cat(paste0("All models with ", file_label, " cutoff weights have been saved as an .rds object in the 'msms/", file_label, "/' folder"),"\n")


  # write.csv(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.csv"))
  # saveRDS(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.rds"))

  cat("USER ALERT: Sensitivity check models using two alternate weight truncation values have been saved in the 'msms/sensitivity checks/' folder")
  cat("\n")
  # cat("A new data file has been saved as a .csv file in the in the 'msms' folder","\n")


  return(fits)

}
