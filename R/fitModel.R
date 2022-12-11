
#' Fit weighted model
#' This code fits a weighted marginal structural model to examine the effects of different exposure histories on outcomes
#' @param object msm object that contains all relevant user inputs
#' @param data_for_model_with_weights_cutoff dataset with truncated weights see truncateWeights
#' @param unbalanced_covariates_for_models unbalanced covariates see assessBalance
#' @importFrom survey svydesign
#' @importFrom survey svyglm
#' @return marginal structural models
#' @export
#' @seealso [truncateWeights()], [asesssBalance()]
#' @examples fitModel(object, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models)

fitModel <- function(object, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models){

  home_dir=object$home_dir
  ID=object$ID
  exposures=object$exposures
  exposure_epochs=object$exposure_epochs
  outcomes=object$outcomes
  outcome_time_pt=object$outcome_time_pt
  factor_covariates=object$factor_covariates
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))


  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)


  if (length(outcome_time_pt)>1){
    stop('This function is designed only for single time point outcomes')}

  names(data_for_model_with_weights_cutoff) <- gsub(".", "_", names(data_for_model_with_weights_cutoff), fixed=TRUE)

  data_for_model_with_weights_cutoff=as.data.frame(data_for_model_with_weights_cutoff)



  #model dictionary:
  #m0 baseline model with only main effects for each exposure epoch
  #m1 full covariate model adding in all unbalanced covariates (if any)
  #m2 covariate model including only significant covariates (retaining all epoch main effects)
  #m3 full interaction model adding in all possible epoch interactions
  #m4 final model including all epoch main effects and only the significant covariates and interactions


  #cycles through each exposure to fit the marginal structural model relating exposure histories to each outcome
  all_models=list()

  cat("USER ALERT: please inspect the following series of models for each exposure-outcome pair (using weights cutoff at the original user-specified percentile value):", "\n")
  cat("\n")
  cat("\n")

  for (c in 1:length(all_cutoffs)){
    cutoff=all_cutoffs[c]

    for (x in 1:length(exposures)){
      exposure=exposures[x]

      #calculates the mean value for each exposure for each exposure epoch
      for (e in 1:nrow(exposure_epochs)){
        epoch=exposure_epochs[e,1]
        temp=data.frame(row.names=1:nrow(data_for_model_with_weights_cutoff))
        new_var=paste0(exposure, "_", epoch)

        #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
        for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
          level=as.numeric(unlist(exposure_epochs[e,2]))[l]
          z=data_for_model_with_weights_cutoff[,which(grepl(paste0(exposure, "_", level), names(data_for_model_with_weights_cutoff)))]
          temp=cbind(temp, z)
        }
        #adds a new variable of the exposure averaged within epoch
        data_for_model_with_weights_cutoff=data_for_model_with_weights_cutoff%>%dplyr::mutate(!!new_var :=rowMeans(temp, na.rm=T))
      }



      #lists out exposure-epoch combos
      exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")


      #cycles through each outcome
      models=list()


      for (y in 1:length(outcomes)){
        outcome=outcomes[y]

        #creates initial survey design type specifying ID, data, and weights to use
        s=svydesign(ids=data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==ID],
                    data=data_for_model_with_weights_cutoff,
                    weights=data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==
                                                                 paste0(exposure,"-", outcome, "_", paste(unlist(strsplit(as.character(cutoff), "\\.")), sep="", collapse="_"), "_weight_cutoff")])


        #creates form for baseline model including averaged exposure at each epoch as a predictor
        f0=paste(paste0(outcome, "_", outcome_time_pt), "~", paste0(exp_epochs, sep="", collapse=" + "))

        cutoff_label=ifelse(cutoff==weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))

        #fits initial baseline model
        m0=svyglm(noquote(f0), design=s)

        if(cutoff==weights_percentile_cutoff){
          cat(paste0("Baseline model results for effects of ", exposure, " on ", outcome, " using ", cutoff_label),"\n")
          print(summary(m0))
          cat("\n")
        }
        models[["m0"]]<-m0


        #adding in covariates that did not fully balance when creating the weights
        covariate_list= unbalanced_covariates_for_models[[paste0(exposure, "-", outcome)]]
        # grepl(paste(factor_covariates, collapse="|"), paste0(unlist(strsplit(noquote(covariate_list), "\\+"))))

        if (covariate_list[1]==""){
          if(cutoff==weights_percentile_cutoff){
            cat(paste0("There are no unbalanced covariates to include in the model of effects of ", exposure, " on ", outcome),"\n")
            cat("\n")
          }
          #f1 and f2 are same as f0
          f1=paste(f0) #no covars
          f2=paste(f0) #no sig covars


        }else{

          # f1=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", covariate_list)
          f1=paste(f0, "+", covariate_list)

          # browser()

          m1=svyglm(noquote(f1), design=s)
          if(cutoff==weights_percentile_cutoff){
            cat(paste0("Covariate model results for effects of ", exposure, " on ", outcome,  " using ", cutoff_label),"\n")
            print(summary(m1))
            cat("\n")
          }
          models[["m1"]]<-m1

          #determining which covariates are significant
          sig_covars=as.data.frame(summary(m1)$coefficients)
          sig_covars=sig_covars[sig_covars$`Pr(>|t|)`<0.05,]
          sig_covars=rownames(sig_covars)
          sig_covars=sig_covars[!grepl(c("Intercept"), sig_covars)]
          sig_covars=sig_covars[!grepl(c(exposure), sig_covars)]
          sig_covars=c(as.character(sig_covars))

          #in the case of no significant covariates
          if (length(sig_covars)!=0){
            f2=paste(f0, "+", paste(sig_covars, sep="", collapse=" + "))
            if(cutoff==weights_percentile_cutoff){
              cat(paste0("The only significant covariate(s) for this model are: ", paste(sig_covars, sep="", collapse=" , ")),"\n")
              cat("\n")
            }

          }else{
            f2=f0
          }

          m2=svyglm(noquote(f2), design=s)
          if(cutoff==weights_percentile_cutoff){
            cat(paste0("Final covariate model results for effects of ", exposure, " on ", outcome," using ", cutoff_label),"\n")
            print(summary(m2))
            cat("\n")
          }
          models[["m2"]]<-m2
        }

        #getting interactions
        interactions=paste(apply(combn(exp_epochs,2), 2, paste, sep="", collapse=":"), sep="", collapse=" + ")
        f3=paste(f2, "+", paste(interactions, sep="", collapse=" + "))

        #testing for interactions
        # f3=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", sig_covars, "+", interactions)
        m3=svyglm(noquote(f3), design=s)
        if(cutoff==weights_percentile_cutoff){
          cat(paste0("Interaction model results for effects of ", exposure, " on ", outcome,  " using ", cutoff_label),"\n")
          print(summary(m3))
          cat("\n")
        }
        models[["m3"]]<-m3


        #determining which interactions are signficant
        sig_ints=as.data.frame(summary(m3)$coefficients)
        sig_ints=sig_ints[sig_ints$`Pr(>|t|)`<0.05,]
        sig_ints=rownames(sig_ints)
        sig_ints=sig_ints[grepl(c(":"), sig_ints)]

        #creating FINAL model
        if (is.na(sig_ints[1])){
          f4=f2 #equal to model with any covariates
        }else{
          if(cutoff==weights_percentile_cutoff){
            cat(paste0("The only significant interactions(s) for this model are: ", paste(sig_ints, sep="", collapse=" , ")), "\n")
            cat("\n")
          }

          f4=paste(f2, "+", paste(sig_ints, sep="", collapse=" + "))

        }

        # cat("The final model is:", "/n")
        m4=svyglm(noquote(f4), design=s)
        if(cutoff==weights_percentile_cutoff){
          cat(paste0("Final model results for effects of ", exposure, " on ", outcome,  " using ", cutoff_label),"\n")
          print(summary(m4))
          cat("\n")
        }
        models[["m4"]]<-m4


        all_models[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]<-models

        models=NULL
        m0=NULL
        m1=NULL
        m2=NULL
        m3=NULL
        m4=NULL

      }

    }

    cat("\n")

    file_label=ifelse(cutoff==weights_percentile_cutoff, "original", "sensitivity checks")

    saveRDS(all_models, file = paste(paste0(home_dir, "msms/", file_label, "/all_exposure-outcome_models_weight_cutoff_", cutoff, ".rds", sep="")))


  }
  # cat(paste0("All models with ", file_label, " cutoff weights have been saved as an .rds object in the 'msms/", file_label, "/' folder"),"\n")
  cat("\n")
  cat("\n")

  write.csv(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.csv"))
  cat("USER ALERT: Sensitivity check models using two alternate weight truncation values have been saved in the 'msms/sensitivity checks/' folder")
  cat("\n")
  cat("A new data file has been saved as a .csv file in the in the 'msms' folder","\n")


  return(all_models)

}
