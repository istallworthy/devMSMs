
#' Fit weighted model
#' This code fits a weighted marginal structural model to examine the effects of different exposure histories on outcomes
#' @param home_dir path to home directory for the project
#' @param ID person-level identifier in your dataset
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param exposure_epochs data frame containing epochs and correponding time points for exposure histories
#' @param outcomes list of variables that represent your outcomes of interest
#' @param outcome_time_pts list of time points for outcomes
#' @param data_for_model_with_weights_cutoff dataset with truncated weights see truncateWeights
#' @param unbalanced_covariates_for_models unbalanced covariates see assessBalance
#' @importFrom survey svydesign
#' @importFrom survey svyglm
#' @return marginal structural models
#' @export
#' @seealso [truncateWeights()], [asesssBalance()]
#' @examples fitModel(ID, exposures, exposure_epochs, outcomes, outcome_time_pts, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models, hi_cutoff=75, lo_cutoff=25)

fitModel <- function(home_dir, ID, exposures, exposure_epochs, outcomes, outcome_time_pts, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models){

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


    #creates initial survey design type specifying ID, data, and weights to use
    s=svydesign(ids=data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==ID],
                data=data_for_model_with_weights_cutoff,
                weights=data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)==paste0(exposure, "_weight_cutoff")])

    #lists out exposure-epoch combos
    exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")


    #cycles through each outcome
    models=list()

    for (y in 1:length(outcomes)){
      outcome=outcomes[y]

      #creates form for baseline model including averaged exposure at each epoch as a predictor
      f0=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "))

      #fits initial baseline model
      m0=svyglm(noquote(f0), design=s)
      print(paste0("Baseline model results for effects of ", exposure, " on ", outcome))
      summary(m0)
      models[["m0"]]<-m0

      #adding in covariates that did not fully balance when creating the weights
      covariate_list= unbalanced_covariates_for_models[[exposure]]

      if (covariate_list[1]==""){
        print(paste0("There are no unbalanced covariates to include in the model of effects of ", exposure, " on ", outcome))

        #f1 and f2 are same as f0
        f1=paste(f0) #no covars
        f2=paste(f0) #no sig covars

        #form for interactions
        # f3=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", interactions)

      }else{

        # f1=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", covariate_list)
        f1=paste(f0, "+", covariate_list)

        m1=svyglm(noquote(f1), design=s)
        print(paste0("Covariate model results for effects of ", exposure, " on ", outcome))
        summary(m1)
        models[["m1"]]<-m1

        #determining which covariates are signfiicant
        sig_covars=as.data.frame(summary(m1)$coefficients)
        sig_covars=sig_covars[sig_covars$`Pr(>|t|)`<0.05,]
        sig_covars=rownames(sig_covars)
        sig_covars=sig_covars[!grepl(c("Intercept"), sig_covars)]
        sig_covars=sig_covars[!grepl(c(exposure), sig_covars)]
        sig_covars=c(as.character(sig_covars))
        print(paste0("The only significant covariate(s) for this model are: ", paste(sig_covars, sep="", collapse=" , ")))

        #final covariate model retaining only significant covariates
        # f2=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", sig_covars)
        f2=paste(f0, "+", paste(sig_covars, sep="", collapse=" + "))

        m2=svyglm(noquote(f2), design=s)
        print(paste0("Final covariate model results for effects of ", exposure, " on ", outcome))
        summary(m2)
        models[["m2"]]<-m2
        # anova(m0, m2)
      }

        #form for interactions
        # f3=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", sig_covars, "+", interactions)

      #getting interactions
      interactions=paste(apply(combn(exp_epochs,2), 2, paste, sep="", collapse=":"), sep="", collapse=" + ")
      f3=paste(f2, "+", paste(interactions, sep="", collapse=" + "))

      #testing for interactions
      # f3=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", sig_covars, "+", interactions)
      m3=svyglm(noquote(f3), design=s)
      summary(m3)
      models[["m3"]]<-m3


      #determining which interactions are signficant
      sig_ints=as.data.frame(summary(m3)$coefficients)
      sig_ints=sig_ints[sig_ints$`Pr(>|t|)`<0.05,]
      sig_ints=rownames(sig_ints)
      sig_ints=sig_ints[grepl(c(":"), sig_ints)]
      print(paste0("The only significant covariate(s) for this model are: ", paste(sig_ints, sep="", collapse=" , ")))

      #creating FINAL model
      if (is.na(sig_ints[1])){
        f4=f2 #equal to model with any covariates
      }else{

        f4=paste(f2, "+", paste(sig_ints, sep="", collapse=" + "))

      }

      m4=svyglm(noquote(f4), design=s)
      summary(m4)
      models[["m4"]]<-m4


      all_models[[paste0(exposure, "-", outcome)]]<-models

      models=NULL
      m0=NULL
      m1=NULL
      m2=NULL
      m3=NULL
      m4=NULL

    }

  }

  saveRDS(all_models, file = paste(paste0(home_dir, "msms/all_exposure-outcome_models.rds", sep="")))
  write.csv(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.csv"))

  return(all_models)

}
