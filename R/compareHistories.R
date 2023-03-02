
#' Compare exposure histories
#' This code uses the best-fitting model for each exposure-outcome pair to compare the effects of user-specified reference and comparison histories of exposure on outcome u sing linear hypothesis testing
#' @param object msm object that contains all relevant user inputs
#' @param best_models best-fitting models outcome from assessModel
#' @importFrom gtools permutations
#' @importFrom car linearHypothesis
#' @seealso [assessModel()]
#' @return history_comparisons lists of linear hypothesis tests
#' @examples compareHistories(object, best_models)
compareHistories <-function(object, data_for_model_with_weights_cutoff, all_models, imp_data_w_t){

  home_dir=object$home_dir
  exposure=object$exposure
  exposure_epochs=object$exposure_epochs
  outcome=object$outcome
  outcome_time_pt=object$outcome_time_pt
  hi_cutoff=object$hi_cutoff
  lo_cutoff=object$lo_cutoff
  reference=object$reference
  comparisons=object$comparisons
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))



  #creating list of all cutoff values to cycle through
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)


  if (length(outcome_time_pt)>1){
    stop('This function is designed only for a single time point outcome')}

  epochs=exposure_epochs$epochs
  # data=read.csv(paste0(home_dir, 'msms/data_for_msms.csv'))
  all_data=readRDS(paste0(home_dir, 'msms/data_for_msms.rds'))
  all_data=unname(all_data)
  a=do.call(rbind.data.frame, all_data)

  data=complete(imp_data_w_t, 1)
  test=mice::as.mids(all_data, .imp=".imp")

  #creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")


  #error checking
  if (hi_cutoff>1 | hi_cutoff<0){
    stop('Please select hi_cutoff between 0 and 1 in the msm object')}
  if (lo_cutoff>1 | lo_cutoff<0){
    stop('Please select lo_cutoff between 0 and 1 in the msm object')}
  if (sum(exposure_levels %in% reference)==0){
    stop(paste0('Please select a valid reference history in the msm object from the following list ', paste0(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"), sep=" ", collapse=" ")))}


  #if no comparison is specified by the user, compare to all histories aside from the reference
  if (comparisons==""){
    comp_histories=exposure_levels[!exposure_levels %in% reference]
  }else {
    if (sum(exposure_levels %in% comparisons)==0){
      stop(paste0('Please select a valid comparison history in the msm object from the following list ', paste0(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"), sep=" ", collapse=" ")))}
    comp_histories=exposure_levels[exposure_levels %in% comparisons]
  }

  #if no reference is set by user, set to low exposure at all time points
  if (reference==""){
    reference=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")[length(apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"))]
  }else {
    if (sum(exposure_levels %in% reference)==0){
      stop('Please select a valid reference history in the msm object')}
  }


  #final list of all comparisons for all exposure-outcome pairs
  history_comparisons=list()

  #final list of all betas and parameters for all comparisons
  parameter_beta_info=list()

  cat("USER ALERT: please inspect the following exposure history (uncorrected) comparisons for each exposure-outcome pair using user-specified original weight cutoff values:", "\n")
  cat("\n")

  for (cut in 1:length(all_cutoffs)){
    cutoff=all_cutoffs[cut]

    comparisons=list()
    param_info=list()

    # exposure=exposure[x]

    #gets best-fitting model formula
    # formula=best_models[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]$formula
    formula=all_models[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]
    formula=formula[[model]]
    # formula=formula$coef.names

    # final_model=best_models[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]
    final_model=as.formula(paste0(paste0(outcome, ".", outcome_time_pt), "~", paste(formula$coef.names[2:length(formula$coef.names)], sep=" ", collapse=" + ")))
    final_model=all_models[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]

    #identifying model parameters
    # parameters=sapply(strsplit(formula, "~"),  "[", 2)
    parameters=formula$coef.names
    parameters=as.data.frame(strsplit(parameters, " +"))
    parameters=parameters[parameters !="+"]

    #gathering epoch information for each exposure for deriving betas
    epoch_info=as.data.frame(rep(exposure, length(epochs)))
    epoch_info$time=epochs
    epoch_info$low=NA
    epoch_info$high=NA


    #cycling through epochs to find hi and lo values of the exposure for each epoch based on user-specified values
    for (t in 1:length(epochs)){
      var_name=paste(exposure, epochs[t], sep="_")
      epoch_info$low[t]=as.numeric(quantile(data[,var_name],probs=lo_cutoff, na.rm=T))
      epoch_info$high[t]=as.numeric(quantile(data[,var_name],probs= hi_cutoff, na.rm=T))
    }

    #DETERMINING REFERENCE HISTORY PARAMETER VALUES
    #finds reference betas based on specified reference history
    epoch_info$reference=sapply(strsplit(reference, "-"),"[")
    epoch_info$ref_betas=NA
    for (b in 1:nrow(epoch_info)){
      if(epoch_info$reference[b]=="l"){
        epoch_info$ref_betas[b]=epoch_info$low[b]
      } else if(epoch_info$reference[b]=="h"){
        epoch_info$ref_betas[b]=epoch_info$high[b]
      }
    }

    #finds reference beta values and formula based on parameters of best-fitting model
    ref_parameters_betas=as.data.frame(parameters[2:length(parameters)])
    names(ref_parameters_betas)="parameter"
    ref_parameters_betas$beta=NA

    #finding the reference beta values for each parameter
    for (z in 1:nrow(ref_parameters_betas)){
      #if it contains the exposure, it must be either a main effect or an interaction term
      if(grepl(exposure, ref_parameters_betas[z,1])){
        if(grepl(":", ref_parameters_betas[z,1])){ #if it is an interaction term
          int_num=length(unlist(strsplit(ref_parameters_betas[z,1], ":"))) #determine interaction order
          terms=as.data.frame(letters[seq(from=1, to=int_num)])

          for (i in 1:nrow(terms)){ #cycle through interaction terms
            terms[i, 2]=sapply(strsplit(ref_parameters_betas[z,1], ":"), "[", i) #gather each interaction term
            # b=sapply(strsplit(ref_parameters_betas[z,1], ":"), "[", 2) #second term
            terms[i, 3]=epoch_info[which(grepl(terms[i, 2], paste(epoch_info[,1], epoch_info[,2], sep="_"))),"ref_betas"] #find beta value
          }

          beta=prod(terms[,3]) #get product of all interaction term betas
          ref_parameters_betas$beta[z]=beta

        }else{
          #populates main effect beta based on history sequence from epoch_info
          ref_parameters_betas$beta[z]=epoch_info[which(ref_parameters_betas[z,1]==paste(epoch_info[,1], epoch_info[,2], sep="_")), "ref_betas"]
        }
      }else{ #if it does not contain exposure, it is a covariate, use grand mean as beta
        ref_parameters_betas$beta[z]=mean(data[,colnames(data)[colnames(data)==ref_parameters_betas[z,1]]], na.rm=T)

      }
    }
    #create reference form with betas and parameters
    ref_form=paste(ref_parameters_betas$beta, ref_parameters_betas$parameter, sep="*", collapse=" + ")




    #DETERMINING COMPARISON HISTORY/HISTORIES PARAMETER VALUES
    #cycles through comparison history list to compare each comparison history to the reference history
    for (c in 1:length(comp_histories)){
      comparison=comp_histories[c]

      #finds comparison betas based on specified comparison history
      epoch_info$comparison=sapply(strsplit(comparison, "-"),"[")
      epoch_info$comp_betas=NA
      for (b in 1:nrow(epoch_info)){
        if(epoch_info$comparison[b]=="l"){
          epoch_info$comp_betas[b]=epoch_info$low[b]
        } else if(epoch_info$comparison[b]=="h"){
          epoch_info$comp_betas[b]=epoch_info$high[b]
        }
      }

      #finds comparison beta values and formula based on parameters of best-fitting model
      comp_parameters_betas=as.data.frame(parameters[2:length(parameters)])
      names(comp_parameters_betas)="parameter"
      comp_parameters_betas$beta=NA

      #finding the comparison beta values for each parameter
      for (z in 1:nrow(comp_parameters_betas)){
        #if it contains the exposure, it must be either a main effect or an interaction term
        if(grepl(exposure, comp_parameters_betas[z,1])){
          if(grepl(":", comp_parameters_betas[z,1])){ #if it is an interaction term
            int_num=length(unlist(strsplit(comp_parameters_betas[z,1], ":"))) #determine interaction order
            terms=as.data.frame(letters[seq(from=1, to=int_num)])

            for (i in 1:nrow(terms)){ #cycle through interaction terms
              terms[i, 2]=sapply(strsplit(ref_parameters_betas[z,1], ":"), "[", i) #gather each interaction term
              terms[i, 3]=epoch_info[which(grepl(terms[i, 2], paste(epoch_info[,1], epoch_info[,2], sep="_"))),"ref_betas"] #find beta value
            }

            beta=prod(terms[,3]) #get product of all interaction term betas
            comp_parameters_betas$beta[z]=beta

          }else{
            #populates main effect beta based on comparison history sequence from epoch_info
            comp_parameters_betas$beta[z]=epoch_info[which(comp_parameters_betas[z,1]==paste(epoch_info[,1], epoch_info[,2], sep="_")), "comp_betas"]
          }
        }else{ #otherwise it is a covariate, use grand mean as beta
          comp_parameters_betas$beta[z]=mean(data[,colnames(data)[colnames(data)==comp_parameters_betas[z,1]]], na.rm=T)

        }
      }
      #create forms with betas and parameters
      comp_form=paste(comp_parameters_betas$beta, comp_parameters_betas$parameter, sep="*", collapse=" + ")



      #linear hypothesis test to compare the comparison and reference values
      linear_hypothesis=car::linearHypothesis(a, paste(ref_form, comp_form, sep=" = "), test="F")
      # mitml::testConstraints(m0,  constraints=paste(ref_form, comp_form, sep=" = "))

      if (cutoff==weights_percentile_cutoff){
        cat(paste0("The uncorrected difference between the effects of ", exposure, " at ", comparison, " compared to ", reference, " on ", outcome, " has a p-value of ", linear_hypothesis$`Pr(>F)`)[2], "\n")
      }
      comparisons[[paste0(reference, " vs. ", comparison)]] <-linear_hypothesis
      param_info[[paste0(reference, " vs. ", comparison)]] <-list(ref=ref_parameters_betas,
                                                                  comp=comp_parameters_betas)

      # write.csv(epoch_info, paste0(home_dir,"msms/history_betas_", exposure, "_", outcome, "_", comparison, "-vs-", reference, ".csv"))
    }

    history_comparisons[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]<- comparisons
    parameter_beta_info[[paste0(exposure, "-", outcome,  "_cutoff_", cutoff)]] <- param_info

    # }

    cat("\n")
    # }

    folder_label=ifelse(cutoff==weights_percentile_cutoff, "original/", "sensitivity checks/")

    saveRDS(history_comparisons, file = paste(paste0(home_dir, "msms/linear hypothesis testing/", folder_label, "all_linear_hypothesis_tests_cutoff_", cutoff, ".rds", sep="")))
    saveRDS(parameter_beta_info, file = paste(paste0(home_dir, "msms/linear hypothesis testing/", folder_label, "all_linear_hypothesis_betas_parameters_cutoff_", cutoff, ".rds", sep="")))
  }

  cat("\n")
  cat("Please see the 'msms/linear hypothesis testing/original' folder for csv files detailing the hi/lo beta values for each comparison","\n")
  cat("Please see the 'msms/linear hypothesis testing/sensitivity checks' folder for sensitivity check values","\n")

  return(history_comparisons)

}
