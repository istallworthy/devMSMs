
#' Compare exposure histories
#' This code uses the best-fitting model for each exposure-outcome pair to compare the effects of user-specified reference and comparison histories of exposure on outcome u sing linear hypothesis testing
#' @param home_dir path to home directory for the project
#' @param exposure_epochs data frame containing epochs and correponding time points for exposure histories
#' @param hi_cutoff integer for percentile considered "high" for exposure
#' @param lo_cutoff integer for percentile considered "low" for exposure
#' @param outcomes list of variables that represent your outcomes of interest
#' @param outcome_time_pts list of time points for outcomes
#' @param best_models best-fitting models outcome from assessModel
#' @param comparisons list of histories to compare to the reference history
#' @param reference reference history to compare comparison histories
#' @importFrom gtools permutations
#' @importFrom car linearHypothesis
#' @seealso [assessModel()]
#' @return history_comparisons lists of linear hypothesis tests
#' @examples compareHistories(home_dir, exposures, exposure_epochs, hi_cutoff=.75, lo_cutoff=.25, outcome, outcome_time_pts, best_models, comparisons, reference)
compareHistories <-function(home_dir, exposures, exposure_epochs, hi_cutoff=.75, lo_cutoff=.25, outcome, outcome_time_pts, best_models, comparisons="", reference=""){

  epochs=exposure_epochs$epochs
  data=read.csv(paste0(home_dir, 'msms/data_for_msms.csv'))

  #creates permutations of high ("h") and low ("l") levels of exposures for each exposure epoch
  exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")


  #error checking
  if (hi_cutoff>1 | hi_cutoff<0){
    stop('Please select hi_cutoff between 0 and 1')}
  if (lo_cutoff>1 | lo_cutoff<0){
    stop('Please select lo_cutoff between 0 and 1')}
  if (sum(exposure_levels %in% reference)==0){
    stop('Please select a valid reference history from the list')}


  #if no comparison is specified by the user, compare to all histories aside from the reference
  if (comparisons==""){
    comp_histories=exposure_levels[!exposure_levels %in% reference]
  }else {
    if (sum(exposure_levels %in% comparisons)==0){
      stop('Please select a valid comparison history from the list')}
    comp_histories=exposure_levels[exposure_levels %in% comparisons]
  }

  #if no reference is set by user, set to low exposure at all time points
  if (reference==""){
    reference=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")[length(apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"))]
  }else {
    if (sum(exposure_levels %in% reference)==0){
      stop('Please select a valid reference history from the list')}
  }


  #final list of all comparisons for all exposure-outcome pairs
  history_comparisons=list()

  #final list of all betas and parameters for all comparisons
  parameter_beta_info=list()



  #cycles through outcomes
  for (y in 1:length(outcomes)){
    outcome=outcomes[y]


    #cycling through exposures
    for (x in 1:length(exposures)){
      comparisons=list()
      param_info=list()

      exposure=exposures[x]

      #gets best-fitting model formula
      formula=best_models[[paste0(exposure, "-", outcome)]]$formula
      final_model=best_models[[paste0(exposure, "-", outcome)]]


      #identifying model parameters
      parameters=sapply(strsplit(formula, "~"),  "[", 2)
      parameters=as.data.frame(strsplit(parameters, " +"))
      parameters=parameters[parameters !="+"]


      #gathering epoch information for each exposure for deriving betas
      epoch_info=as.data.frame(rep(exposure, length(epochs)))
      epoch_info$time=epochs
      epoch_info$low=NA
      epoch_info$high=NA


      #cycling through epochs to find hi and lo values of the exposure for each epoch
      for (t in 1:length(epochs)){
        var_name=paste(exposure, epochs[t], sep="_")
        epoch_info$low[t]=as.numeric(quantile(data[,var_name],probs= lo_cutoff, na.rm=T))
        epoch_info$high[t]=as.numeric(quantile(data[,var_name],probs= hi_cutoff, na.rm=T))
      }


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
            a=sapply(strsplit(ref_parameters_betas[z,1], ":"), "[", 1) #first term
            b=sapply(strsplit(ref_parameters_betas[z,1], ":"), "[", 2) #second term
            beta=epoch_info[which(grepl(a, paste(epoch_info[,1], epoch_info[,2], sep="_"))),"ref_betas"]* #find product of a*b
              epoch_info[which(grepl(b, paste(epoch_info[,1], epoch_info[,2], sep="_"))),"ref_betas"]
            ref_parameters_betas$beta[z]=beta

          }else{
            #populates main effect beta based on history sequence from epoch_info
            ref_parameters_betas$beta[z]=epoch_info[which(ref_parameters_betas[z,1]==paste(epoch_info[,1], epoch_info[,2], sep="_")), "ref_betas"]
          }
        }else{ #otherwise it is a covariate, use grand mean as beta
          ref_parameters_betas$beta[z]=mean(data[,colnames(data)[colnames(data)==ref_parameters_betas[z,1]]], na.rm=T)

        }
      }
      #create reference form with betas and parameters
      ref_form=paste(ref_parameters_betas$beta, ref_parameters_betas$parameter, sep="*", collapse=" + ")





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
              a=sapply(strsplit(comp_parameters_betas[z,1], ":"), "[", 1) #first term
              b=sapply(strsplit(comp_parameters_betas[z,1], ":"), "[", 2) #second term
              beta=epoch_info[which(grepl(a, paste(epoch_info[,1], epoch_info[,2], sep="_"))),"comp_betas"]* #find product of a*b
                epoch_info[which(grepl(b, paste(epoch_info[,1], epoch_info[,2], sep="_"))),"comp_betas"]
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
        linear_hypothesis=car::linearHypothesis(final_model, paste(ref_form, comp_form, sep=" = "), test="F")
        print(paste0("The difference between the effects of ", exposure, " at ", comparison, " compared to ", reference, " on ", outcome, " has a p-value of ", linear_hypothesis$`Pr(>F)`)[2])

        comparisons[[paste0(reference, " vs. ", comparison)]] <-linear_hypothesis
        param_info[[paste0(reference, " vs. ", comparison)]] <-list(ref=ref_parameters_betas,
                                                              comp=comp_parameters_betas)

        # write.csv(epoch_info, paste0(home_dir,"msms/history_betas_", exposure, "_", outcome, "_", comparison, "-vs-", reference, ".csv"))
      }

      history_comparisons[[paste0(exposure, "-", outcome)]]<- comparisons
      parameter_beta_info[[paste0(exposure, "-", outcome)]] <- param_info
    }
  }

  saveRDS(history_comparisons, file = paste(paste0(home_dir, "msms/all_linear_hypothesis_tests.rds", sep="")))
  saveRDS(parameter_beta_info, file = paste(paste0(home_dir, "msms/all_linear_hypothesis_betas_parameters.rds", sep="")))

  print("Please see the 'msms' folder for csv files detailing the hi/lo beta values for each comparison")

  return(history_comparisons)

}