
#' Compare exposure histories
#' This code uses the best-fitting model for each exposure-outcome pair to compare the effects of user-specified reference and comparison histories of exposure on outcome u sing linear hypothesis testing
#' @param object msm object that contains all relevant user inputs
#' @param best_models best-fitting models outcome from assessModel
#' @importFrom gtools permutations
#' @importFrom car linearHypothesis
#' @seealso [assessModel()]
#' @return history_comparisons lists of linear hypothesis tests
#' @examples compareHistories(object, best_models)
compareHistories <-function(object, data_for_model_with_weights_cutoff, all_models){

  home_dir=object$home_dir
  exposure=object$exposure
  exposure_epochs=object$exposure_epochs
  outcome=object$outcome
  outcome_time_pt=object$outcome_time_pt
  hi_cutoff=object$hi_cutoff
  lo_cutoff=object$lo_cutoff
  factor_covariates=object$factor_covariates
  reference=object$reference
  comps=object$comparisons
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))
  method=object$mc_method
  dose_level=object$dose_level

  #creating list of all cutoff values to cycle through
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)


  if (length(outcome_time_pt)>1){
    stop('This function is designed only for a single time point outcome')}

  epochs=exposure_epochs$epochs
  #creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")
  exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")

  #error checking
  if (hi_cutoff>1 | hi_cutoff<0){
    stop('Please select hi_cutoff between 0 and 1 in the msm object')}
  if (lo_cutoff>1 | lo_cutoff<0){
    stop('Please select lo_cutoff between 0 and 1 in the msm object')}

  if (!is.na(reference)){
    if (sum(exposure_levels %in% reference)==0){
      stop(paste0('If you wish to specify a reference event in the msm object, please select a valid exposure history from the following list ', paste0(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"), sep=" ", collapse=" ")))}
  }
  # reference=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")[length(apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"))]


  #if no comparison is specified by the user, compare to all histories aside from the reference
  # if (comparisons==""){
  #   comp_histories=exposure_levels[!exposure_levels %in% reference]
  # }else {
  if (!is.na(comps)){
    if (sum(exposure_levels %in% comps)==0){
      stop(paste0('If you wish to specify comparison(s) in the msm object, please select a valid history/valid histories from the following list ', paste0(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"), sep=" ", collapse=" ")))}
    comp_histories=exposure_levels[exposure_levels %in% comps]
  }

  # #if no reference is set by user, set to low exposure at all time points
  # if (reference==""){
  # }else {
  #   if (sum(exposure_levels %in% reference)==0){
  #     stop('Please select a valid reference history in the msm object')}
  # }


  #final list of all comparisons for all exposure-outcome pairs
  # history_comparisons=list()
  #
  # #final list of all betas and parameters for all comparisons
  # parameter_beta_info=list()

  cat("USER ALERT: please inspect the following exposure history comparisons using user-specified original weight cutoff values:", "\n")
  cat("\n")

  fits=all_models
  ints=gsub(" ", "", as.character(unlist(strsplit(as.character(unlist(fits[[1]][[1]]$terms)), "\\+"))))
  ints=ifelse(sum(grepl(":", ints))>0, 1, 0)
  #gathering epoch information for each exposure for deriving betas
  epoch_info=as.data.frame(rep(exposure, length(epochs)))
  epoch_info$time=epochs
  epoch_info$low=NA
  epoch_info$high=NA
  #stacking imputed data to compute hi/lo vals
  long_dat<-lapply(fits, function(x){lapply(x, function(y){y$data})})[[1]]
  long_dat=do.call("rbind", long_dat)
  #cycling through epochs to find hi and lo values of the exposure for each epoch based on user-specified values
  for (t in 1:length(epochs)){
    var_name=paste(exposure, epochs[t], sep="_")
    epoch_info$low[t]=as.numeric(quantile(long_dat[,var_name],probs=lo_cutoff, na.rm=T))
    epoch_info$high[t]=as.numeric(quantile(long_dat[,var_name],probs= hi_cutoff, na.rm=T))
  }


  # #gather high and low values for each exposure epoch
  # if (ints==0){
    d=data.frame(e=paste(epoch_info[[1]],epoch_info[[2]], sep="_"),
               l=epoch_info$low,
               h=epoch_info$high)

    d$v=paste(d$l, d$h, sep=",")
    d$z=lapply(1:nrow(d), function(x){c(as.numeric(unlist(strsplit(unlist(strsplit(d$v[x], " ")), "\\,")))) })
    args=d$z
    names(args)=(c(d$e))
    # }




  # #if the user specified reference and comparison groups: have to specify both (make warning)
  # if (!is.na(reference) & !is.na(comps)){
  #   if (length(unlist(strsplit(reference, "-")))!=nrow(d)){
  #     stop(paste0("If you wish to provide a reference event, please provide a single string of 'h' or 'l' (separated by a '-') equal in length to the number of exposure epochs (",
  #                 nrow(d), "). For example: ", paste(rep("l", nrow(d)), collapse="-")))
  #   }
  #   #gets reference values
  #   ref_vals=sapply(1:length(unlist(strsplit(reference, "-"))), function(x){d[x,unlist(strsplit(reference, "-"))[x]]})
  #   comp_vals=t(sapply(1:length(comps), function(y){
  #     lapply(1:length(unlist(strsplit(comps[y], "-"))), function(x){d[x,unlist(strsplit(comps[y], "-"))[x]]})
  #     }))
  #
  #   to_pred=as.data.frame(rbind(ref_vals, comp_vals))
  #   to_pred=sapply(to_pred, as.numeric)
  #   to_pred=as.data.frame(to_pred)
  #   colnames(to_pred)<-d$e
  #   # makes custom set of values to predict
  #   args=as.list(as.data.frame(to_pred))
  #
  #   marginaleffects::avg_predictions(f, newdata=datagridcf(args),
  #                                    wts= f$weights)
  #   #matching ref w/each comp
  #   all_vals=expand.grid(paste(ref_vals, collapse=" "), apply(comp_vals, 1, paste, collapse=" ", sep=" "))
  #
  # }

  #STEP 1: Estimated marginal predictions
  #gets esimated marginal predictions
  preds <- lapply(1:length(fits), function(y){ #goes through different weight truncation cutoff values (for sensitivity tests)
    c=fits[[y]]
    lapply(c, function(f){ #goes through imputations
      final_model=f
      p=marginaleffects::avg_predictions(f, variables = args,
                                         wts= "(weights)")
      class(p) <- c("pred_custom", class(p))
      p
    })
  })
  names(preds)<-all_cutoffs


  # paste0(paste(names(args), collapse=" = %s, "), " = %s")
  # out[,names(args)]
  #
  # temp=out[,names(args)]
  # test=unclass(temp)
  # a=sapply(test, function(x){
  #   test[[1]]
  # })

  # out$term=sprintf(paste0(paste(names(args), collapse=" = %s, "), " = %s"),
  #                  a)
  #                  # unclass(temp)#c(asplit(temp,2))
  #                  # )
  # # a=do.call(sprintf, c(fmt=paste0(paste(names(args), collapse=" = %s, "), " = %s"),
  # #                    as.list(test)))
  # a=do.call(sprintf, c(fmt=paste0(paste(names(args), collapse=" = %s, "), " = %s"),
  #                      as.list(unclass(out[,names(args)]))))
  #
  # browser()

  #STEP 2: Create custom tidy method--need to automate this
  # ; not called directly but used in mice::pool()
  #this does not work
  tidy.pred_custom <<- function(x, ...) {
    # args <<-args
    out <- NextMethod("tidy", x)
    out$term <- do.call(sprintf, c(paste0(paste(names(args), collapse=" = %s, "), " = %s"),
                                  as.list(unclass(out[,names(args)]))))
    out
  }


  # tidy.pred_custom <- function(x, ...) {
  #   out <- NextMethod("tidy", x)
  #   # out$term <- sprintf("treat = %s, married = %s", out$treat, out$married)
  #   out$term <- sprintf("HOMEETA1_Infancy = %s, HOMEETA1_Toddlerhood = %s, HOMEETA1_Childhood = %s",
  #                       out$HOMEETA1_Infancy, out$HOMEETA1_Toddlerhood, out$HOMEETA1_Childhood)
  #   out
  # }
  # tidy.pred_custom <- function(x, ...) {
  #   out <- NextMethod("tidy", x)
  #   # out$terrm<-do.call(sprintf, c(fmt=paste0(paste(names(args), collapse=" = %s, "), " = %s"),
  #   #                               as.list(unclass(out[,names(args)]))))
  #   out$term<-sprintf(paste0(paste(names(args), collapse=" = %s, "), " = %s"),
  #                     out[,"ESETA1_Infancy"], out[,"ESETA1_Toddlerhood"], out[,"ESETA1_Childhood"])
  #   out
  # }

  # library(tidyr)
  #this works
  # tidy.pred_custom <- function(x, ...) {
  #   out <- NextMethod("tidy", x)
  #   # out$term <- sprintf("treat = %s, married = %s", out$treat, out$married)
  #   out$term <- sprintf("InRatioCor_Infancy = %s, InRatioCor_Toddlerhood = %s, InRatioCor_Childhood = %s",
  #                       out$InRatioCor_Infancy, out$InRatioCor_Toddlerhood, out$InRatioCor_Childhood)
  #   out
  # }


  #STEP 3: pooling predicted estimates --seem to be the same for all cutoff values
  # Pool results
  preds_pool<- lapply(preds, function(y){ mice::pool(mice::as.mira(y)) |> summary()
  })

  # test=preds_pool[[1]]
  # test[1]=round(test[1],2)


  #adding histories to preds_pool
  history=matrix(data=NA, nrow=nrow(preds_pool[[1]]), ncol=1) #get histories from first element
  for (i in 1:nrow(preds_pool[[1]])){
    vals=as.numeric(sapply(strsplit(unlist(strsplit(as.character(preds_pool[[1]]$term[i]), "\\,")), "="), "[",2))
    history[i]=as.data.frame(paste(ifelse(round(vals,3)==round(d$l,3), "l", "h"), collapse="-"))
  }
  history=unlist(history)
  history=rep(list(history), length(preds_pool))
  preds_pool<-Map(cbind, preds_pool, history = history)

  #adding dose to preds_pool
  doses=stringr::str_count(history[[1]], dose_level)
  doses=rep(list(doses), length(preds_pool))
  preds_pool<-Map(cbind, preds_pool, dose = doses)

  cat("Below are the pooled average predictions by history:")
  print(preds_pool)


  #makes table of average estimates and saves out per cutoff value
  lapply(1:length(preds_pool), function(x){
    y=preds_pool[[x]]
    sink(paste0(home_dir, "msms/estimated means/estimated_means_", names(preds_pool)[x],".html"))
    stargazer::stargazer(as.data.frame(y),type="html", digits=2, column.labels = colnames(y),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/estimated means/estimated_means_", names(preds_pool)[x],".html"))
    sink()
  })



  #STEP 4: THIS CODE APPEARS TO RUN OK
  # Pairwise comparisons; don't need to use custom class
  comps <- lapply(preds, function(y){
    lapply(y, function(p) {
      p |> marginaleffects::hypotheses("pairwise")
    })
  })


  #STEP 5: pool comparison values
  #pool summary
  comps_pool<-lapply(comps, function(x){mice::pool(x) |> summary()})


  # browser()
  #apply user-specified correction to pooled comparisons
  corr_p<-lapply(comps_pool, function(x){
    stats::p.adjust(x$p.value, method=method)
  })
  # corr_p=list(corr_p)
  comps_pool<-Map(cbind, comps_pool, p.value_corr= corr_p)

  #adding histories to comps_pool
  history=matrix(data=NA, nrow=nrow(comps_pool[[1]]), ncol=1)
  for (i in 1:nrow(comps_pool[[1]])){
    temp=comps_pool[[1]]$term[i]
    temp=as.character(temp)
    pair=lapply (1:2, function(y){
      a=sapply(strsplit(temp, " - "), "[", y)
      his=lapply(1:nrow(d), function(z){
        ifelse(round(as.numeric(gsub("[^0-9.-]", "", sapply(strsplit(a, "\\,"), "[", z))),3) == round(d[z,"l"],3), "l", "h")
      })
    })
    history[i,1]=paste(sapply(pair, paste, collapse = "-"), collapse=" vs ")
  }
  history=rep(list(history), length(comps_pool))
  comps_pool<-Map(cbind, comps_pool, history = history)


  #adding dose to comps_pool
  dose_a=stringr::str_count(sapply(strsplit(history[[1]], "vs"), "[", 1), dose_level)
  dose_b=stringr::str_count(sapply(strsplit(history[[1]], "vs"), "[", 2), dose_level)
  dose_count=data.frame(dose=gsub(" ", " vs ", paste(dose_a, dose_b)))
  dose_count=rep(list(dose_count), length(comps_pool))
  comps_pool<-Map(cbind, comps_pool, dose_count = dose_count)

  cat("Below are the pooled comparisons by history which are saved in the 'msms/constrasts folder':")
  print(comps_pool)

  #makes table of comparisons and save out for each cutoff value
  lapply(1:length(comps_pool), function(x){
    y=comps_pool[[x]]
    sink(paste0(home_dir, "msms/contrasts/contrasts_", names(comps_pool)[x],".html"))
    stargazer::stargazer(as.data.frame(y),type="html", digits=2, column.labels = colnames(y),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/contrasts/contrasts_", names(comps_pool)[x],".html"))
    sink()
  })



  # #filters by user spc
  # if(!is.na(reference)){
  #   comparisons=comparisons[gsub(" ", "", sapply(strsplit(comparisons$history, "vs"), "[",1)) ==reference,]
  # }
  # if(!is.na(comps)){
  #   comparisons=comparisons[gsub(" ", "", sapply(strsplit(comparisons$history, "vs"), "[",2)) %in% comps,]
  # }

  # folder_label=ifelse(cutoff==weights_percentile_cutoff, "original/", "sensitivity checks/")

  # browser()
  cat(paste0("USER ALERT: please inspect the following comparisons for models with weights truncated at ", weights_percentile_cutoff, " :"), "\n")
  cat(knitr::kable(comps_pool[names(comps_pool)==weights_percentile_cutoff], format='pipe', digits=2), sep="\n")
  cat("\n")

  cat("\n")
  # cat("Please see the 'msms/linear hypothesis testing/estimated means' folder for csv files detailing the hi/lo beta values for each comparison","\n")
  # cat("Please see the 'msms/linear hypothesis testing/sensitivity checks' folder for sensitivity check values","\n")

  return(preds_pool)
}



