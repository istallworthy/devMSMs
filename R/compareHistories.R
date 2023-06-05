
#' Compare exposure histories
#' This code uses the best-fitting model for each exposure-outcome pair to compare the effects of user-specified reference and comparison histories of exposure on outcome u sing linear hypothesis testing
#' @param object msm object that contains all relevant user inputs
#' @param data_for_model_with_weights_cutoff imputed datasets with truncated weights
#' @param all_models fitted models for each imputed dataset
#' @param reference optional reference event for custom comparison
#' @param compare optional coparison event(s) for custom comparison
#' @importFrom gtools permutations
#' @importFrom marginaleffects avg_predictions
#' @importFrom marginaleffects hypotheses
#' @importFrom stringr str_count
#' @importFrom dplyr filter
#' @importFrom stargazer stargazer
#' @importFrom mice pool
#' @importFrom stats p.adjust
#' @importFrom knitr kable
#' @return preds_pool
#' @examples compareHistories(object, data_for_model_with_weights_cutoff, all_models, reference=NA, compare=NA)
compareHistories <-function(object, data_for_model_with_weights_cutoff, all_models){

  home_dir=object$home_dir
  exposure=object$exposure
  exposure_time_pts=object$exposure_time_pts
  exposure_epochs=object$exposure_epochs
  outcome=object$outcome
  outcome_time_pt=object$outcome_time_pt
  hi_cutoff=object$hi_cutoff
  lo_cutoff=object$lo_cutoff
  factor_covariates=object$factor_covariates
  reference=object$reference
  compare=object$comparisons
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))
  method=object$mc_method
  dose_level=object$dose_level

  #creating list of all cutoff values to cycle through
  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  #gathering epoch info
  epochs=exposure_epochs$epochs
  #creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")
  exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")

  exposure_type=ifelse(length(unique(data_for_model_with_weights_cutoff[[1]][[1]][,paste0(exposure,".", exposure_time_pts[1])]))<3, "binary", "continuous")

  #error checking
  if (hi_cutoff>1 | hi_cutoff<0){
    stop('In the msmObject, please select hi_cutoff between 0 and 1 in the msm object')}
  if (lo_cutoff>1 | lo_cutoff<0){
    stop('In the msmObject, please select lo_cutoff between 0 and 1 in the msm object')}

  if (!is.na(reference)){
    if (sum(exposure_levels %in% reference)==0){
      stop(paste0('If you wish to conduct custom comparisons, in the msmObject please select a valid reference history from the following list ', paste(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep=",", collapse="-"), sep=", ", collapse=", ")))}
    if(sum(is.na(compare))>0){
      stop(paste0("If you wish to conduct custom comparisons, in the msmObject please specify at least one valid comparison history from the following list ",  paste0(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep=",", collapse="-"), sep=" ", collapse=" "),
                  " otherwise, set the 'reference' and 'comparison' fields in the msmObject to 'NA' to conduct all comparisons."))
    }else if(reference %in% compare){
        stop("If you wish to make a custom comparison, please provide unique reference and comparison events in the msmObject")
      }
    }

  if (sum(is.na(compare))==0){
    if (sum(exposure_levels %in% compare)==0){
      stop(paste0('If you wish to specify comparison(s), please provide at least one comparison history from the following list ', paste0(apply(gtools::permutations(2, nrow(object$exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-"), sep=", ", collapse=" "),
           " otherwise, set the 'reference' and 'comparison' fields in the msmObject to 'NA' to conduct all comparisons."))}
    comp_histories=exposure_levels[exposure_levels %in% compare]
  }else{comp_histories=NA}


  fits=all_models
  ints=gsub(" ", "", as.character(unlist(strsplit(as.character(unlist(fits[[1]][[1]]$terms)), "\\+"))))
  ints=ifelse(sum(grepl(":", ints))>0, 1, 0)
  #gathering epoch information for each exposure for deriving betas
  epoch_info=as.data.frame(rep(exposure, length(epochs)))
  epoch_info$time=epochs
  epoch_info$low=NA
  epoch_info$high=NA


  #cycling through epochs to find hi and lo values of the exposure for each epoch based on user-specified values
  if (exposure_type=="continuous"){
    #stacking imputed data to compute hi/lo vals across all
    long_dat<-lapply(fits, function(x){lapply(x, function(y){y$data})})[[1]]
    long_dat=do.call("rbind", long_dat)
    for (t in 1:length(epochs)){
      var_name=paste(exposure, epochs[t], sep="_")
      epoch_info$low[t]=as.numeric(quantile(long_dat[,var_name],probs=lo_cutoff, na.rm=T))
      epoch_info$high[t]=as.numeric(quantile(long_dat[,var_name],probs= hi_cutoff, na.rm=T))
    }}
  if(exposure_type=="binary"){
    for (t in 1:length(epochs)){
      var_name=paste(exposure, epochs[t], sep="_")
      epoch_info$low[t]=0
      epoch_info$high[t]=1
    } }
  # #gather high and low values for each exposure epoch based on quantile values
    d=data.frame(e=paste(epoch_info[[1]],epoch_info[[2]], sep="_"),
               l=epoch_info$low,
               h=epoch_info$high)
    d$v=paste(d$l, d$h, sep=",")
    d$z=lapply(1:nrow(d), function(x){c(as.numeric(unlist(strsplit(unlist(strsplit(d$v[x], " ")), "\\,")))) })
    args=d$z #creating vector of each epoch and each corresponding h/l value
    names(args)=(c(d$e))



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


  #STEP 2: Create custom tidy method--need to automate this
  # ; not called directly but used in mice::pool()
  tidy.pred_custom <<- function(x, ...) {
    # args <<-args
    out <- NextMethod("tidy", x)
    out$term <- do.call(sprintf, c(paste0(paste(names(args), collapse=" = %s, "), " = %s"),
                                  as.list(unclass(out[,names(args)]))))
    out
  }


  #STEP 3: pooling predicted estimates --seem to be the same for all cutoff values
  # Pool results
  preds_pool<- lapply(preds, function(y){ mice::pool(mice::as.mira(y)) |> summary(digits=3)
  })

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


  #if the user specified reference and comparison groups, subset pred_pool for inspection and plotting
  if (!is.na(reference) & sum(is.na(comp_histories))==0){
    preds_pool<-lapply(preds_pool, function(y){
      as.data.frame(y)%>%dplyr::filter(history %in% c(reference, comp_histories))
    })
  }

  cat("\n")
  cat("Below are the pooled average predictions by user-specified history:") #not sure if we need to print this?
  print(preds_pool)

  #makes table of average estimates and saves out per cutoff value
  lapply(1:length(preds_pool), function(x){
    y=preds_pool[[x]]
    y$term<-sapply(y$term, function(x){ gsub(paste0(exposure, "_"), "", x)})
    sink(paste0(home_dir, "msms/estimated means/",exposure, "-", outcome, "_estimated_means_",
                names(preds_pool)[x],"_hi=", hi_cutoff, "_lo=",lo_cutoff, ".html"))
    stargazer::stargazer(as.data.frame(y),type="html", digits=4, column.labels = colnames(y),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/estimated means/", exposure, "-", outcome, "_estimated_means_",
                                    names(preds_pool)[x],"_hi=", hi_cutoff, "_lo=",lo_cutoff, ".html"))
    sink()
  })



  #STEP 4: CONDUCT COMPARISONS
  #if the user specified reference and comparison groups, create custom contrasts
  if (!is.na(reference) & sum(is.na(comp_histories))==0){
    #gets reference values
    ref_vals=sapply(1:length(unlist(strsplit(reference, "-"))), function(x){d[x,unlist(strsplit(reference, "-"))[x]]})
    comp_vals=t(sapply(1:length(comp_histories), function(y){
      lapply(1:length(unlist(strsplit(comp_histories[y], "-"))), function(x){d[x,unlist(strsplit(comp_histories[y], "-"))[x]]})
    }))
    # custom comparisons to reference; don't need to use custom class
    cus_comps=matrix(ncol=nrow(comp_vals), nrow=nrow(as.data.frame(preds[1][[1]])))
    #finds ref position in preds
    ref_pos=which(apply(as.data.frame(preds[[1]][1])[grepl(exposure, colnames(as.data.frame(preds[[1]][1])))], 1,paste, collapse=",")==paste0(ref_vals, collapse=","))
    cus_comps[ref_pos,]<- -1 #setting ref position to 1
    #finds comp position in preds
    for (x in 1:nrow(comp_vals)){
      c=paste(comp_vals[x,], collapse=",")
      cus_comps[which(apply(as.data.frame(preds[[1]][1])[grepl(exposure, colnames(as.data.frame(preds[[1]][1])))], 1,paste, collapse=",")
                      ==paste0(c, collapse=",")), x]<-1
    }
    if (nrow(comp_vals>1)){
      colnames(cus_comps)=paste0("(", paste0(paste0(ref_vals, collapse=","), ") - (", apply(comp_vals, 1, paste,collapse=",")), ")")
    } else{colnames(cus_comps)=paste0("(", paste0(paste0(ref_vals, collapse=","), ") - (", paste0(paste0(comp_vals, collapse=",")), ")"))}
    cus_comps[is.na(cus_comps)]<-0
    comps <- lapply(preds, function(y){ #custom contrasts
      lapply(y, function(p) {
        p |> marginaleffects::hypotheses(cus_comps)
      })
    })
  }

  #conduct all pairwise comparisons if no ref/comparisons were specified by user
  if (sum(is.na(reference) & is.na(comp_histories))==1){
    # Pairwise comparisons; don't need to use custom class
    comps <- lapply(preds, function(y){
      lapply(y, function(p) {
        p |> marginaleffects::hypotheses("pairwise")
      })
    })
  }


  #STEP 5: pool comparison values
  #pool summary
  comps_pool<-lapply(comps, function(x){mice::pool(x) |> summary()})

  #apply user-specified MC correction to pooled comparisons >1
  if (sum((is.na(reference) & is.na(comp_histories)) | length(comp_histories)>1)>0){
    cat("\n")
    cat(paste0("Conducting multiple comparison correction using the ", method, " method."), "\n")
    cat("\n")
    corr_p<-lapply(comps_pool, function(x){
      stats::p.adjust(x$p.value, method=method)
    })
    comps_pool<-Map(cbind, comps_pool, p.value_corr= corr_p)
  } else {
    cat(paste0("The user specified comparison only between ", reference, " and ", comp_histories, " so no correction for multiple comparisons will be implemented."), "\n")
  }

  #adding history labels to comps_pool
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

  cat("USER ALERT: please inspect the following exposure history comparison(s) for the user-specified original weight cutoff", "\n")
  cat("\n")
  cat("Below are the pooled comparisons by history which are saved in the 'msms/constrasts/' folder:", "\n")
  print(comps_pool)

  #makes table of comparisons and save out for each cutoff value
  lapply(1:length(comps_pool), function(x){
    y=comps_pool[[x]]
    sink(paste0(home_dir, "msms/contrasts/",  exposure, "-", outcome,"_contrasts_",
                names(comps_pool)[x],"_hi=", hi_cutoff, "_lo=",lo_cutoff,".html"))
    stargazer::stargazer(as.data.frame(y),type="html", digits=4, column.labels = colnames(y),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/", exposure, "-", outcome,"_contrasts_",
                                    names(comps_pool)[x],"_hi=", hi_cutoff, "_lo=",lo_cutoff,".html"))
    sink()
  })

  cat(paste0("USER ALERT: please inspect the following comparisons for models with weights truncated at ", weights_percentile_cutoff, " :"), "\n")
  cat(knitr::kable(comps_pool[names(comps_pool)==weights_percentile_cutoff], format='pipe', digits=2), sep="\n")
  cat("\n")
  cat("\n")

  return(preds_pool)
}



