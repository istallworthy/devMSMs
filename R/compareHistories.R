
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
  # # data=read.csv(paste0(home_dir, 'msms/data_for_msms.csv'))
  # all_data=readRDS(paste0(home_dir, 'msms/data_for_msms.rds'))
  # all_data=unname(all_data)
  # a=do.call(rbind.data.frame, all_data)

  # data=readRDS(paste0(home_dir, "final weights/values/imp_data_w_t.rds"))
  # data=complete(data, 1)
  # test=mice::as.mids(all_data, .imp=".imp")

  #creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")


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

  cat("USER ALERT: please inspect the following exposure history (uncorrected) comparisons for each exposure-outcome pair using user-specified original weight cutoff values:", "\n")
  cat("\n")

  #gets esimated marginal predictions
  preds <- lapply(1:length(fits), function(y){ #cutoffs
    d=fits[[y]]
    cutoff=all_cutoffs[y]
    cutoff_label=ifelse(cutoff==weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))

    lapply(d, function(f){
      final_model=f
      formula=f$formula
      parameters=gsub(")", "", formula[3])
      parameters=gsub(" ", "", paste0(unlist(strsplit(parameters, "\\+"))))

      #gathering epoch information for each exposure for deriving betas
      epoch_info=as.data.frame(rep(exposure, length(epochs)))
      epoch_info$time=epochs
      epoch_info$low=NA
      epoch_info$high=NA
      #cycling through epochs to find hi and lo values of the exposure for each epoch based on user-specified values
      for (t in 1:length(epochs)){
        var_name=paste(exposure, epochs[t], sep="_")
        epoch_info$low[t]=as.numeric(quantile(final_model$data[,var_name],probs=lo_cutoff, na.rm=T))
        epoch_info$high[t]=as.numeric(quantile(final_model$data[,var_name],probs= hi_cutoff, na.rm=T))
      }

      #gather high and low values for each exposure epoch
      d=data.frame(e=exp_epochs,
                   l=epoch_info$low,
                   h=epoch_info$high)
      d$v=paste(d$l, d$h, sep=",")
      d$z=lapply(1:nrow(d), function(x){c(as.numeric(unlist(strsplit(unlist(strsplit(d$v[x], " ")), "\\,")))) })
      args=d$z
      names(args)=(c(d$e)) #, "model"))

      #from Noah
      # Compute predictions; need to give custom class
      p=marginaleffects::avg_predictions(f,
                                         # newdata=marginaleffects::datagridcf(args),
                                         variables = args,
                                         # vcov="HC3",
                                         wts= "(weights)")
      class(p) <- c("pred_custom", class(p))
      p
    })
  })
  names(preds)<-all_cutoffs

  # sprintf(list(paste(paste0(d$e, " = %s"), collapse=", "),
  #      noquote(paste(paste0("out$", d$e), collapse=", "))))

  # Create custom tidy method; not called directly but used in mice::pool()
  tidy.pred_custom <- function(x, ...) {
    out <- NextMethod("tidy", x)
    # out$term <- sprintf("treat = %s, married = %s", out$treat, out$married)
    out$term <- sprintf("HOMEETA1_Infancy = %s, HOMEETA1_Toddlerhood = %s, HOMEETA1_Childhood = %s",
                        out$HOMEETA1_Infancy, out$HOMEETA1_Toddlerhood, out$HOMEETA1_Childhood)
    out
  }

  # Pool results; issue getting NA
  preds_pool<- lapply(preds, function(y){ mice::pool(mice::as.mira(y)) |> summary()
  })

  # Pairwise comparisons; don't need to use custom class
  comps <- lapply(preds, function(y){
    lapply(y, function(p) {
      browser()
    p |> marginaleffects::hypotheses("pairwise")
  })
  })


    #datagrd cf (from Noah 4/7)
    #datagridcf() generates "counter-factual" data frames, by replicating the entire dataset once for every combination of predictor values supplied by the user.
    # mod=lm(`EF_avg_perc.58`~ HOMEETA1_Infancy+ HOMEETA1_Toddlerhood + HOMEETA1_Childhood , weights = `HOMEETA1-EF_avg_perc_0.95_weight_cutoff`, data=as.data.frame(data_for_model_with_weights_cutoff[[1]]))
    #


    ####
    # marginaleffects::datagrid(model=final_model,
    #                          HOMEETA1_Infancy=c(-0.87125,0.48375),
    #                          HOMEETA1_Toddlerhood=c(-0.403875,1.123375),
    #                          HOMEETA1_Childhood=c(-0.02475,1.60475),
    #                          # wts=`HOMEETA1-EF_avg_perc_0.95_weight_cutoff`
    #                          grid_type = "counterfactual"
    #           )







    # #
    # g=do.call(marginaleffects::datagrid, args, quote=F)
    #
    #
    # #new data weights based on real weights
    # g$weights=seq(from=min(data_for_model_with_weights_cutoff[[1]][, paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff")]), #`HOMEETA1-EF_avg_perc_0.95_weight_cutoff`
    #               to=max(data_for_model_with_weights_cutoff[[1]][, paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff")]),
    #               by=max(data_for_model_with_weights_cutoff[[1]][, paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff")])/(nrow(g)))
    #
    # # datagrid(model=mice::as.mira(final_model),
    # #          g)
    #
    # #average predicted marginal means for all combos of high and low
    # pred_mean <- marginaleffects::avg_predictions(mice::as.mira(final_model),
    #                                               newdata=g,
    #                                               by=d$e,
    #                                               wts="weights"
    # )
    #




    comparisons=hypotheses(pred_mean, "pairwise")




    pred_mean=as.data.frame(pred_mean)

    #getting histories
    history=matrix(data=NA, nrow=nrow(pred_mean), ncol=nrow(d))
    for (i in 1:ncol( pred_mean[,colnames(pred_mean)[colnames(pred_mean) %in% d$e]])){
      temp=pred_mean[,colnames(pred_mean)[colnames(pred_mean) %in% d$e]]
      history[,i]=ifelse(round(temp[,i],3) == round(d[i,"l"],3), "l", "h")
    }
    history=as.data.frame(history)
    colnames(history)=gsub(paste0(exposure, "_"), "", d$e)
    pred_mean=cbind(history, pred_mean)
    pred_mean=pred_mean[,!colnames(pred_mean) %in% c("rowid", "type")]



    #makes table of average estimates
    folder_label=ifelse(cutoff==weights_percentile_cutoff, "original/", "sensitivity checks/")
    sink(paste0(home_dir, "msms/estimated means/", folder_label, "estimated_mean_", cutoff, ".html"))
    stargazer::stargazer(pred_mean,type="html", digits=2, column.labels = colnames(pred_mean),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/estimated means/", folder_label, "estimated_mean_", cutoff, ".html"))
    sink()



    #apply user-specified correction
    comparisons$p_vals_corr=stats::p.adjust(comparisons$p.value, method=method)

    #getting histories
    history=matrix(data=NA, nrow=nrow(comparisons), ncol=1)
    for (i in 1:nrow(comparisons)){
      temp=comparisons$term[i]
      pair=lapply (1:2, function(y){
        a=sapply(strsplit(temp, " - "), "[", y)
        his=lapply(1:nrow(d), function(z){
          ifelse(round(as.numeric(gsub("[^0-9.-]", "", sapply(strsplit(a, "\\,"), "[", z))),3) == round(d[z,"l"],3), "l", "h")
        })
      })
      history[i,1]=paste(sapply(pair, paste, collapse = "-"), collapse=" vs ")
    }
    comparisons=cbind(history, comparisons)


    #filters by user spc
    if(!is.na(reference)){
      comparisons=comparisons[gsub(" ", "", sapply(strsplit(comparisons$history, "vs"), "[",1)) ==reference,]
    }
    if(!is.na(comps)){
      comparisons=comparisons[gsub(" ", "", sapply(strsplit(comparisons$history, "vs"), "[",2)) %in% comps,]
    }


    #add dose
    dose_a=stringr::str_count(sapply(strsplit(comparisons$history, "vs"), "[", 1), dose_level)
    dose_b=stringr::str_count(sapply(strsplit(comparisons$history, "vs"), "[", 2), dose_level)
    comparisons=cbind(data.frame(dose=gsub(" ", " vs ", cbind(paste(dose_a, dose_b)))), comparisons)


    #finding significant comparisons
    sig_comparisons=comparisons[comparisons$p_vals_corr<0.05,]
    significant_comparisons[[cutoff]] <- sig_comparisons


    folder_label=ifelse(cutoff==weights_percentile_cutoff, "original/", "sensitivity checks/")

    # browser()
    cat(paste0("USER ALERT: please inspect the followig significant comparisons for models with weights truncated at ", cutoff, " :"), "\n")
    cat(knitr::kable(sig_comparisons, format='pipe'), sep="\n")
    cat("\n")

    #save out table for all contrasts with old and corrected p-values
    sink(paste0(home_dir, "msms/estimated means/", folder_label, "cutoff_",cutoff , "_comparisons.doc", sep=""))
    stargazer::stargazer(comparisons, type="html", digits=2, column.labels = colnames(comparisons),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/estimated means/", folder_label, "cutoff_",cutoff , "_comparisons.doc", sep=""))
    sink()



    # sink(paste0(home_dir, "msms/estimated means/", folder_label, "cutoff_",cutoff , "_comparisons.doc", sep=""))
    # stargazer::stargazer(comparisons, type="html", digits=2, column.labels = colnames(comparisons),summary=FALSE, rownames = FALSE, header=FALSE,
    #                      out=paste0(home_dir, "msms/linear hypothesis testing/", folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_presentation_table.doc", sep=""))
    # sink()
    # # marginaleffects::avg_comparisons(mice::as.mira(final_model),
    #                               # variables=list(HOMEETA1_Infancy=c(-0.87125,0.48375),
    #                               #                HOMEETA1_Toddlerhood=c(-0.403875,1.123375),
    #                               #                HOMEETA1_Childhood=c(-0.02475,1.60475)),
    #                               by=colnames(g)[1:3],
    #                               hypotheses(as.matrix(pred_mean), "pairwise"))
    #                               # cross=F,
    #                               # comparison= "difference",
    #                               # newdata =g)


    # marginaleffects::avg_predictions(mice::as.mira(fits),
    #                                  # variables="HOMEETA1_Infancy", by="HOMEETA1_Infancy",
    #                                  newdata=datagrid(HOMEETA1_Infancy=c(-0.87125,0.48375),
    #                                                     HOMEETA1_Toddlerhood=c(-0.403875,1.123375),
    #                                                     HOMEETA1_Childhood=c(-0.02475,1.60475)),
    #                                                     # grid_type = "counterfactual"),
    #                                  # type="response",
    #                                  by=c("HOMEETA1_Infancy","HOMEETA1_Toddlerhood","HOMEETA1_Childhood")
    #                                  )

    # #Frrom Noah
    #     data("lalonde_mis", package ="cobalt")
    #     m <- mice::mice(lalonde_mis, print = F)
    #     W <- MatchThem::weightthem(treat ~ age + educ + race + re74 + re75, m)
    #     #> Estimating weights     | dataset: #1 #2 #3 #4 #5
    #     Wd <- MatchThem::complete(W, "all")
    #     s <- survey::svydesign(~1, weights = ~weights, data = mitools::imputationList(Wd))
    #     fits <- with(s, survey::svyglm(re78 ~ treat + re75))
    #     marginaleffects::avg_predictions(mice::as.mira(fits), variables = "treat", by = "treat")
    #     #> Warning in get.dfcom(object, dfcom): Infinite sample size assumed.
    #     #>
    #     #>  treat Estimate Std. Error     t Pr(>|t|) 2.5 % 97.5 %
    #     #>      0     6411        343 18.68   <0.001  5738   7084
    #     #>      1     7127        731  9.75   <0.001  5694   8560
    #     #>
    #     #> Columns: treat, estimate, std.error, df, statistic, p.value, conf.low, conf.high
    #



    #
    #      pooled_pred <- pool(as.mira(test))
    #      summary(pooled_pred)
    #
    #
    #
    #
    #
    #     pred= marginaleffects::avg_predictions(fits,
    #                                        newdata=datagridcf(HOMEETA1_Infancy=c(-0.87125,0.48375) ,
    #                                                               HOMEETA1_Toddlerhood=c(-0.403875,1.123375),
    #                                                               HOMEETA1_Childhood=c(-0.02475,1.60475)),
    #                                        by=c("HOMEETA1_Infancy","HOMEETA1_Toddlerhood","HOMEETA1_Childhood"))
    #
    #      pred |> marginaleffects::hypotheses("pairwise")

    #
    #       if (cutoff==weights_percentile_cutoff){
    #         cat(paste0("The uncorrected difference between the effects of ", exposure, " at ", comparison, " compared to ", reference, " on ", outcome, " has a p-value of ", linear_hypothesis$`Pr(>F)`)[2], "\n")
    #       }
    #       comparisons[[paste0(reference, " vs. ", comparison)]] <-linear_hypothesis
    #       param_info[[paste0(reference, " vs. ", comparison)]] <-list(ref=ref_parameters_betas,
    #                                                                   comp=comp_parameters_betas)

    #   # write.csv(epoch_info, paste0(home_dir,"msms/history_betas_", exposure, "_", outcome, "_", comparison, "-vs-", reference, ".csv"))
    # }

    # history_comparisons[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]<- comparisons
    # parameter_beta_info[[paste0(exposure, "-", outcome,  "_cutoff_", cutoff)]] <- param_info

    # }

    # cat("\n")
    # }


    # saveRDS(history_comparisons, file = paste(paste0(home_dir, "msms/linear hypothesis testing/", folder_label, "all_linear_hypothesis_tests_cutoff_", cutoff, ".rds", sep="")))
    # saveRDS(parameter_beta_info, file = paste(paste0(home_dir, "msms/linear hypothesis testing/", folder_label, "all_linear_hypothesis_betas_parameters_cutoff_", cutoff, ".rds", sep="")))


    cat("\n")
    cat("Please see the 'msms/linear hypothesis testing/original' folder for csv files detailing the hi/lo beta values for each comparison","\n")
    cat("Please see the 'msms/linear hypothesis testing/sensitivity checks' folder for sensitivity check values","\n")

    return(history_comparisons)
  }



