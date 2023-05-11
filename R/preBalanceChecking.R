
#examining initial imbalance for each exposure, each time point, each covariate for each history

preBalanceChecking <- function(object, wide_long_datasets, all_forms, histories=1){

  exposure=object$exposure
  outcome=object$outcome
  m=object$m
  exposure_time_pts=object$exposure_time_pts
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates
  home_dir=object$home_dir

  .call <- match.call()
  form_type=.call[[4]]
  form_name=as.character(form_type)

  exposure_type=ifelse(length(unique(wide_long_datasets[[1]][,paste0(exposure,".", exposure_time_pts[1])]))<3, "binary", "continuous")



  cat("The following statistics display covariate imbalance for each exposure at each exposure time point prior to weighting, using full formulas that reflect covariates at all lagged time points.
      Best practice for assessing balance for time-varying exposures (Jackson, 2016) first assesses imbalance for each exposure history
      (note that these histories are different from the user-specified ones used in the final model as they reflect median split high (1) and low (0) values at each time point) up until the exposure time point and then averages the balance statistics across histories.
      Below you will see the distributions of individuals in each history for each exposure time point that contribbute to the average balance statistics used to determine balance in the final table.", "\n")

  # forms=all_forms

  #running balance stats function
  bal_stats<- lapply(1:m, function(k){
    # browser()
    calcBalStats(object, wide_long_datasets, all_forms, form_name, exposure, outcome, k=k, weighted=0, histories)
  })
  # return(bal_summary)

  # browser()

  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(bal_stats, "[", "std_bal_stats"))))
  # browser()
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=bal_stats[[1]]$exp_time,
                               covar_time=bal_stats[[1]]$covar_time,
                               covariate=bal_stats[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0))

  #getting totals
  tot_covars=sapply(strsplit(bal_stats[[1]]$covariate, "\\."), "[", 1)
  tot_cons=tot_covars[!duplicated(tot_covars)]

  x_lab=ifelse(exposure_type=="continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")


  #make love plot per exposure time point
  sapply(seq_along(exposure_time_pts), function(i){
    exposure_time_pt=exposure_time_pts[i]

    temp=unbalanced_covars%>%dplyr::filter(exp_time==exposure_time_pt)

    labels=ifelse(temp$balanced_avg==0, temp$covariate, "")
    min_val=ifelse(min(temp$avg_bal)<0, min(temp$avg_bal)-0.1, balance_thresh-0.05)
    max_val=ifelse(max(temp$avg_bals)>0, max(temp$avg_bal)+0.1, balance_thresh+0.05)

    lp <- ggplot2::ggplot(aes(x =  temp$avg_bal, y = temp$covariate), data = temp) +
      ggplot2::geom_point(aes(y = as.factor(temp$covariate),
                              x = temp$avg_bal,
                              fill = "white", na.rm = TRUE,
                              alpha = 1))+
      ggplot2::geom_text(aes(label=labels, hjust=-0.2, vjust=0.2), size=1.5, color="red")+
      ggplot2::geom_vline(xintercept = balance_thresh,linetype = "dashed", color = "red")+
      ggplot2::geom_vline(xintercept = -balance_thresh,linetype = "dashed", color = "red")+
      ggplot2::xlab(x_lab)+
      ggplot2::ylab("Covariate")+
      ggplot2::xlim(min_val, max_val)+
      ggplot2::ggtitle(paste0(exposure, "(t=", exposure_time_pt, ") Balance"))+
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                     axis.text.x = ggplot2::element_text(color = "black"),
                     axis.text.y = ggplot2::element_text(color = "black"),
                     axis.text = ggplot2::element_text(size=8),
                     # axis.ticks.length.y = unit(0.85, "cm"),
                     panel.border = ggplot2::element_rect(fill = NA, color = "black"),
                     plot.background = ggplot2::element_blank(),
                     plot.title=ggplot2::element_text(size=10),
                     legend.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank(),
                     legend.position="none")+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if(nrow(temp)>40){ #stagger covariate labels if there are many
      lp <- lp+ ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge=2))
    }
    # browser()
    suppressMessages(ggplot2::ggsave(lp, filename=paste0(home_dir, "pre balance/plots/", form_name, "_", exposure, "_all_imps_", exposure_time_pt, "_summary_pre_balance_plot.jpeg"),
    ))
    cat(paste0("Balance summary plot for ", form_name, "_", exposure,  " all imputations at time ", exposure_time_pt, " has now been saved in the '", "pre balance/plots/' folder."), "\n")
  })



  #save out correlations/std mean differences
  sink(paste0(home_dir, "pre balance/",  exposure, "-", outcome, "_all_pre-balance_associations.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "pre balance/",  exposure, "-", outcome,"_all_pre-balance_assocations.html"))
  sink()

  cat(paste0("USER ALERT: Check 'pre balance/comparison values/' folder for a table of all pre-balance correlations or standardized mean differences and averaged across imputed datasets"),"\n")
  # cat("\n")

  unbalanced_covars= unbalanced_covars%>%
    dplyr::filter(balanced_avg==0)
  # browser()

  #renames factor covariates
  unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]

  unbalanced_constructs=sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1))]

  cat("\n")
  cat(paste0("USER ALERT: Before weighting, averaging across all imputed datasets for exposure ", exposure, ", the following ", nrow(unbalanced_covars) ," covariates across time points (out of ",
             length(tot_covars), " total) spanning ",
             length(unbalanced_constructs), " unique constructs (out of ", length(tot_cons), ") are imbalanced with an average absolute value correlation/std mean difference in relation to ",
             exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)),2), " (range=",
             round(min(unbalanced_covars$avg_bal),2), "-", round(max(unbalanced_covars$avg_bal),2), "), based on the ", form_name, " : "), "\n")
  # print(unbalanced_covars)
  cat(knitr::kable(unbalanced_covars, caption="Imbalanced Covariates Before Weighting", format='pipe'),  sep="\n")

  #save out only imbalanced covariates
  sink(paste0(home_dir, "pre balance/",  exposure, "-", outcome,"_all_imbalanced_covariates.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "pre balance/",  exposure, "-", outcome,"_all_imbalanced_covariates.html"))
  sink()

  write.csv(unbalanced_covars, paste0(home_dir, "pre balance/", exposure, "_prebalance_stat_summary.csv"))
  cat(paste0("Imbalanced covariates for ", form_name, " ",exposure, " prior to balancing have been saved in the 'pre balance/' folder"), "\n")


}

