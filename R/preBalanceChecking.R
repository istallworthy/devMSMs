
#' Examining initial imbalance for each exposure, each time point, each covariate for each history using the Jackson method
#'
#' @param object msm object that contains all relevant user inputs
#' @param  wide_long_datasets from formatForWeights()
#' @param  all_forms from createForms()
#' @param histories binary indicator of whether to print history sample distributions
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 guide_axis
#' @importFrom ggplot2 ggsave
#' @importFrom stargazer stargazer
#' @importFrom dplyr bind_rows
#' @export
#' @examples preBalanceChecking(object, wide_long_datasets, all_forms, histories=1)
#'
preBalanceChecking <- function(object, wide_long_datasets, all_forms, histories=1){

  exposure=object$exposure
  outcome=object$outcome
  m=object$m
  exposure_time_pts=object$exposure_time_pts
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates
  home_dir=object$home_dir

  #getting form type from input
  .call <- match.call()
  form_type=.call[[4]]
  form_name=as.character(form_type)

  exposure_type=ifelse(length(unique(wide_long_datasets[[1]][,paste0(exposure,".", exposure_time_pts[1])]))<3, "binary", "continuous")

  cat("USER ALERT: The following statistics display covariate imbalance at each exposure time point prior to weighting, using full formulas that reflect covariates at all lagged time points.", "\n")
  cat("\n")

  #running balance stats function, unweighted, on each imputed dataset
  bal_stats<- lapply(1:m, function(k){
    calcBalStats(object, wide_long_datasets, all_forms, form_name, exposure, outcome, k=k, weighted=0, histories)
  })

  #save out balance stats for each imputed dataset
  bal_stats_all_imp <-bal_stats %>%
    dplyr::bind_rows(.id = "imputation")
  bal_stats_all_imp=bal_stats_all_imp[order(bal_stats_all_imp$covariate),]
  write.csv(bal_stats_all_imp, paste0(home_dir, "pre balance/", exposure, "_all_imps_balance_stat_summary.csv"), row.names = F)
  cat("Pre balance statistics for each imputed dataset have now been saved in the 'pre balance/' folder", "\n")

  #gathering imbalanced covariate statistics to average across imputed datasets for final list/assessment of imbalanced covariates
  #averaging across imputed datases
  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(bal_stats, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=bal_stats[[1]]$exp_time,
                               covar_time=bal_stats[[1]]$covar_time,
                               covariate=bal_stats[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0)) #assessing new averaged bal stat in relation to balance threshold

  #getting totals
  tot_covars=sapply(strsplit(bal_stats[[1]]$covariate, "\\."), "[", 1)
  tot_cons=tot_covars[!duplicated(tot_covars)] #total domains/constructs
  x_lab=ifelse(exposure_type=="continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")

  cat("\n")
  cat(paste0("*** Averaging Across All Imputations ***"), "\n")

  #make love plot to summarize imbalance at each exposure time point
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
                     panel.border = ggplot2::element_rect(fill = NA, color = "black"),
                     plot.background = ggplot2::element_blank(),
                     plot.title=ggplot2::element_text(size=10),
                     legend.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank(),
                     legend.position="none")+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if(nrow(temp)>40){ #stagger covariate labels if there are many to ease viewing
      lp <- lp+ ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge=2))
    }
    suppressMessages(ggplot2::ggsave(lp,
                                     filename=paste0(home_dir, "pre balance/plots/", form_name, "_", exposure, "_all_imps_", exposure_time_pt, "_summary_pre_balance_plot.jpeg")))
    cat(paste0("USER ALERT: A balance summary plot for ", form_name, " ", exposure,  " averaged across all imputations at time ", exposure_time_pt, " has now been saved in the '", "pre balance/plots/' folder."), "\n")
  })


  #save out all prebal correlations/std mean differences
  sink(paste0(home_dir, "pre balance/",  exposure, "-", outcome, "_all_pre-balance_associations.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "pre balance/",  exposure, "-", outcome,"_all_pre-balance_assocations.html"))
  sink()
  cat(paste0("USER ALERT: Check 'pre balance/comparison values/' folder for a table of all pre-balance correlations or standardized mean differences and averaged across imputed datasets."),"\n")


  #renames factor covariates
  unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]
  unbalanced_constructs=sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1))]

  #saving out all prebal associations
  write.csv(unbalanced_covars, paste0(home_dir, "pre balance/", exposure, "_prebalance_stat_summary.csv"), row.names = F)
  cat(paste0("All associations between exposure and covariates for ", form_name, " ",exposure, "-", outcome, " prior to balancing have been saved in the 'pre balance/' folder"), "\n")
  cat("\n")

  #finding all imbalanced variables
  unbalanced_covars= unbalanced_covars%>%
    dplyr::filter(balanced_avg==0)
  cat("\n")

  cat(paste0("USER ALERT: Before weighting, averaging across all imputed datasets for exposure ", exposure, " using the ", form_name, ", the following ", nrow(unbalanced_covars) ," covariates across time points out of ",
             length(tot_covars), " total (",round(nrow(unbalanced_covars)/length(tot_covars)*100,2)  ,",%) spanning ",
             length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round(length(unbalanced_constructs)/length(tot_cons)*100,2),"%) are imbalanced with a remaining average absolute value correlation/std mean difference in relation to ",
             exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)),2), " (range=",
             round(min(unbalanced_covars$avg_bal),2), "-", round(max(unbalanced_covars$avg_bal),2), ") : "), "\n")
  cat(knitr::kable(unbalanced_covars, caption="Imbalanced Covariates Before Weighting", format='pipe'),  sep="\n")

  #save out only imbalanced covariates
  sink(paste0(home_dir, "pre balance/",  exposure, "-", outcome,"_all_imbalanced_covariates.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "pre balance/",  exposure, "-", outcome,"_all_imbalanced_covariates.html"))
  sink()

}

