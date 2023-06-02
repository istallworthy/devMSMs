#' Assess success of covariate balancing
#' #Assesses how well balance was achieved for each of the covariates/potential confounds in relation to each of the exposure,
#' and returns a list of unbalanced covariates for each exposure to add to future models
#'
#' @param object msm object that contains all relevant user inputs
#' @param forms formula for assessing balance
#' @param data_for_model_with_weights imputed data with weights
#' @param histories optional binary indicator of whether to print histories for bal stasts
#' @return list of unbalanced_covariates_for_models for each exposure
#' @seealso [msmObject()] for more on the weights_models param
#' @seealso [createWeights()] for more on the weights_models param
#' @export
#' @importFrom knitr kable
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
#' @importFrom dplyr bind_rows
#' @examples assessBalance(object, forms, data_for_model_with_weights, histories=1)
assessBalance <- function (object, forms, data_for_model_with_weights, histories=1){

  home_dir=object$home_dir
  m=object$m
  exposure=object$exposure
  exposure_time_pts=object$exposure_time_pts
  outcome=object$outcome
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates
  weights_method=object$weights_method


  #getting form type from input
  .call <- match.call()
  form_type<-.call[[3]]
  form_name<-as.character(form_type)
  assign(as.character(form_type), forms)

  exposure_type=ifelse(length(unique(data_for_model_with_weights[[1]][,paste0(exposure,".", exposure_time_pts[1])]))<3, "binary", "continuous")

  #gathering weighted balance stats for all imputed datasets
  bal_stats= lapply(seq(m), function(i){
    calcBalStats(object, data_for_model_with_weights, get(form_type), form_name, exposure, outcome, k=i, weighted=1, histories)
  })

  #save out balance stats for each imputed dataset
  bal_stats_all_imp <-bal_stats %>%
    dplyr::bind_rows(.id = "imputation")
  bal_stats_all_imp=bal_stats_all_imp[order(bal_stats_all_imp$covariate),]
  write.csv(bal_stats_all_imp, paste0(home_dir, "balance/",form_name, "_", exposure, "_all_imps_balance_stat_summary.csv"), row.names = F)
  cat("Balance statistics for each imputed dataset have now been saved in the 'balance/' folder", "\n")

  #gather imbalanced covariates from imputed datasets and averaging across imputations
  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(bal_stats, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=bal_stats[[1]]$exp_time,
                               covar_time=bal_stats[[1]]$covar_time,
                               covariate=bal_stats[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0)) #determining balance based on avg bal stats

  #getting totals
  tot_covars=sapply(strsplit(bal_stats[[1]]$covariate, "\\."), "[", 1)
  tot_cons=tot_covars[!duplicated(tot_covars)]

  x_lab=ifelse(exposure_type=="continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")

  cat("\n")
  cat(paste0("*** Averaging Across All Imputations ***"), "\n")

  #make love plot per exposure time point for bal stats averaged across imputed datasets
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
      ggplot2::ggtitle(paste0(exposure, "(t=", exposure_time_pt, ") Balance Averaged Across Imputed Datasets"))+
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
    if(nrow(temp)>40){ #stagger covariate labels if there are many
      lp <- lp+ ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge=2))
    }
    suppressMessages(ggplot2::ggsave(lp, filename=paste0(home_dir, "balance/plots/", form_name, "_", exposure, "_all_imps_", exposure_time_pt,"_",
                                                         weights_method, "_summary_balance_plot.jpeg")))
    cat(paste0("USER ALERT: Balance summary plot for ", form_name, "_", exposure,  " all imputations at time ",
               exposure_time_pt, " for weighting method ",weights_method," has now been saved in the 'balance/plots/' folder."), "\n")
  })

  #save out correlations/std mean differences
  sink(paste0(home_dir, "balance/comparison values/", form_name,"_",weights_method, "_all_post-balance_associations.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/comparison values/", form_name, "_",weights_method,"_all_post-balance_assocations.html"))
  sink()
  cat(paste0("USER ALERT: Check 'balance/comparison values/' folder for a table of all post-balance correlations/standardized mean differences and averaged across imputed datasets"),"\n")

  # browser()

  #comparing balance to prebalance stats if using full balance forms (which were used for prebalance assessment)
  if (form_type=="full_forms"){
    prebal<-read.csv(paste0(home_dir, "pre balance/", exposure, "_prebalance_stat_summary.csv"), row.names = NULL)
    colnames(prebal)[colnames(prebal)=="avg_bal"]<-"avg_pre_bal"
    colnames(prebal)[colnames(prebal)=="balanced_avg"]<-"pre_balanced_avg"
    comp=merge(unbalanced_covars, prebal, by=c("exposure", "exp_time","covar_time", "covariate"))
    cat("\n")
    cat(paste0("Prior to weighting using the full balancing formula, ", sum(comp$pre_balanced_avg), " out of ",nrow(comp), " (", round((sum(comp$pre_balanced_avg)/nrow(comp))*100,2), "%) covariates were balanced compared to ",
               sum(comp$balanced_avg),  " out of ",nrow(comp), " (", round((sum(comp$balanced_avg)/nrow(comp))*100,2), "%) after weighting. ",
               sum(comp$balanced_avg-comp$pre_balanced_avg==-1), " covariates were balanced prior to weighting and are now imbalanced after weighting:"), "\n")
    # colnames(comp)<-c("exposure", "exposure time point", "covariate time point", "covariate", "avg_bal_stat", "balanced", "avg_pre_bal_stat", "pre-balanced")
    cat(knitr::kable(comp[comp$balanced_avg-comp$pre_balanced_avg==-1,], format='pipe'),sep="\n")
  }

  #getting just imbalanced covars
  unbalanced_covars= unbalanced_covars%>%
    dplyr::filter(balanced_avg==0)

  #renames factor covariates
  unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]
  unbalanced_constructs=sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1))]

  cat("\n")
  cat(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ", weights_method, " weighting method and the ", form_name,
             ", the following ", nrow(unbalanced_covars) ," covariates across time points out of ",
             length(tot_covars), " total (", round((nrow(unbalanced_covars)/length(tot_covars))*100,2), "%) spanning ",
             length(unbalanced_constructs), " domains out of ", length(tot_cons), " (", round((length(unbalanced_constructs)/length(tot_cons))*100,2),
             " %) remain imbalanced with a remainingg average absolute value correlation/std mean difference in relation to ",
             exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)),2), " (range=",
             round(min(unbalanced_covars$avg_bal),2), "-", round(max(unbalanced_covars$avg_bal),2), "): "), "\n")
  cat("\n")
  cat(knitr::kable(unbalanced_covars, caption=paste0("Imbalanced Covariates using ", weights_method, " and ",form_name), format='pipe'),  sep="\n")
  cat("\n")

  #save out only imbalanced covariates
  sink(paste0(home_dir, "balance/comparison values/", form_name, "_",weights_method, "_all_imbalanced_covariates.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/comparison values/",  "_",weights_method, "_", form_name, "_all_imbalanced_covariates.html"))
  sink()

  write.csv(unbalanced_covars, paste0(home_dir, "balance/comparison values/", form_name, "_", exposure, "_unbalanced_stat_summary.csv"))
  cat(paste0("Imbalanced covariates for ", exposure, " using the ", form_name, " and ", weights_method, " weights method have been saved in the 'balance/comparison values/' folder"), "\n")

  return(bal_stats)
}






