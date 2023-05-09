#' Assess success of covariate balancing
#'
#' #Assesses how well balance was achieved for each of the covariates/potential confounds in relation to each of the exposure,
#' and returns a list of unbalanced covariates for each exposure to add to future models
#'
#' @param object msm object that contains all relevant user inputs
#' @param weights_models optional if createWeights has just run or will read in weights models from local storage
#' @param just_made_weights="no" no=weights are not assigned in global environment and are instead saved locally
#' @return list of unbalanced_covariates_for_models for each exposure
#' @seealso [CBPS::balance()] for more on the balance function
#' @seealso [msmObject()] for more on the weights_models param
#' @seealso [createWeights()] for more on the weights_models param
#' @export
#' @importFrom CBPS balance
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @importFrom cobalt love.plot
#' @importFrom dplyr group_by
#' @examples assessBalance(object, weights_models=list(), just_made_weights="no")
assessBalance <- function (object, forms, data_for_model_with_weights, weighted, histories=1){

  home_dir=object$home_dir
  m=object$m
  exposure=object$exposure
  exposure_time_pts=object$exposure_time_pts
  outcome=object$outcome
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates
  weights_method=object$weights_method

  # browser()

  .call <- match.call()
  form_type<-.call[[3]]
  form_name<-as.character(form_type)

  assign(as.character(form_type), forms)
  # browser()

  # forms=short_forms
  exposure_type=ifelse(length(unique(data_for_model_with_weights[[1]][,paste0(exposure,".", exposure_time_pts[1])]))<3, "binary", "continuous")

  #gathering balance stats for all imputed datasets
  bal_stats= lapply(seq(m), function(i){
    # browser()
    calcBalStats(object, data_for_model_with_weights, get(form_type), form_name, exposure, outcome, k=i, weighted, histories)
  })


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
    suppressMessages(ggplot2::ggsave(lp, filename=paste0(home_dir, "balance/plots/", form_name, "_", exposure, "_all_imps_", exposure_time_pt,"_",weights_method, "_summary_balance_plot.jpeg"),
                                      ))
    cat(paste0("Balance summary plot for ", form_name, "_", exposure,  " all imputations at time ", exposure_time_pt, " for weighting method ",weights_method," has now been saved in the '", "balance/plots/' folder."), "\n")
  })



  #save out correlations/std mean differences
  sink(paste0(home_dir, "balance/comparison values/", form_name,"_",weights_method, "_all_post-balance_associations.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/comparison values/", form_name, "_",weights_method,"_all_post-balance_assocations.html"))
  sink()

  cat(paste0("USER ALERT: Check 'balance/comparison values/' folder for a table of all post-balance correlations or standardized mean differences and averaged across imputed datasets"),"\n")
  # cat("\n")

  unbalanced_covars= unbalanced_covars%>%
    dplyr::filter(balanced_avg==0)
  # browser()

  #renames factor covariates
  unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]

  unbalanced_constructs=sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1)[!duplicated(sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1))]

  cat("\n")
  cat(paste0("USER ALERT: Averaging across all imputed datasets for exposure ", exposure, " using the ", weights_method, " weighting method and the ", form_name,
             ", the following ", nrow(unbalanced_covars) ," covariates across time points out of ",
             length(tot_covars), " total (", round((nrow(unbalanced_covars)/length(tot_covars))*100,2), "%) spanning ",
             length(unbalanced_constructs), " unique constructs (out of ", length(tot_cons), ") remain imbalanced with an average absolute residual value correlation/std mean difference in relation to ",
             exposure, " of ", round(mean(abs(unbalanced_covars$avg_bal)),2), " (range=",
             round(min(unbalanced_covars$avg_bal),2), "-", round(max(unbalanced_covars$avg_bal),2), ") : "), "\n")
  # print(unbalanced_covars)
  cat("\n")
  cat(knitr::kable(unbalanced_covars, caption=paste0("Imbalanced Covariates using ", weights_method, " and ",form_name), format='pipe'),  sep="\n")
  cat("\n")

  #save out only imbalanced covariates
  sink(paste0(home_dir, "balance/comparison values/", "_",weights_method, "_", form_name, "_all_imbalanced_covariates.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/comparison values/",  "_",weights_method, "_", form_name, "_all_imbalanced_covariates.html"))
  sink()

  write.csv(unbalanced_covars, paste0(home_dir, "balance/comparison values/", exposure, "_unbalanced_stat_summary.csv"))
  cat(paste0("Imbalanced covariates for ", exposure, " using the ", form_name, " and ", weights_method, " weights method have been saved in the 'balance/comparison values/' folder"), "\n")


  return(bal_stats)

}






