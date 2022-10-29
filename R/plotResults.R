#' Plot results from history comparisons
#' Code to plot predicted values of the different exposure histories by history and colored by dosage of exposure
#' @param object msm object that contains all relevant user inputs
#' @param best_models outcome from assessModels
#' @seealso [assessModels()]
#' @return plots
#' @importFrom ggplot2 ggplot
#' @examples plotResults(object, best_models)
plotResults <- function(object, best_models){

  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  exposure_labels=object$exposure_labels
  outcome_labels=object$outcome_labels

  if (exposure_labels==""){
    exposure_labels=exposures #default is to use
  }else {
    if (length(exposure_labels) != length(exposures)){
      stop('Please provide labels for all exposures (in the same order as listed in the exposures field) in the msm object')
    }
  }

  if (outcome_labels==""){
    outcome_labels=outcomes #default is to use
  }else {
    if (length(outcome_labels) != length(outcomes)){
      stop('Please provide labels for all outcomes (in the same order as listed in the outcome field) in the msm object')
    }
  }


  for (x in 1:length(best_models)){
    exp_out=names(best_models)[x]
    exposure=sapply(strsplit(exp_out, "-"), "[",1)
    outcome=sapply(strsplit(exp_out, "-"), "[",2)

    final_model=best_models[[exp_out]]

    parameter_beta_info=readRDS(paste0(home_dir, "msms/linear hypothesis testing/all_linear_hypothesis_betas_parameters.rds"))

    parameters=parameter_beta_info[[exp_out]][[1]]$ref$parameter


    #getting predicted value info for ref and comparison histories
    param_betas_ref=lapply(parameter_beta_info[[exp_out]], function(x) x$ref)
    param_betas_ref=as.data.frame(t(as.data.frame(param_betas_ref[[1]]$beta))) #pulls first one bc ref is always the same for all contrasts
    colnames(param_betas_ref)=parameters
    param_betas_ref=cbind(sapply(strsplit(names(lapply(parameter_beta_info[[exp_out]], function(x) x$ref))[1], "vs."), "[",1),  param_betas_ref)
    colnames(param_betas_ref)=c("seq", parameters)
    param_betas_comp=lapply(parameter_beta_info[[exp_out]], function(x) x$comp)
    param_betas_comp=as.data.frame(t(as.data.frame(lapply(param_betas_comp, function(x) x$beta))))
    param_betas_comp=cbind(sapply(strsplit(row.names(param_betas_comp), "vs."), "[",2), param_betas_comp)
    colnames(param_betas_comp)=c("seq", parameters)

    all_param_betas=rbind(param_betas_ref, param_betas_comp)
    all_param_betas[,2:ncol(all_param_betas)]=lapply(all_param_betas[,2:ncol(all_param_betas)], function(x) as.numeric(x))
    temp=as.data.frame(predict(final_model, newdata= all_param_betas[,2:ncol(all_param_betas)], interval="confidence"))
    all_param_betas=cbind(all_param_betas, temp)
    all_param_betas$low_ci=all_param_betas$link-(1.96*all_param_betas$SE) #finds CIs
    all_param_betas$upp_ci=all_param_betas$link+(1.96*all_param_betas$SE) #finds CIs

    all_param_betas=as.data.frame(cbind(sapply(strsplit(exp_out, "-"), "[",1), sapply(strsplit(exp_out, "-"), "[",2),
                          all_param_betas))

    colnames(all_param_betas)=c("exposure", "outcome", "seq", parameters, "fit", "se", "lwr", "upr")

    #create dataset for predicting
    plot_data=all_param_betas
    plot_data$seq=gsub(".", "", as.character(plot_data$seq), fixed = T)
    plot_data$seq=gsub("-", "", as.character(plot_data$seq), fixed = T)

    #finds doses of "h"
    plot_data$dose=as.numeric(lapply(regmatches(plot_data$seq, gregexpr("h", plot_data$seq)), function(x) length(x)))
    plot_data$dose=as.factor(plot_data$dose)
    plot_data$seq=as.factor(plot_data$seq)
    plot_data=plot_data[order(plot_data$dose),]

    ggplot2::ggplot(aes(x=fit, y=seq, color=dose), data=plot_data)+
      ggplot2::geom_point()+
      ggplot2::scale_color_manual(values=c("blue4", "darkgreen", "darkgoldenrod", "red2"))+
      ggplot2::geom_errorbarh(xmin = plot_data$fit-plot_data$se, xmax = plot_data$fit+plot_data$se, height=0.6)+
      ggplot2::xlab(paste0("Predicted ", outcome_labels[which(outcome %in% outcomes)], " Value"))+
      ggplot2::ylab(paste0(exposure_labels[which(exposure %in% exposures)], " Exposure History"))+
      ggplot2::xlim(min(plot_data$fit-plot_data$se)-2*sd(plot_data$fit-plot_data$se), max(plot_data$fit+plot_data$se)+2*sd(plot_data$fit+plot_data$se))+
      ggplot2::theme(text = ggplot2::element_text(size=18))+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
    ggplot2::ggsave(paste0(home_dir, "results figures/", sapply(strsplit(exp_out, "-"), "[",1), "-",sapply(strsplit(exp_out, "-"), "[",2), ".jpeg"), plot=ggplot2::last_plot())

  }
  cat("See the 'results figures' folder for graphical representations of results","\n")
}
