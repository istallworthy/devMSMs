#' Plot results from history comparisons
#' Code to plot predicted values of the different exposure histories by history and colored by dosage of exposure
#' @param object msm object that contains all relevant user inputs
#' @param history_estimates pooled predicted values for marginal model
#' @return plots
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_errorbarh
#' @importFrom ggplot2 ggsave
#' @examples plotResults(object, history_estimates)
plotResults <- function(object, history_estimates){

  home_dir=object$home_dir
  exposure=object$exposure
  outcome=object$outcome
  exposure_labels=object$exposure_labels
  outcome_labels=object$outcome_labels
  colors=object$colors
  weights_percentile_cutoff=object$weights_percentile_cutoff
  dose_level=object$dose_level

  #error checking
  if (length(colors)>1 & length(colors)!=nrow(object$exposure_epochs)+1){
    stop(paste0('Please provide either: ',nrow(object$exposure_epochs)+1, ' different colors, a color palette, or leave this entry blank in the msm object'))}

  #cycling through pooled estimates for each exposure history implicated in comparisons
  lapply(1:length(history_estimates), function(x){
    cutoff=names(history_estimates)[x]
    comparisons=history_estimates[[x]]
    comparisons$term=gsub(paste0(exposure, "_"), "",comparisons$term)
    comparisons$low_ci=comparisons$estimate-(1.96*comparisons$std.error) #finds CIs
    comparisons$high_ci=comparisons$estimate+(1.96*comparisons$std.error) #finds CIs
    comparisons$history=as.factor(comparisons$history)
    comparisons$dose=as.factor(comparisons$dose)
    comparisons=comparisons[order(comparisons$dose),]

    if (length(colors)>1){ #if user input a list of colors
      ggplot2::ggplot(aes(x=estimate, y=history, color=dose), data=comparisons)+
        ggplot2::geom_point(size=5)+
        ggplot2::scale_color_manual(values=colors)+
        ggplot2::scale_y_discrete(limits=c(as.character(comparisons$history)), expand=c(0, 0.2))+
        ggplot2::geom_errorbarh(xmin = comparisons$estimate-comparisons$std.error,
                                xmax = comparisons$estimate+comparisons$std.error, height=0.6)+
        ggplot2::xlab(paste0("Predicted ", outcome_labels, " Value"))+
        ggplot2::ylab(paste0(exposure_labels, " Exposure History"))+
        ggplot2::xlim(min(comparisons$estimate-comparisons$std.error)-1*sd(comparisons$estimate-comparisons$std.error),
                      max(comparisons$estimate+comparisons$std.error)+1*sd(comparisons$estimate+comparisons$std.error))+
        ggplot2::theme(text = ggplot2::element_text(size=14))+
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
      suppressMessages(ggplot2::ggsave(paste0(home_dir, "results figures/", exposure, "-", outcome, " cutoff= ", cutoff, ".jpeg"), plot=ggplot2::last_plot()))
    }else{ #if user lists a palette
      ggplot2::ggplot(aes(x=estimate, y=history, color=dose), data=comparisons)+
        ggplot2::geom_point(size=5)+
        ggplot2::scale_color_brewer(palette=colors)+
        ggplot2::scale_y_discrete(limits=c(as.character(comparisons$history)), expand=c(0, 0.2))+
        ggplot2::geom_errorbarh(xmin = comparisons$estimate-comparisons$std.error,
                                xmax = comparisons$estimate+comparisons$std.error, height=0.6)+
        ggplot2::xlab(paste0("Predicted ", outcome_labels, " Value"))+
        ggplot2::ylab(paste0(exposure_labels, " Exposure History"))+
        ggplot2::xlim(min(comparisons$estimate-comparisons$std.error)-1*sd(comparisons$estimate-comparisons$std.error),
                      max(comparisons$estimate+comparisons$std.error)+1*sd(comparisons$estimate+comparisons$std.error))+
        ggplot2::theme(text = ggplot2::element_text(size=14))+
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
      suppressMessages(ggplot2::ggsave(paste0(home_dir, "results figures/", exposure, "-", outcome, " cutoff= ", cutoff, ".jpeg"), plot=ggplot2::last_plot()))
    }

    cat("\n")
    cat("See the 'results figures/' folder for graphical representations of results","\n")
  })

}
