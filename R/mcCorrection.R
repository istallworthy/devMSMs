
#' Multiple comparison correction
#' Adjusts p-values for multiple comparisons by exposure-outcome pairing
#' @param object msm object that contains all relevant user inputs
#' @param history_comparisons output from compareHistories
#' @importFrom stats p.adjust
#' @importFrom stargazer stargazer
#' @seealso [stats::p.adjust()] [comparisonHistories()]
#' @return significant_comparisons
#' @examples mcCorrection(object, history_comparisons)
mcCorrection <- function(object, history_comparisons){

  home_dir=object$home_dir
  method=object$mc_method
  exposure_epochs=object$exposure_epochs
  weights_percentile_cutoff=object$weights_percentile_cutoff
  dose_level=object$dose_level


  significant_comparisons=list()

  for (x in 1:length(history_comparisons)){
    exp_out=names(history_comparisons)[x]

    #get p-values from saved models
    uncorr_p_vals=lapply(history_comparisons[[exp_out]], function(x) x$`Pr(>F)`[2])
    uncorr_p_vals <- data.frame(contrast = rep(names(uncorr_p_vals), sapply(uncorr_p_vals, length)),
                     p_val = unlist(uncorr_p_vals))
    comparisons=uncorr_p_vals[order(uncorr_p_vals$p_val),]

    #apply user-specified correction
    comparisons$p_vals_corr=stats::p.adjust(comparisons$p_val, method=method)

    sig_comparisons=comparisons[comparisons$p_vals_corr<0.05,]

    significant_comparisons[[exp_out]] <- sig_comparisons


    folder_label=ifelse(as.numeric(sapply(strsplit(exp_out, "_"), "[", 3))==weights_percentile_cutoff, "original/", "sensitivity checks/")

    # browser()
    cat(paste0("USER ALERT: please inspect the followig uncorrected and corrected contrast comparisons for ", exp_out, " :"), "\n")
    cat(knitr::kable(comparisons, format='pipe'), sep="\n")
    cat("\n")

    #save out table for all contrasts with old and corrected p-values
    sink(paste0(home_dir, "msms/linear hypothesis testing/", folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_table.doc", sep=""))
    stargazer::stargazer(comparisons, type="html", digits=2, column.labels = colnames(comparisons),summary=FALSE, rownames = FALSE, header=FALSE,
                         out=paste0(home_dir, "msms/linear hypothesis testing/", folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_table.doc", sep=""))
    sink()

    # print(stargazer::stargazer(comparisons, type="html", digits=2, column.labels = colnames(comparisons),summary=FALSE, rownames = FALSE, header=FALSE,
    #                            out=paste0(home_dir, "msms/linear hypothesis testing/", folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_table.doc", sep=""))
    # )

    #make fancy table
   histories={}
   ref=data.frame(matrix(nrow=1, ncol=nrow(exposure_epochs)+1))
   ref[1,1:(ncol(ref)-1)]=as.data.frame(t(rep("l", nrow(exposure_epochs))))
   colnames(ref)=c(exposure_epochs$epochs, "p-value of Comparison to Reference After MC Correction")
   #adding in levels (h v l) at each exposure epoch
   for (e in 1:nrow(exposure_epochs)){
       histories[[as.symbol(exposure_epochs$epochs[e])]] <-sapply(strsplit(sapply(strsplit(comparisons$contrast, "vs."), "[", 2), "-"), "[", e)
       histories[[as.symbol(exposure_epochs$epochs[e])]]=gsub(" ", "", histories[[as.symbol(exposure_epochs$epochs[e])]])
   }
   histories$`p-value of Comparison to Reference After MC Correction`=comparisons$p_vals_corr
   histories=rbind(ref, histories)
   exposure_history=data.frame(`Exposure History`=LETTERS[seq(from=1, to=nrow(comparisons)+1)])
   exposure_dose=data.frame(`Exposure Dose`= rowSums(histories[1:nrow(exposure_epochs)]==dose_level))
   tab=cbind(exposure_history, exposure_dose, histories)

   sink(paste0(home_dir, "msms/linear hypothesis testing/", folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_presentation_table.doc", sep=""))
   stargazer::stargazer(tab, type="html", digits=2, column.labels = colnames(tab),summary=FALSE, rownames = FALSE, header=FALSE,
                        out=paste0(home_dir, "msms/linear hypothesis testing/", folder_label, sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_presentation_table.doc", sep=""))
   sink()
  }


  cat("See 'msms/linear hypothesis testing/original/' folder for tables of likelihood ratio tests","\n")
  cat("See 'msms/linear hypothesis testing/sensitivity checks/' folder for sensitivity checks","\n")

  return(significant_comparisons)
}
