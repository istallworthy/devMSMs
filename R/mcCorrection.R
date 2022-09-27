
#' Multiple comparison correction
#' Adjusts p-values for multiple comparisons by exposure-outcome pairing
#' @param history_comparisons output from compareHistories
#' @param method abbreviation for multiple comparison correction method (see p.adjust)
#' @importFrom stats p.adjust
#' @importFrom stargazer stargazer
#' @seealso [stats::p.adjust()] [comparisonHistories()]
#' @return significant_comparisons
#' @examples mcCorrection(history_comparisons, method="BH")
mcCorrection <- function(history_comparisons, method="BH"){

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

    stargazer::stargazer(comparisons, type="html", digits=2, column.labels = colnames(comparisons),summary=FALSE, rownames = FALSE,
                         out=paste0(home_dir, "msms/", sapply(strsplit(exp_out, "-"), "[",1), "_", sapply(strsplit(exp_out, "-"), "[",2), "_lht_table.doc", sep=""))



  }

  print("See 'msms' folder for tables of likelihood ratio tests")
  return(significant_comparisons)
}
