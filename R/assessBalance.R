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
assessBalance <- function (object, forms, data_for_model_with_weights, weighted){

  home_dir=object$home_dir
  m=object$m
  exposure=object$exposure
  outcome=object$outcome
  balance_thresh=object$balance_thresh

  # forms=short_forms

  #gathering balance stats for all imputed datasets
  bal_stats= lapply(seq(m), function(i){
    calcBalStats(object, data_for_model_with_weights, forms, exposure, outcome, k=i, weighted)
  })


  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(bal_stats, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=bal_stats[[1]]$exp_time,
                               covar_time=bal_stats[[1]]$covar_time,
                               covariate=bal_stats[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0))



  #save out correlations/std mean differences
  sink(paste0(home_dir, "balance/comparison values/all_post-balance_associations.html"))
  stargazer::stargazer(unbalanced_covars,type="html", digits=2, column.labels = colnames(unbalanced_covars),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/comparison values/all_post-balance_assocations.html"))
  sink()

  cat(paste0("USER ALERT: Check 'balance/comparison values/' folder for a table of all post-balance correlations or standardized mean differences and averaged across imputed datasets"),"\n")
  # cat("\n")

  unbalanced_covars= unbalanced_covars%>%
    dplyr::filter(balanced_avg==0)

  cat(paste0("Averaging across all imputed datasets, the following ", nrow(unbalanced_covars) ," covariates remain imbalanced: "), "\n")
  # print(unbalanced_covars)
  cat(knitr::kable(unbalanced_covars, caption="Imbalanced Covariates", format='pipe'),  sep="\n")


  return(bal_stats)

}






