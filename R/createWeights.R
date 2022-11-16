#' Creates balancing weights for potential confounding covariates in relation to exposures at each time point
#'
#' Creates balancing weights using the CBPS package (See CBPS documentation for more detail: https://cran.r-project.org/web/packages/CBPS/CBPS.pdf) for more
#' returns a list of weights_models for each exposure-outcome pair
#' @param object msm object that contains all relevant user inputs
#' @param wide_long_datasets from formatForWeights
#' @param forms from createForms
#' @param ATT CBPS parameter 1= average treatment effect on the treated, 0=average treatment effect
#' @param iterations CBPS parameter for maximum number of number of iterations for optimization
#' @param standardize CBPS parameter for whether weights should be normalized to sum to 1; F returns Horvitz-Thompson weights
#' @param method CBPS parameter "exact"= fit a model that contains only covariate balancing conditions, "over"= over-identified model
#' @param twostep CBPS parameter for two-step estimator, runs fater than continuous updating, set to F for continuous updating from Imai & Ratkovic 2014
#' @param sample.weights CBPS parameter for survey sampling weights for obs, NULL defaults to sampling weight of 1 for each observation
#' @param baseline.formula CBPS parameter for iCBPS, only works with binary treatments
#' @param diff.formula CBPS parameter for iCBPS, only works with binary treatments
#' @return weights_models
#' @export
#' @importFrom CBPS CBPS
#' @importFrom CBPS balance
#' @importFrom data.table setDT
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggsave
#' @seealso [CBPS::CBPS()], [formatForWeights()], [createForms()]
#' @examples createWeights(object, wide_long_datasets,forms,read_in_from_file="no", ATT=0, iterations=1000, standardize=FALSE, method="exact", twostep=TRUE, sample.weights=NULL, baseline.forumula=NULL, diff.formula=NULL)
#'
createWeights <-function(object, wide_long_datasets, forms, read_in_from_file="no", ATT=0, iterations=1000, standardize=FALSE, method="exact", twostep=TRUE, sample.weights=NULL, baseline.formula=NULL, diff.formula=NULL){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  exposure_time_pts=object$exposure_time_pts
  m=object$m
  balance_thresh=object$balance_thresh

  if (read_in_from_file=="yes"){ #reads in saved weights instead of re-creating
    weights_models=readRDS(paste0(home_dir, "original weights/weights_models.rds"))

    return(weights_models)

  }else{ #otherwise, calculate weights

    weights_models=list()
    #Cycles through imputed datasets
    for (k in 1:m){

      MSMDATA=wide_long_datasets[[paste("imp", k, "_widelong", sep="")]]

      for (z in seq(length(outcomes))){
        outcome=outcomes[z]

        #cycles through exposures
        for (y in 1:length(exposures)){
          exposure=exposures[y]

          #Cycles through all time points
          for (x in 1:length(exposure_time_pts)){
            time=exposure_time_pts[x]

            #references newly ordered time variable (1:n(exposure_time_pts))
            new_time=as.numeric(c(1:length(exposure_time_pts))[x])

            form=forms[[paste("form_", exposure, "-", outcome, "_", time, sep="")]]

            #Order by time
            MSMDATA <- MSMDATA[order(MSMDATA$WAVE),]

            #select only data at the appropriate time point
            MSMDATA_temp<- MSMDATA %>%
              dplyr::filter(WAVE==as.numeric(new_time))

            #Use the CBPS continuous version of this function to generate weights
            fit <- CBPS::CBPS(form, data=MSMDATA_temp, ATT=ATT, iterations = iterations,
                              standardize = standardize, method = method, twostep = twostep, #twostep=false does continuous updating
                              sample.weights = sample.weights, baseline.formula = NULL,
                              diff.formula = NULL)

            weights_models[[paste("fit_", k, "_", exposure, "-", outcome, "_", time, sep="")]] <-fit

            #saves out fit objects made by CBPS so user can re-load at later time
            # saveRDS(fit, file = paste(paste0(home_dir, "weight fits/fit_imp=", k, "_exp=", exposure, "_t=", time, ".rds", sep="")))

            # print(paste0("Weights model objects for exposure ", exposure, " at time ", time, " for imputation ", k,  " have been saved in the 'weight fits' folder"))

            #get weights from output
            weights=as.data.frame(cbind(MSMDATA_temp[,colnames(MSMDATA_temp)==ID], fit$weights))

            #Save weights
            # write.csv(x=as.data.frame(fit$weights), file=paste(home_dir, "weights/weights_exp=", exposure, "-", outcome, "_t=", time, "_imp=", k, ".csv", sep=""))
            #Save weights merged with ID variable
            write.csv(x=as.data.frame(weights), file=paste(home_dir, "original weights/values/weights_id_exp=", exposure, "-", outcome, "_t=", time,"_imp=",k, ".csv", sep=""))

            cat(paste0("Weights for exposure ", exposure, "-", outcome," at time ", time, " for imputation ", k,  " have now been saved into the 'original weights/values/' folder"),"\n")

            # #Writes image of histogram of weights to assess heavy tails
            ggplot2::ggplot(data=as.data.frame(fit$weight), ggplot2::aes(x = fit$weight)) +
              ggplot2::geom_histogram(color = 'black', bins = 15)
            ggplot2::ggsave(paste("Hist_exp=", exposure, "-", outcome, "_t=", time, "_imp=", k, ".png", sep=""), path=paste0(home_dir, "original weights/histograms"), height=8, width=14)

            cat(paste0("A weights histogram for exposure ", exposure,  "-", outcome," at time ", time, " for imputation ", k,  " has now been saved in the 'original weights/histograms/' folder --likely has heavy tails"),"\n")

            #Writes image to balance check
            jpeg(filename=paste(home_dir, "balance/post-balance correlation plots/Balance_exp=", exposure, "-", outcome, "_t=", time, "_imp=",k, ".jpg", sep=""), width=480,height=480)
            plot(fit, covars = NULL, silent = FALSE, boxplot = TRUE)
            dev.off()

            cat(paste0("USER ALERT: Balancing figures for exposure ", exposure,  "-", outcome," at time ", time, " for imputation ", k,  " have now been saved into the 'balance/post-balance correlation plots/' folder for future inspection", "\n"))

          }
        }
      }
    }

    saveRDS(weights_models, file = paste(paste0(home_dir, "original weights/weights_models.rds", sep="")))
    cat("Weights models have been saved as an .rds object in the 'original weights' folder","\n")

    return(weights_models)
  }
}




