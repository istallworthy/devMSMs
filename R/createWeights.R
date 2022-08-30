#See CBPS documentation for more detail: https://cran.r-project.org/web/packages/CBPS/CBPS.pdf

#' Title
#'
#' @return
#' @export
#' @importFrom CBPS CBPS
#' @importFrom CBPS balance
#' @importFrom data.tabe setDT
#'
#' @param exposures list of exposures
#' @param time_pts list of time points
#' @param balance_thresh correlation value above which covariates are not considered balanced with respect to expsure
#' @param ATT CBPS parameter 1= average treatment effect on the treated, 0=average treatment effect
#' @param iteraetions CBPS parameter for maximum number of number of iterations for optimization
#' @param standardize CBPS parameter for whether weights should be normalized to sum to 1; F returns Horvitz-Thompson weights
#' @param method CBPS parameter "exact"= fit a model that contains only covariate balancing conditions, "over"= over-identified model
#' @param twostep CBPS parameter for two-step estimator, runs fater than continuous updating, set to F for continuous updating from Imai & Ratkovic 2014
#' @param sample.weights CBPS parameter for survey sampling weights for obs, NULL defaults to sampling weight of 1 for each observation
#' @param baseline.formula CBPS parameter for iCBPS, only works with binary treatments
#' @param diff.formula CBPS parameter for iCBPS, only works with binary treatments


#' @examples
#'
createWeights <-function(exposures, time_pts, balance_thresh=0.12, ATT=0, iterations=1000, standardize=FALSE, method="exact", twostep=TRUE, sample.weights=NULL, baseline.forumula=NULL, diff.formula=NULL){

  ####Fit CBPS for all treatments and time points using forms created above for each imputed dataset, each tx, and each time point. Creates variables fit_k_treat for each imputed dataset. Also saves out balance checking image, histogram, and weights when PLOT_CODE is set to "yes" to 'weights' folder.

  #Cycles through imputed datasets
  for (k in 1:5){

    MSMDATA=get(paste("imp", k, "_widelong", sep="")) #reads in widelong imputed data

    #cycles through exposures
    for (y in 1:length(exposures)){
      exposure=exposures[y]

      #Cycles through all time points
      for (x in 1:length(time_pts)){
        time=time_pts[x]

        #references newly ordered time variable (1:n(time_pts))
        new_time=as.numeric(c(1:length(time_pts))[x])

        form=get(paste("form_", exposure, "_", time, sep="")) #pulls appropriate form

        #Order by time
        MSMDATA <- MSMDATA[order(MSMDATA$WAVE),]

        #select only data at the appropriate time point
        MSMDATA_temp<- MSMDATA %>%
          dplyr::filter(WAVE==as.numeric(new_time))

        #Use the CBPS continuous version of this function to generate weights
        fit <- CBPS(form, data=MSMDATA_temp, ATT=ATT, iterations = iterations,
                    standardize = standardize, method = method, twostep = twostep, #twostep=false does continuous updating
                    sample.weights = sample.weights, baseline.formula = NULL,
                    diff.formula = NULL)

        return(assign(paste("fit_", k, "_", exposure, "_", time, sep=""), fit))

        #saves out fit objects made by CBPS so user can re-load at later time
        saveRDS(get(paste("fit_", k, "_", treat, "_", time, sep="")), file = paste(paste0(home_dir, "weight fits/fit_", k, "_", treat, "_", time, ".rds", sep="")))
        print("Weights objects have been saved in the 'weight fits' folder so you can access them later")

        #get weights from output
        weights=as.data.frame(cbind(MSMDATA_temp$s_id, fit$weights))

        #Writes image to balance check in imp folder
        jpeg(filename=paste(homedir, "balance/Balance_", exposure, "_", time, "_",k, ".jpg", sep=""), width=480,height=480)
        plot(fit, covars = NULL, silent = FALSE, boxplot = TRUE)
        dev.off()

        print("Please check for successful balancing using the images in the 'balance' folder")


        #added by IS 3.21.20 for more formal balance checking (prints out which variables are unbalanced)

        #assessing balance and determining unbalanced covariates (i.e., those still correlated above 0.12)

        balance=balance(fit)
        balance=as.data.frame(balance$balanced)
        balance=as.data.frame(balance)
        balance=setDT(balance, keep.rownames = TRUE)[]
        unbalanced_vars=balance[abs(balance$`Pearson Correlation`)>balance_thresh,] #0.12 corr default threshold

        #returns covariates that remain associated with exposure
        return(assign(paste("unbalanced_vars", exposure, "_", time, "_", k, sep=""), unbalanced_vars))

        #writes out list of unbalanced variables
        write.csv(x=as.data.frame(unbalanced_vars), file=paste(home_dir, "balance/unbalanced_vars_", exposure, "_", time,"_",k, ".csv", sep=""))
        unbalanced_vars=NULL


        # #Writes image of histogram of weights to assess heavy tails
        ggplot(data=as.data.frame(fit$weight), aes(x = fit$weight)) +
          geom_histogram(color = 'black', bins = 15)
        ggsave(paste("Hist_", exposure, "_", time, "_", k, ".png", sep=""), path=paste0(home_dir, "weights/", height=8, width=14))

        #Save weights
        write.csv(x=as.data.frame(fit$weights), file=paste(home_dir, "weights/weights_", exposure, "_", time, "_", k, ".csv", sep=""))

        #Save weights merged with ID variable
        write.csv(x=as.data.frame(wights), file=paste(home_dir, "weights/weights_id_", exposure, "_", time,"_",k, ".csv", sep=""))

      }
    }
  }
}




