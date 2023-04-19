#' Creates balancing weights for potential confounding covariates in relation to exposure at each time point
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
createWeights <-function(object, wide_long_datasets, short_forms, read_in_from_file="no"){

  ID=object$ID
  home_dir=object$home_dir
  exposure=object$exposure
  outcome=object$outcome
  exposure_time_pts=object$exposure_time_pts
  m=object$m
  balance_thresh=object$balance_thresh
  time_varying_covariates=object$time_varying_variables
  weights_method=object$weights_method

  library(cobalt)

  library(dbarts)

  # browser()
  data_for_model_with_weights=list()

  if (read_in_from_file=="yes"){ #reads in saved weights instead of re-creating
    # weights_models=readRDS(paste0(home_dir, "original weights/weights_models.rds"))
    dat_w <-readRDS(paste(paste0(home_dir, "original weights/dat_w.rds", sep="")))

    return(dat_w)

  }else{ #otherwise, calculate weights


    if(!weights_method %in% c("cbps", "npcbps", "bart", "ps")){
      stop("Please enter a weighting method in the msmObject from this list: cbps, npcbps, bart, ps")
    }

    weights_models=list()
    # all_weights=data.frame(expand.grid(wide_long_datasets[[paste("imp", 1, "_widelong", sep="")]][,ID], 1:m))
    # colnames(all_weights)<-c(ID,".imp")


    # all_weights_exp=data.frame()
    # name=paste0(exposure, "-", outcome)

    #list of forms for each time point
    form=short_forms[which(grepl(paste("form_", exposure, "-", outcome, sep=""), names(short_forms)))]
    form=unname(form)

    cat("For the ", weights_method, " weighting method:", "\n")

    # for (k in 1:m){
    dat_w<-lapply(1:length(wide_long_datasets), function(i){
      k=i
      d=wide_long_datasets[[i]]
      d=as.data.frame(d)

      # browser()

      # MSMDATA=wide_long_datasets[[paste("imp", k, "_widelong", sep="")]]

      #determining exposure type
      exposure_type=ifelse(length(unique(d[,colnames(d)[colnames(d)==paste0(exposure,".", exposure_time_pts[1])]]))<3, "binary", "continuous")

      set.seed(1234)
      fit <- WeightIt::weightitMSM(form, data=d,
                                   method=weights_method,
                                   stabilize=T, #recommended stabilized weights
                                   density="dt_2", #using t-distribution w/ 2 DF for bart
                                   use.kernel=T,
                                   over=F) #similar to exact=T in cbsp package

      d$weights <- fit$weights

      bal=cobalt::bal.tab(fit, stats = c("corr"), continuous="std", binary="std",
                          un=T,thresholds = c(corr =balance_thresh), which.time=.none)

      cat(paste0("Preliminary summary balance statistics for imputation ", k, " using the ", weights_method, " weighting method, ", bal$Balanced.correlations$count[1], " covariates were balanced and ",
                 bal$Balanced.correlations$count[2], " covariates remain imbalanced."), "\n")
      # cat('\n')

      cobalt::bal.tab(fit, stats = c("corr"), continuous="std", binary="std",
                      un=T,thresholds = c(corr =balance_thresh), which.time=.all)

      weights_models[[k]] <-fit

      cat(paste0("For ", k, "_", exposure, "-", outcome, " , the median weight value is ", round(median(fit$weights),2) ,
                 " (SD= ", round(sd(fit$weights),2), "; range= ", round(min(fit$weights),2), " - ", round(max(fit$weights),2), ")"), "\n")
      cat('\n')

      #get weights from output and cbind to data
      # weights=as.data.frame(cbind(MSMDATA_temp[,colnames(MSMDATA_temp)==ID], fit$weights))
      # w=data.frame(imp=k,
      #              ID=MSMDATA[,ID],
      #              name=fit$weights)
      # colnames(w)=c(".imp",ID, paste0(name, "_weights"))
      #
      # all_weights_exp=rbind(all_weights_exp,w)
      #
      #
      # weights=merge(MSMDATA, w%>%dplyr::select(c(ID,paste0(name, "_weights"))), by=ID)
      # data_for_model_with_weights[[paste("fit_", k, "_", exposure, "-", outcome, sep="")]]<-weights

      # browser()

      #Save weights merged with ID variable
      write.csv(x=as.data.frame(d), file=paste(home_dir, "original weights/values/weights_method=",weights_method, "_",k, ".csv", sep=""))

      # cat("\n")

      # #Writes image of histogram of weights to assess heavy tails
      ggplot2::ggplot(data=as.data.frame(fit$weight), ggplot2::aes(x = fit$weight)) +
        ggplot2::geom_histogram(color = 'black', bins = 15)
      suppressMessages(ggplot2::ggsave(paste("Hist_weights_method=",weights_method, "_", k, ".png", sep=""), path=paste0(home_dir, "original weights/histograms"), height=8, width=14))

      # cat("\n")

      # cat(paste0("USER ALERT: Balancing figures for ", exposure,  "-", outcome," for imputation ", k,  " have now been saved into the 'balance/plots/' folder for future inspection", "\n"))
      # cat("\n")

      fit=NULL
      d
    }) #ends k loop

    cat(paste0("Weights for each imputation have now been saved into the 'original weights/values/' folder"),"\n")
    cat(paste0("Weights histogram for each imputation have now been saved in the 'original weights/histograms/' folder --likely has heavy tails"),"\n")


    #gathering all weights
    # all_weights<-merge(all_weights, all_weights_exp, by=c(ID, ".imp"), all.x=T)

    #   } #ends exposure loop
    #
    # } #ends outcome loop

    # }

    # #adds weights to mids object
    # # colnames(all_weights)=c(ID, "weights")
    # a=complete(imputed_datasets, action="long", include = TRUE)
    # # a$weights=NA
    # b=merge(a, all_weights, by=c(ID, ".imp"), all.x=T)
    # imp_data_w=mice::as.mids(b, .imp = ".imp")
    #
    #
    #
    #
    # saveRDS(data_for_model_with_weights, file = paste(paste0(home_dir, "original weights/data_for_model_with_weights.rds", sep="")))
    saveRDS(dat_w, file = paste(paste0(home_dir, "original weights/weights_method=", weights_method, "_dat_w.rds", sep="")))

    saveRDS(weights_models, file = paste(paste0(home_dir, "original weights/weights_method=", weights_method, "_weights_models.rds", sep="")))
    cat("\n")
    cat("Weights models have been saved as an .rds object in the 'original weights' folder","\n")

    # return(weights_models)
    return(dat_w)

  }


}



