#' Creates balancing weights for potential confounding covariates in relation to exposure at each time point
#'
#' returns a list of weights_models f
#' @param object msm object that contains all relevant user inputs
#' @param wide_long_datasets from formatForWeights
#' @param short_forms from createShortForms
#' @param read_in_from_file optional binary indicator of whether weights should be read in from local file
#' @return weights_models
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 ggsave
#' @importFrom WeightIt weightitMSM
#' @importFrom cobalt bal.tab
#' @seealso [CBPS::CBPS()], [formatForWeights()], [createShortForms()]
#' @examples createWeights(object, wide_long_datasets, short_forms, read_in_from_file="no")
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

  #getting form type from input
  .call <- match.call()
  form_type<-.call[[4]]
  form_name<-as.character(form_type)

  library(cobalt) #still not sure how to call this without loading lib

  data_for_model_with_weights=list()

  if (read_in_from_file=="yes"){ #reads in saved weights instead of re-creating
    tryCatch({dat_w <- readRDS(paste(paste0(home_dir, "original weights/", exposure, "-", outcome, "_", form_name, "_", weights_method, "_dat_w.rds", sep="")))},
             error=function(x){cat("These weights have not previously been saved locally. Please re-run with read_in_from_file='no'")})
    cat("Reading in balancing weights from local folder.")
    cat("\n")
    return(dat_w)

  }else{ #otherwise, calculate weights


    if(!weights_method %in% c("cbps", "npcbps", "bart", "ps", "entropy")){ #will remove entropy as has not been validated for longit exp
      stop("Please enter a weighting method in the msmObject from this list: cbps, npcbps, bart, ps")
    }

    cat(paste0("Creating longitudinal balancing weights using the ", weights_method, " weighting method and the ", form_name, " formulas"), "\n")

    weights_models=list()

    #list of forms for each time point
    form=short_forms[which(grepl(paste("form_", exposure, "-", outcome, sep=""), names(short_forms)))]
    form=unname(form)

    #not sure if we should output the old bal stats from cobalt or nah?
    cat("You will see some preliminary balance statistics (output from the cobalt::bal_tab() function. These are only preliminary. Please use the assesBalance function for final balance statistics.", "\n")

    #cycling through imputed datasets
    dat_w<-lapply(1:length(wide_long_datasets), function(i){
      k=i
      d=wide_long_datasets[[i]]
      d=as.data.frame(d)

      #determining exposure type
      exposure_type=ifelse(length(unique(d[,colnames(d)[colnames(d)==paste0(exposure,".", exposure_time_pts[1])]]))<3, "binary", "continuous")

      set.seed(1234)
      #creating weights
      fit <- WeightIt::weightitMSM(form, data=d,
                                   method=weights_method,
                                   stabilize=T, #recommended stabilized weights
                                   density="dt_2", #using t-distribution w/ 2 DF for bart
                                   use.kernel=T,
                                   weightit.force = TRUE, #to allow entropy (temp)
                                   over=F) #similar to exact=T in cbsp package
      d$weights <- fit$weights

      #prelim bal stats
      if (exposure_type=="continuous"){
        bal=cobalt::bal.tab(fit, stats = c("corr"), continuous="std", binary="std",
                            un=T,thresholds = c(corr =balance_thresh), which.time=.none)
        cat(paste0("Preliminary summary balance statistics for imputation ", k, " using the ", weights_method, " weighting method and the ", form_name, " formulas,",
                   bal$Balanced.correlations$count[1], " covariates were balanced and ",
                   bal$Balanced.correlations$count[2], " covariates remain imbalanced."), "\n")
        cobalt::bal.tab(fit, stats = c("corr"), continuous="std", binary="std",
                        un=T,thresholds = c(corr =balance_thresh), which.time=.all)
      }
      if(exposure_type=="binary"){ #need to test
        bal=cobalt::bal.tab(fit, continuous="std", binary="std",
                            un=T,thresholds = c(m =balance_thresh), which.time=.none)
        cat(paste0("Preliminary summary balance statistics for imputation ", k, " using the ", weights_method, " weighting method and the ", form_name, " formulas, ",
                   bal$Balanced.correlations$count[1], " covariates were balanced and ",
                   bal$Balanced.correlations$count[2], " covariates remain imbalanced."), "\n")
        cobalt::bal.tab(fit, continuous="std", binary="std",
                        un=T,thresholds = c(m =balance_thresh), which.time=.all)
      }


      weights_models[[k]] <-fit

      cat(paste0("For ", k, "_", exposure, "-", outcome, ", the median weight value is ", round(median(fit$weights),2) ,
                 " (SD= ", round(sd(fit$weights),2), "; range= ", round(min(fit$weights),2), "-", round(max(fit$weights),2), ")."), "\n")
      cat('\n')


      #Save weights merged with ID variable
      write.csv(x=as.data.frame(d), file=paste(home_dir, "original weights/values/", exposure, "-", outcome,"_", form_name, "_", weights_method, "_",k, ".csv", sep=""))

      # #Writes image of histogram of weights to assess heavy tails
      ggplot2::ggplot(data=as.data.frame(fit$weight), ggplot2::aes(x = fit$weight)) +
        ggplot2::geom_histogram(color = 'black', bins = 15)
      suppressMessages(ggplot2::ggsave(paste("Hist_", exposure,"-", outcome, "_", form_name,"_",weights_method, "_", k, ".png", sep=""),
                                       path=paste0(home_dir, "original weights/histograms/"), height=8, width=14))

      fit=NULL
      d
    }) #ends k loop

    cat(paste0("Weights for each imputation have now been saved into the 'original weights/values/' folder"),"\n")
    cat(paste0("Weights histograms for each imputation have now been saved in the 'original weights/histograms/' folder --likely has heavy tails"),"\n")


    # saveRDS(data_for_model_with_weights, file = paste(paste0(home_dir, "original weights/data_for_model_with_weights.rds", sep="")))
    saveRDS(dat_w, file = paste(paste0(home_dir, "original weights/", exposure, "-", outcome,"_", form_name, "_", weights_method, "_dat_w.rds", sep="")))

    saveRDS(weights_models, file = paste(paste0(home_dir, "original weights/",  exposure, "-", outcome,"_", form_name, "_", weights_method, "_weights_models.rds", sep="")))
    cat("\n")
    cat("Weights models have been saved as an .rds object in the 'original weights' folder","\n")

    return(dat_w)
  }
}



