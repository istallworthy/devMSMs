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
createWeights <-function(object, wide_long_datasets, forms, read_in_from_file="no"){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  outcomes=object$outcomes
  exposure_time_pts=object$exposure_time_pts
  m=object$m
  balance_thresh=object$balance_thresh
  time_varying_covariates=object$time_varying_variables
  weights_method=object$weights_method

  # browser()
  data_for_model_with_weights=list()

  if (read_in_from_file=="yes"){ #reads in saved weights instead of re-creating
    # weights_models=readRDS(paste0(home_dir, "original weights/weights_models.rds"))
    imp_data_w <-readRDS(paste(paste0(home_dir, "original weights/imp_data_w.rds", sep="")))

    return(imp_data_w)

  }else{ #otherwise, calculate weights

    weights_models=list()
    all_weights=data.frame(expand.grid(wide_long_datasets[[paste("imp", 1, "_widelong", sep="")]][,ID], 1:m))
    colnames(all_weights)<-c(ID,".imp")
    #Cycles through imputed datasets

    for (z in seq(length(outcomes))){
      outcome=outcomes[z]

      #cycles through exposures
      for (y in 1:length(exposures)){
        exposure=exposures[y]

        # all_weights_exp=data.frame(expand.grid(wide_long_datasets[[paste("imp", 1, "_widelong", sep="")]][,ID], 1:m))
        # colnames(all_weights_exp)<-c(ID,".imp")
        # all_weights_exp=all_weights_exp%>%dplyr::mutate(!!paste0(exposure, "-", outcome, "_weights"):=NA)
        #
        all_weights_exp=data.frame()
        name=paste0(exposure, "-", outcome)

        for (k in 1:m){

          MSMDATA=wide_long_datasets[[paste("imp", k, "_widelong", sep="")]]
          MSMDATA$PmMrSt2=as.numeric(MSMDATA$PmMrSt2)
          #Order by time
          # MSMDATA <- MSMDATA[order(MSMDATA$WAVE),]

          #determining exposure type
          exposure_type=ifelse(length(unique(MSMDATA[,colnames(MSMDATA)[colnames(MSMDATA)==paste0(exposure,".", exposure_time_pts[1])]]))<3, "binary", "continuous")



          #Cycles through all time points
          # for (x in 1:length(exposure_time_pts)){
          #   time=exposure_time_pts[x]

          #references newly ordered time variable (1:n(exposure_time_pts))
          # new_time=as.numeric(c(1:length(exposure_time_pts))[x])

          # form=forms[[paste("form_", exposure, "-", outcome, "-", time, sep="")]]

          #gathers all forms for exposure and outcome (across all time points)
          form=forms[which(grepl(paste("form_", exposure, "-", outcome, sep=""), names(forms)))]
          # form=form[1:2]
          # form=unname(form)
          # f=as.formula(InRatioCor.15 ~BioDadInHH2 + KFASTScr + PEmply2 + PmAge2 + PmBlac2 +
          #                PmEd2 + PmMrSt2 + ScEd2 + state + TcBlac2 + TcSex2+
          #                IBRReac.6 + Cotin_rep1ngmL.6+ pcx_pos.6
          #                )


          #FIT 1 CBPS
          #Use the CBPS continuous version of this function to generate weights
          # fit.2 <- CBPS::CBPS(as.formula(forms[[1]]), data=MSMDATA_temp, ATT=ATT, iterations = iterations,
          # fit.2b <- CBPS::CBPS(as.formula(form[[2]]), data=MSMDATA, ATT=0, iterations = 1000,
          #                   standardize = TRUE, method = "exact", twostep = F, #twostep=false does continuous updating
          #                   )
          # CBPS::balance(fit.2b)
          # plot(fit.2b, boxplot = T)
          # cobalt::bal.tab(fit.2b, stats = c("corr"), continuous="std", binary="std",
          #         un=T,thresholds = c(cor =0.1))
          # cobalt::love.plot(fit.2b,
          #                   var.order="unadjusted",
          #                   line=TRUE,
          #                   stat=c("corr"),
          #                   continuous="std",
          #                   binary="std",
          #                   thresholds=c(0.1), #save out this one plotting std m diffs for all covariates
          #                   drop.distance=TRUE,
          #                   colors = c("red", "blue"),
          #                   title= "Covariate Balance Summary",
          #                   shapes = c("triangle filled", "circle filled"),
          #                   sample.names = c("Unweighted", "Weighted"))
          #
          # # form=list(unlist(form))
          #
          # form=as.formula(form)

          # #Order by time
          # MSMDATA <- MSMDATA[order(MSMDATA$WAVE),]

          #select only data at the appropriate time point
          # MSMDATA_temp<- MSMDATA %>%
          #   dplyr::filter(WAVE==as.numeric(new_time))

          # browser()

          library(cobalt)

          # # #using weightit
          # fit.1 <- WeightIt::weightit(as.formula(form[[1]]), data=MSMDATA, method="cbps", stabilize=TRUE, int=T, over=F)
          # cobalt::bal.tab(fit.1, stats = c("corr"), continuous="std", binary="std",
          #         un=T,thresholds = c(corr =0.1))
          # cobalt::love.plot(fit.1,
          #                   var.order="unadjusted",
          #                   line=TRUE,
          #                   stat=c("corr"),
          #                   continuous="std",
          #                   binary="std",
          #                   thresholds=c(0.1), #save out this one plotting std m diffs for all covariates
          #                   drop.distance=TRUE,
          #                   colors = c("red", "blue"),
          #                   title= "Covariate Balance Summary",
          #                   shapes = c("triangle filled", "circle filled"),
          #                   sample.names = c("Unweighted", "Weighted"))
          #
          # cor(fit.1$weights, fit.2$weights)


          # library(SuperLearner)
          library(dbarts)

          fit.3b <- WeightIt::weightitMSM(form, data=MSMDATA, method=weights_method, stabilize=T, int=T, density="dt_2", use.kernel=T, over=F)

          cobalt::bal.tab(fit.3b, stats = c("corr"), continuous="std", binary="std",
                          un=T,thresholds = c(corr =0.1), which.time=1)
          cobalt::love.plot(fit.3b,
                            var.order="unadjusted",
                            line=TRUE,
                            stat=c("corr"),
                            continuous="std",
                            binary="std",
                            thresholds=c(0.1), #save out this one plotting std m diffs for all covariates
                            drop.distance=TRUE,
                            colors = c("red", "blue"),
                            title= "Covariate Balance Summary",
                            shapes = c("triangle filled", "circle filled"),
                            sample.names = c("Unweighted", "Weighted"))

          # cor(fit.1$weights, fit.2$weights)


          # forms_times=list(as.formula(forms[[1]]))
          # forms_times=list(as.formula(forms[[1]]), as.formula(forms[[2]]))

          # fit.2 <- WeightIt::weightit(as.formula(forms[[1]]), data=MSMDATA_temp, method="npcbps", stabilize=TRUE, int=TRUE, over=F)
          # bal.tab(fit, stats = c("corr", "v"), thresholds = c(m = .05))

          # f=forms[[1]]
          # #setting balance stats based on exposure type
          # bal_vars=ifelse(exposure_type=="binary", c("m", "ks"), c("c","m"))

          # fit_1_time <- WeightIt::weightit(f, data=MSMDATA, method="cbps", stabilize=TRUE, int=TRUE,
          #                                  over=F)
          # bal.tab(fit_1_time, stats = bal_vars, thresholds = c(cor=0.1), un=T,continuous="std", binary="std",
          #         which.time = .none)
          #
          # bal.plot(fit_1_time, "KFASTScr", which = "both")
          # love.plot(fit_1_time, stats = c("c", "ks"),
          #           thresholds = c(cor = .1),
          #           abs = TRUE, wrap = 20,
          #           var.order = "unadjusted", line = TRUE)


          #temp
          # # MSMDATA$InRatioCor.15=NULL
          # MSMDATA$InRatioCor.24=NULL
          # MSMDATA$InRatioCor.35=NULL
          # MSMDATA$InRatioCor.58=NULL

          #
          # f2=as.formula(InRatioCor.15 ~ BioDadInHH2 + KFASTScr + PEmply2 + PmAge2 + PmBlac2 +
          #                 PmEd2 + PmMrSt2 + ScEd2 + state + TcBlac2 + TcSex2)
          # form=as.list(c(f,f2))
          #
          # fit_all_time <- WeightIt::weightitMSM(form, data=MSMDATA, method="cbps", stabilize=T, int=T,
          #                                       over=F)
          #
          # fit_all_time
          # summary(fit_all_time)
          #
          # bal.tab(fit_all_time, stats = c("c"), thresholds = c(cor=0.1), un=T,continuous="std", binary="std",
          #         which.time = .all)
          #
          # love.plot(fit_all_time, stats = c("c"),
          #           thresholds = c(cor = .1),
          #           abs = TRUE, wrap = 20,
          #           var.order = "unadjusted", line = TRUE)




          # fit_all_time
          # summary(fit_all_time)
          # bal.tab(fit_all_time, stats = c("correlation", "ks"), thresholds = c(m = .05),
          #         which.time = .all)

          # #conducts weighting on all imputed datasets; no need to loop
          # fit <- MatchThem::weightthem(form, data=imputed, method="cbps", approach="within", over=F)
          # bal.tab(fit, stats = c("m", "v"), thresholds = c(m = .05))







          fit=fit.3b
          # weights_models[[paste("fit_", k, "_", exposure, "-", outcome, "-", time, sep="")]] <-fit
          weights_models[[paste("fit_", k, "_", exposure, "-", outcome, sep="")]] <-fit


          cat(paste0("For ", "fit_", k, "_", exposure, "-", outcome, " , the median weight is ", median(fit$weights) ,
                     " (SD= ",sd(fit$weights), "; range= ", min(fit$weights), " - ", max(fit$weights), ")"), "\n")




          #get weights from output and cbind to data
          # weights=as.data.frame(cbind(MSMDATA_temp[,colnames(MSMDATA_temp)==ID], fit$weights))
          w=data.frame(imp=k,
                       ID=MSMDATA[,ID],
                       name=fit$weights)
          colnames(w)=c(".imp",ID, paste0(name, "_weights"))

          all_weights_exp=rbind(all_weights_exp,w)




          weights=merge(MSMDATA, w%>%dplyr::select(c(ID,paste0(name, "_weights"))), by=ID)
          data_for_model_with_weights[[paste("fit_", k, "_", exposure, "-", outcome, sep="")]]<-weights


          #Save weights merged with ID variable
          write.csv(x=as.data.frame(weights), file=paste(home_dir, "original weights/values/weights_id_exp=", exposure, "-", outcome, "_imp=",k, ".csv", sep=""))

          cat(paste0("Weights for ", exposure, "-", outcome," for imputation ", k,  " have now been saved into the 'original weights/values/' folder"),"\n")
          # cat("\n")

          # #Writes image of histogram of weights to assess heavy tails
          ggplot2::ggplot(data=as.data.frame(fit$weight), ggplot2::aes(x = fit$weight)) +
            ggplot2::geom_histogram(color = 'black', bins = 15)
          suppressMessages(ggplot2::ggsave(paste("Hist_exp=", exposure, "-", outcome, "_imp=", k, ".png", sep=""), path=paste0(home_dir, "original weights/histograms"), height=8, width=14))

          cat(paste0("A weights histogram for ", exposure,  "-", outcome, " for imputation ", k,  " has now been saved in the 'original weights/histograms/' folder --likely has heavy tails"),"\n")
          # cat("\n")

          # #Writes image to balance check
          # jpeg(filename=paste(home_dir, "balance/plots/Balance_exp=", exposure, "-", outcome, "_imp=",k, ".jpg", sep=""), width=480,height=480)
          # plot(fit, covars = NULL, silent = FALSE, boxplot = TRUE)
          # dev.off()

          cat(paste0("USER ALERT: Balancing figures for ", exposure,  "-", outcome," for imputation ", k,  " have now been saved into the 'balance/plots/' folder for future inspection", "\n"))
          cat("\n")

          fit=NULL

        } #ends k loop

        #gathering all weights
        all_weights<-merge(all_weights, all_weights_exp, by=c(ID, ".imp"), all.x=T)

      } #ends exposure loop

    } #ends outcome loop

  }

  #adds weights to mids object
  # colnames(all_weights)=c(ID, "weights")
  a=complete(imputed_datasets, action="long", include = TRUE)
  # a$weights=NA
  b=merge(a, all_weights, by=c(ID, ".imp"), all.x=T)
  imp_data_w=mice::as.mids(b, .imp = ".imp")




  saveRDS(data_for_model_with_weights, file = paste(paste0(home_dir, "original weights/data_for_model_with_weights.rds", sep="")))


  saveRDS(imp_data_w, file = paste(paste0(home_dir, "original weights/imp_data_w.rds", sep="")))

  saveRDS(weights_models, file = paste(paste0(home_dir, "original weights/weights_models.rds", sep="")))
  cat("\n")
  cat("Weights models have been saved as an .rds object in the 'original weights' folder","\n")

  # return(weights_models)
  return(imp_data_w)

}





