
#' Fit weighted model
#' This code fits a weighted marginal structural model to examine the effects of different exposure histories on outcome
#' @param object msm object that contains all relevant user inputs
#' @param data_for_model_with_weights_cutoff dataset with truncated weights see truncateWeights
#' @param unbalanced_covariates_for_models unbalanced covariates see assessBalance
#' @importFrom survey svydesign
#' @importFrom survey svyglm
#' @return marginal structural models
#' @export
#' @seealso [truncateWeights()], [asesssBalance()]
#' @examples fitModel(object, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models)

fitModel <- function(object, data_for_model_with_weights_cutoff, balance_stats_final, model="m3"){

  home_dir=object$home_dir
  ID=object$ID
  exposure=object$exposure
  exposure_epochs=object$exposure_epochs
  outcome=object$outcome
  outcome_time_pt=object$outcome_time_pt
  factor_covariates=object$factor_covariates
  weights_percentile_cutoff=object$weights_percentile_cutoff
  weights_percentile_cutoffs_sensitivity=c(as.numeric(unlist(strsplit(object$weights_percentile_cutoffs_sensitivity, " "))))
  m=object$m
  balance_thresh=object$balance_thresh

  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)

  unbalanced_covars=as.data.frame(rowMeans(do.call(cbind, lapply(balance_stats_final, "[", "std_bal_stats"))))
  unbalanced_covars=data.frame(exposure=exposure,
                               exp_time=balance_stats_final[[1]]$exp_time,
                               covar_time=balance_stats_final[[1]]$covar_time,
                               covariate=balance_stats_final[[1]]$covariate,
                               avg_bal=unname(unbalanced_covars))%>%
    dplyr::mutate(balanced_avg=ifelse(abs(avg_bal)<balance_thresh,1,0))%>%
    dplyr::filter(balanced_avg==0)%>%
    dplyr::select(covariate)
  #makes wide for modeling
  unbalanced_covars=sapply(strsplit(unbalanced_covars$covariate, "\\."), "[",1)

  #renaming any factor covariates
  # unbalanced_covars$covariate[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]
  # unbalanced_covars$covariate[unbalanced_covars$covariate %in% factor_covariates] <-sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(unbalanced_covars$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]

  if (!model %in% c("m0", "m1", "m2", "m3")){
    stop('Please provide a valid model "m" from 0-3 (e.g., "m1")')
  }

  if (length(outcome_time_pt)>1){
    stop('This function is designed only for single time point outcome')}

  names(data_for_model_with_weights_cutoff) <- gsub(".", "_", names(data_for_model_with_weights_cutoff), fixed=TRUE)



  #cycles through each exposure to fit the marginal structural model relating exposure histories to each outcome
  all_models=list()


  #creates matrix of k and m for sequencing through
  # cut_m=expand.grid(all_cutoffs, 1:m)
  # cut_m=as.data.frame(cut_m)
  # colnames(cut_m)=c("cutoff", "k")
  #
  library(MatchThem)

  #creating lists of models for each truncation cutoff value
  all_models <-lapply(1:length(all_cutoffs), function(i){
    cutoff=all_cutoffs[i]
    # k=cut_m$k[i]

    #cycles through each outcome
    models=list()

    imp_models=list()

    #lists out exposure-epoch combos
    exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")

    # data_all_imp=data_for_model_with_weights_cutoff[which(grepl(paste(exposure, "-", outcome, "_", cutoff, sep=""), names(data_for_model_with_weights_cutoff)))]

    #getting mids object of all imputed data with truncated weights
    data_all_imp= readRDS(paste0(home_dir, "final weights/values/imp_data_w_t.rds"))

    weights_name=colnames( data_all_imp[[1]])[grepl(paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff"), colnames(data_all_imp[[1]]))]

    #renaming for ease
    a=complete(data_all_imp, action="long", include = TRUE)
    colnames(a)[colnames(a) %in% paste0(exposure, ".", outcome, "_", cutoff, "_weight_cutoff")]<-"weights"
    b=mice::as.mids(a, .imp = ".imp")

    # names(data_all_imp$imp)[names(data_all_imp$imp)==weights_name]  <- "weights" #renaming for ease
    # data_all_imp=lapply(data_all_imp, function(x) {names(x$imp)[names(x$imp)==weights_name] <- "weights";names(x$imp)})
    # data_all_imp=lapply(data_all_imp, function(x) {names(x)[names(x)==ID] <- "ID";x})

    #creating large mild file with all imputed datasets
    Wd2 <- MatchThem::complete(b, "all")

    # #renaming to make it easier for svydesign
    # # weights_name=colnames( data_all_imp[[1]])[grepl(paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff"), colnames(data_all_imp[[1]]))]
    # Wd2=lapply( Wd2, function(x) {names(x)[names(x)==weights_name] <- "weights";x[]})
    # Wd2 <- MatchThem::complete(mice::as.mids(Wd2), "all")

    # data_all_imp=lapply( Wd2, function(x) {names(x)[names(x)==ID] <- "ID";x})

   a= lapply(1:length(dat_w_t), function(y){
     browser()
     d=dat_w_t[[y]]
     lapply(d, function(x){
       d=d[[x]]
    survey::svydesign(~ID,#~1, #
                  # data=data_for_model_with_weights_cutoff,
                  data=d, #adds list of imputation data
                  # data=data, #adds list of imputation data
                  weights=~weights)

    })
   })


    s=svydesign(~1, #id=~ID,#~1, #
                # data=data_for_model_with_weights_cutoff,
                data=mitools::imputationList(Wd2), #adds list of imputation data
                # data=data, #adds list of imputation data
                weights=~weights)

    #
    #     s=svydesign(~1, #id=~ID,#~1, #
    #                 # data=data_for_model_with_weights_cutoff,
    #                 data=mitools::imputationList(Wd2), #adds list of imputation data
    #                 # data=data, #adds list of imputation data
    #                 weights=~`HOMEETA1-EF_avg_perc_0.95_weight_cutoff`)

    #https://r-survey.r-forge.r-project.org/survey/svymi.html
    f0=paste("outcome", "~", paste0(exp_epochs, sep="", collapse=" + "))
    m0=with(s, survey::svyglm(as.formula(f0))) #list of model fitted to all imputed datasets

    # #renaming to make it easier for svydesign
    # weights_name=colnames( data_all_imp[[1]])[grepl(paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff"), colnames(data_all_imp[[1]]))]
    # data_all_imp=lapply(data_all_imp, function(x) {names(x)[names(x)==weights_name] <- "weights";x})
    # data_all_imp=lapply(data_all_imp, function(x) {names(x)[names(x)==ID] <- "ID";x})

    # data=data_all_imp[[k]]

    #creates initial survey design type specifying ID, data, and weights to use
    # s=svydesign(id=~ID,
    #             data=data, #adds list of imputation data
    #             weights=~weights)


    cutoff_label=ifelse(cutoff==weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))

    model_ref=data.frame()

    #fits initial baseline model to all imps
    # f0=paste(paste0(outcome, ".", outcome_time_pt), "~", paste0(exp_epochs, sep="", collapse=" + "))
    # m0=svyglm(as.formula(f0), design=s)

    #baseline model (main effects) only
    if(model=="m0"){
      models[["m0"]]<-m0 #save model
    }else{

      #listing any imbalanced covariates
      covariate_list= paste(as.character(unlist(unbalanced_covars)), sep="", collapse=" + ")

      ##gathers imbalanced covariates if they exist
      #baseline + sig covar model OR baseline + sig covar + int model
      if(model=="m1" | model =="m3"){
        #if there are no imbalanced covariates
        if(covariate_list[1]==""){
          stop("You have selected a covariate model but there are no imbalanced covariates for one or more imputed datasets. Please choose another model.")
        }else{ #there are imbalanced covariates
          cat("The following imbalanced covariates are included in the covariate model: ", "\n")
          cat(covariate_list)
          f1.temp=paste(f0, "+", covariate_list) #baseline + covariate model
          # m1.temp=svyglm(noquote(f1.temp), design=s)
          # m1.temp=with(s, survey::svyglm(as.formula(f1.temp)))

          # #determining which covariates are significant
          # sig_covars=as.data.frame(summary(m1.temp)$coefficients)
          # sig_covars=sig_covars[sig_covars$`Pr(>|t|)`<0.05,]
          # sig_covars=rownames(sig_covars)
          # sig_covars=sig_covars[!grepl(c("(Intercept)"), sig_covars)]
          # sig_covars=sig_covars[!grepl(c(exposure), sig_covars)]
          # sig_covars=c(as.character(sig_covars))

          #in the case of significant covariates
          # if (length(sig_covars)!=0){
          # f1=paste(f0, "+", paste(sig_covars, sep="", collapse=" + ")) #baseline + sig covar model
          f1=f1.temp
          # cat(paste0("The only significant covariate(s) for this model are: ", paste(sig_covars, sep="", collapse=" , ")),"\n")
          # cat("\n")
          # m1=svyglm(noquote(f1), design=s)
          m1=with(s, survey::svyglm(as.formula(f1)))

          #baseline + imbalanced covars
          if(model=="m1"){
            models[["m1"]]<-m1
          }
          # }else{
          #   stop("You have selected a covariate model but there are no imbalanced covariates that are statistically significant for one or more imputed datasets. Please choose another model.")
          # }
        }
      }



      ##Gathers interactions if they exist
      #baseline + interaction OR baseline + covars + interactions
      if (model=="m2" | model=="m3"){
        #getting interactions
        interactions=paste(apply(combn(exp_epochs,2), 2, paste, sep="", collapse=":"), sep="", collapse=" + ")
        #baseline + interactions
        if (model=="m2"){
          f2.temp=paste(f0, "+", paste(interactions, sep="", collapse=" + "))
          # m2.temp=svyglm(noquote(f2.temp), design=s)

          # #determining which interactions are signficant
          # sig_ints=as.data.frame(summary(m2.temp)$coefficients)
          # sig_ints=sig_ints[sig_ints$`Pr(>|t|)`<0.05,]
          # sig_ints=rownames(sig_covars)
          # sig_ints=sig_ints[grepl(c(":"), sig_ints)] #makes sure to only get interactions

          # if (length(sig_ints)==0){
          #   stop("You have selected an interaction (m2) model but there are no interactions are statistically significant for. Please choose another model.")
          # # }else{ #there are sig ints
          # cat(paste0("The only significant interactions(s) for this model are: ", paste(sig_ints, sep="", collapse=" , ")), "\n")
          # cat("\n")
          # f2=paste(f0, "+", paste(sig_ints, sep="", collapse=" + "))
          f2=f2.temp
          # m2=svyglm(noquote(f2), design=s)
          m2=with(s, survey::svyglm(as.formula(f2)))
          #baseline + interactions
          models[["m2"]]<-m2
          # }
        }

        #baseline + covars + interactions
        if (model=="m3"){
          f3.temp=paste(f2, "+", paste(interactions, sep="", collapse=" + "))
          # m3.temp=svyglm(noquote(f3.temp), design=s)
          #
          # #determining which interactions are signficant
          # sig_ints=as.data.frame(summary(m3.temp)$coefficients)
          # sig_ints=sig_ints[sig_ints$`Pr(>|t|)`<0.05,]
          # sig_ints=rownames(sig_covars)
          # sig_ints=sig_ints[grepl(c(":"), sig_ints)] #makes sure to only get interactions

          # if (length(sig_ints)==0){
          #   stop("You have selected an interaction (m3) model but no interactions are statistically significant. Please choose another model.")
          # }else{ #there are sig ints
          # cat(paste0("The only significant interactions(s) for this model are: ", paste(sig_ints, sep="", collapse=" , ")), "\n")
          # cat("\n")
          # f3=paste(f2, "+", paste(sig_ints, sep="", collapse=" + "))
          f3=f3.temp
          # m3=svyglm(noquote(f3), design=s)
          m3=with(s, survey::svyglm(as.formula(f3)))
          #baseline + covars+ interactions
          models[["m3"]]<-m3
          # }
        }
      }

    }


    cat(paste0("The marginal model selected and run for each imputed dataset using the weights truncation cutoff ", cutoff, " is summarized here:"), "\n")
    # print(summary(models[[model]]))
    print(lapply(models[[model]], summary))

    # all_models[[paste0("fit", k, "_", exposure, "-", outcome, "_cutoff_", cutoff)]]<-models

    # models=NULL
    m0=NULL
    m1=NULL
    m2=NULL
    m3=NULL


    cat("\n")

    file_label=ifelse(cutoff==weights_percentile_cutoff, "original", "sensitivity checks")
    saveRDS(all_models, file = paste(paste0(home_dir, "msms/", file_label, "/all_exposure-outcome_models_weight_cutoff_", cutoff, ".rds", sep="")))

    models


  })

  # browser()
  names(all_models)<-paste0(exposure, "-", outcome, "_cutoff_", all_cutoffs)

  # cat(paste0("All models with ", file_label, " cutoff weights have been saved as an .rds object in the 'msms/", file_label, "/' folder"),"\n")


  # write.csv(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.csv"))
  saveRDS(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.rds"))

  cat("USER ALERT: Sensitivity check models using two alternate weight truncation values have been saved in the 'msms/sensitivity checks/' folder")
  cat("\n")
  cat("A new data file has been saved as a .csv file in the in the 'msms' folder","\n")


  return(all_models)

}
