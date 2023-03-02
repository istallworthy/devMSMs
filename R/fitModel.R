
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

fitModel <- function(object, data_for_model_with_weights_cutoff, unbalanced_covariates_for_models=NULL, model="m4"){

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

  #adapted from from Francis Huang (https://francish.net/post/multiple-imputation-in-r-with-regression-output/)
  library(psych)
  library(texreg)
  extract.svymi <- function(z){
    #require(norm)
    m2 <- length(z)
    ns <- nobs(z[[1]]) #first dataset
    betas <- lapply(z, FUN = function(x){coef(x)} )
    sess <- lapply(z, FUN = function(x){vcov(x)}) #the vcov
    #ses <- lapply(sess, FUN = function(x){sqrt(diag(x))}) #the se, not the complete vcov

    dm <- lapply(z, FUN = function(x){model.matrix(x)})
    r.1 <- sapply(1:m2, function(row) cor(dm[[row]] %*% betas[[row]], z[[row]]$y))
    fisher.r2z <- function(r) {0.5 * (log(1+r) - log(1-r))}
    zz <- mean(fisher.r2z(r.1))
    r2 <- psych::fisherz2r(zz)^2 #in psych
    R2 <- mean(r2)

    test1 <- as.data.frame(summary(miceadds::pool_mi(betas, sess)))

    test1$term <- row.names(test1)
    test2 <- test1[,c(8,1,2,4)] #4 is the pvalue

    row.names(test2) <- c()
    names(test2) <- c('term','estimate','std.error', 'p.value')

    # tr <- texreg::createTexreg(
    #   coef.names = test2$term,
    #   coef = test2$estimate,
    #   se = test2$std.error,
    #   pvalues = test2$p.value,
    #   gof.names = c("R2", "Nobs"),
    #   gof = c(R2, ns),
    #   gof.decimal = c(T,F))

    # tr=as.data.frame(tr)
    #
    tr2 <- list(
      coef.names = test2$term,
      coef = test2$estimate,
      se = test2$std.error,
      pvalues = test2$p.value,
      gof.names = c("R2", "Nobs"),
      gof = c(R2, ns),
      gof.decimal = c(T,F)
    )

  }


  all_cutoffs=c(weights_percentile_cutoff, weights_percentile_cutoffs_sensitivity)


  if (length(outcome_time_pt)>1){
    stop('This function is designed only for single time point outcome')}

  names(data_for_model_with_weights_cutoff) <- gsub(".", "_", names(data_for_model_with_weights_cutoff), fixed=TRUE)

  # data_for_model_with_weights_cutoff=as.data.frame(data_for_model_with_weights_cutoff)


  #model dictionary:
  #m0 baseline model with only main effects for each exposure epoch
  #m1 full covariate model adding in all unbalanced covariates (if any)
  #m2 covariate model including only significant covariates (retaining all epoch main effects)
  #m3 full interaction model adding in all possible epoch interactions
  #m4 final model including all epoch main effects and only the significant covariates and interactions


  #cycles through each exposure to fit the marginal structural model relating exposure histories to each outcome
  all_models=list()

  cat("USER ALERT: please inspect the following series of models for each exposure-outcome pair (using weights cutoff at the original user-specified percentile value):", "\n")
  cat("\n")
  cat("\n")

  for (c in 1:length(all_cutoffs)){
    cutoff=all_cutoffs[c]

    # for (y in 1:length(outcome)){
    #   outcome=outcome[y]


    #cycles through each outcome
    models=list()

    #
    #     for (x in 1:length(exposure)){
    #       exposure=exposure[x]
    imp_models=list()


    #lists out exposure-epoch combos
    exp_epochs= apply(expand.grid(exposure, as.character(exposure_epochs[,1])), 1, paste, sep="", collapse="_")



    data_all_imp=data_for_model_with_weights_cutoff[which(grepl(paste(exposure, "-", outcome, "_", cutoff, sep=""), names(data_for_model_with_weights_cutoff)))]
    # data_all_imp=unname(data_all_imp)
    # data_all_imp=as.mira(data_all_imp)
    # data_all_imp=lapply( 1:m , function( n ) complete( data_for_model_with_weights_cutoff , action = n ) )

    #renaming to make it easier for svydesign
    # weights_name=colnames(tidyr::complete(data_for_model_with_weights_cutoff , 1))[grepl(paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff"), colnames(complete(data_for_model_with_weights_cutoff, 1)))]
    weights_name=colnames( data_all_imp[[1]])[grepl(paste0(exposure, "-", outcome, "_", cutoff, "_weight_cutoff"), colnames(data_all_imp[[1]]))]
    # weights_name="InRatioCor.EF_avg_perc_0.98_weight_cutoff"
    data_all_imp=lapply(data_all_imp, function(x) {names(x)[names(x)==weights_name] <- "weights";x})
    data_all_imp=lapply(data_all_imp, function(x) {names(x)[names(x)==ID] <- "ID";x})


    a=complete(imp_data_w_t, action="long", include = TRUE)
    colnames(a)[colnames(a)==weights_name]="weights"
    colnames(a)[colnames(a)==ID]="ID"
    # a=a%>%dplyr::filter(.imp>0)
    c=mice::as.mids(a, .imp = ".imp")

    library(mitools)
    #creates initial survey design type specifying ID, data, and weights to use
    s=svydesign(id=~ID,
                # data=data_for_model_with_weights_cutoff,
                data=mitools::imputationList(data_all_imp), #adds list of imputation data
                # data=data_all_imp, #adds list of imputation data
                weights=~weights)
    # data_for_model_with_weights_cutoff[,colnames(data_for_model_with_weights_cutoff)=="cutoff_weights"])
    # paste0(exposure,"-", outcome, "_", paste(unlist(strsplit(as.character(cutoff), "\\.")), sep="", collapse="_"), "_weight_cutoff")])

    #creates forms
    #creates form for baseline model including averaged exposure at each epoch as a predictor
    f0=paste(paste0(outcome, ".", outcome_time_pt), "~", paste0(exp_epochs, sep="", collapse=" + "))

    cutoff_label=ifelse(cutoff==weights_percentile_cutoff, paste0("original weight cutoff value (", cutoff, ")"), paste0("sensitivity test weight cutoff value (", cutoff, ")"))



    #fits initial baseline model to all imps
    # m0=svyglm(noquote(f0), design=s)
    m0=with(s, svyglm(noquote(f0), design=s))
    # a=mitools::MIcombine(m0)

    #see vignette: https://cran.r-project.org/web/packages/mitml/vignettes/Analysis.html
    m0.pool=mitml::testEstimates(m0) #df.com=10)


    # test=mitools::MIcombine(m0)
    # #combine results across imps
    # m0.comb=summary(mitools::MIcombine(m0))
    a<-MIcombine(m0)
    # pt( abs(coef(a)/SE(a)),df=a$df,lower.tail=FALSE)*2 #p-values

    #OR
    sink(tempfile())
    m0.pool <- extract.svymi(m0)
    sink()
    # m0.res=texreg::screenreg( m0.pool, digits = 3, custom.model.names = "using svyglm")


    models[["m0"]]<-m0


    #adding in covariates that did not fully balance when creating the weights
    covariate_list= unbalanced_covariates_for_models[[paste0(exposure, "-", outcome)]]
    # grepl(paste(factor_covariates, collapse="|"), paste0(unlist(strsplit(noquote(covariate_list), "\\+"))))

    if (covariate_list[1]==""){
      if(cutoff==weights_percentile_cutoff){
        cat(paste0("There are no unbalanced covariates to include in the model of effects of ", exposure, " on ", outcome),"\n")
        cat("\n")
      }
      #f1 and f2 are same as f0
      f1=paste(f0) #no covars
      f2=paste(f0) #no sig covars


    }else{

      # f1=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", covariate_list)
      f1=paste(f0, "+", covariate_list)

      # browser()

      # m1=svyglm(noquote(f1), design=s)
      m1=with(s, svyglm(noquote(f1), design=s))

      sink(tempfile())
      m1.pool <- extract.svymi(m1)
      sink()
      # m1.res=texreg::screenreg(m1.pool, digits = 3, custom.model.names = "using svyglm")

      # if(cutoff==weights_percentile_cutoff){
      #   cat(paste0("Covariate model results for effects of ", exposure, " on ", outcome,  " using ", cutoff_label),"\n")
      #   print(summary(m1))
      #   cat("\n")
      # }
      models[["m1"]]<-m1

      #determining which covariates are significant
      # sig_covars=as.data.frame(summary(m1)$coefficients)
      # sig_covars=as.data.frame(coefs=m1.pool$coef)

      # sig_covars=sig_covars[sig_covars$`Pr(>|t|)`<0.05,]
      sig_covars=m1.pool$coef.names[m1.pool$pvalues<0.05]
      # sig_covars=rownames(sig_covars)
      sig_covars=sig_covars[!grepl(c("(Intercept)"), sig_covars)]
      sig_covars=sig_covars[!grepl(c(exposure), sig_covars)]
      sig_covars=c(as.character(sig_covars))

      #in the case of no significant covariates
      if (length(sig_covars)!=0){
        f2=paste(f0, "+", paste(sig_covars, sep="", collapse=" + "))
        if(cutoff==weights_percentile_cutoff){
          cat(paste0("The only significant covariate(s) for this model are: ", paste(sig_covars, sep="", collapse=" , ")),"\n")
          cat("\n")
        }

      }else{
        f2=f0
      }

      # m2=svyglm(noquote(f2), design=s)
      m2=with(s, svyglm(noquote(f2), design=s))

      sink(tempfile())
      m2.pool <- extract.svymi(m2)
      sink()

      # if(cutoff==weights_percentile_cutoff){
      #   cat(paste0("Final covariate model results for effects of ", exposure, " on ", outcome," using ", cutoff_label),"\n")
      #   print(summary(m2))
      #   cat("\n")
      # }
      models[["m2"]]<-m2
    }


    #getting interactions
    interactions=paste(apply(combn(exp_epochs,2), 2, paste, sep="", collapse=":"), sep="", collapse=" + ")
    f3=paste(f2, "+", paste(interactions, sep="", collapse=" + "))

    #testing for interactions
    # f3=paste(paste0(outcome, "_", outcome_time_pts), "~", paste0(exp_epochs, sep="", collapse=" + "), "+", sig_covars, "+", interactions)
    # m3=svyglm(noquote(f3), design=s)
    m3=with(s, svyglm(noquote(f3), design=s))

    sink(tempfile())
    m3.pool <- extract.svymi(m3)
    sink()


    models[["m3"]]<-m3

    #determining which interactions are signficant
    # sig_ints=as.data.frame(summary(m3)$coefficients)
    sig_covars=m3.pool$coef.names[m3.pool$pvalues<0.05]
    # sig_ints=sig_ints[sig_ints$`Pr(>|t|)`<0.05,]
    sig_ints=rownames(sig_covars)
    sig_ints=sig_ints[grepl(c(":"), sig_ints)] #makes sure to only get interactions

    #creating FINAL model
    if (length(sig_ints)==0){
      f4=f2 #equal to model with any covariates
    }else{
      if(cutoff==weights_percentile_cutoff){
        cat(paste0("The only significant interactions(s) for this model are: ", paste(sig_ints, sep="", collapse=" , ")), "\n")
        cat("\n")
      }

      f4=paste(f2, "+", paste(sig_ints, sep="", collapse=" + "))

    }

    # cat("The final model is:", "/n")
    # m4=svyglm(noquote(f4), design=s)
    m4=with(s, svyglm(noquote(f4), design=s))

    sink(tempfile())
    m4.pool <-invisible(extract.svymi(m4))
    sink()

    models[["m4"]]<-m4

    cat(paste0("The marginal model selected for ", exposure, "-", outcome, "_cutoff_", cutoff, " is summarized here:"), "\n")
    extract.svymi(get(model))


    all_models[[paste0(exposure, "-", outcome, "_cutoff_", cutoff)]]<-models

    models=NULL
    m0=NULL
    m1=NULL
    m2=NULL
    m3=NULL
    m4=NULL

    # }

    # }

    cat("\n")

    file_label=ifelse(cutoff==weights_percentile_cutoff, "original", "sensitivity checks")

    saveRDS(all_models, file = paste(paste0(home_dir, "msms/", file_label, "/all_exposure-outcome_models_weight_cutoff_", cutoff, ".rds", sep="")))


  }
  # cat(paste0("All models with ", file_label, " cutoff weights have been saved as an .rds object in the 'msms/", file_label, "/' folder"),"\n")
  cat("\n")
  cat("\n")

  # write.csv(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.csv"))
  saveRDS(data_for_model_with_weights_cutoff, paste0(home_dir, "msms/data_for_msms.rds"))

  cat("USER ALERT: Sensitivity check models using two alternate weight truncation values have been saved in the 'msms/sensitivity checks/' folder")
  cat("\n")
  cat("A new data file has been saved as a .csv file in the in the 'msms' folder","\n")


  return(all_models)

}
