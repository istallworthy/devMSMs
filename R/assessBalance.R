#' Assess success of covariate balancing
#'
#' #Assesses how well balance was achieved for each of the covariates/potential confounds in relation to each of the exposures,
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
#' @importFrom ggplot2 scale_x_continuous, theme
#' @importFrom knitr kable
#' @importFrom cobalt love.plot
#' @importFrom dplyr group_by
#' @examples assessBalance(object, weights_models=list(), just_made_weights="no")
assessBalance <- function (object, data, weights_models){

  # knitr::opts_chunk$set(include= FALSE, warning = FALSE, message = FALSE)
  require(cobalt) #currently does not work to call using 2 colons; Noah is aware and working on it

  home_dir=object$home_dir
  m=object$m
  exposures=object$exposures
  outcomes=object$outcomes
  exposure_time_pts=object$exposure_time_pts
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates

  #assessing balance and determining unbalanced covariates (i.e., those still correlated above 0.12)
  unbalanced_variables=list()
  all_post_balance_corrs=list()

  #iterating through saved weights models ("fits") to gather all post balance correlations
  for (f in 1:length(weights_models)){

    balance=CBPS::balance(weights_models[f][[1]])

    #pre and post balance comparisons
    temp=do.call(rbind, lapply(balance, as.data.frame))
    temp <- tibble::rownames_to_column(temp, "rn")
    temp=as.data.frame(temp)
    status=sapply(strsplit(temp[,1], "\\."), "[", 1)
    covariate=sapply(strsplit(temp[,1], "\\."), "[", 2)
    temp=as.data.frame(cbind(status, covariate, temp[,2:ncol(temp)]))

    # browser()

    #comparing pre and post correlations for continuous exposures
    if (sum(grepl("Correlation", colnames(balance$balanced)))>0){
      temp=reshape(temp, idvar="covariate", timevar="status", direction="wide")
      colnames(temp)=c("Covariate", "Balanced Corr", "Unweighted Corr")
      temp$`Balanced Corr`=as.numeric(temp$`Balanced Corr`)
      temp$`Unweighted Corr`=as.numeric(temp$`Unweighted Corr`)
      cat(paste0("For ", names(weights_models[f]), ", the average correlation with exposure across all confounders was ", abs(mean(temp$`Unweighted Corr`)) ," prior to balancing and ", abs(mean(temp$`Balanced Corr`)) ," following balancing"), "\n")
      sum=as.data.frame(summary(as.data.frame(temp)))
      sum=sum[7:nrow(sum),2:ncol(sum)]
      temp=rbind(temp, cbind(data.frame("Covariate"=sapply(strsplit(sum[sum$Var2=="Balanced Corr", 2], "\\:"), "[", 1),
                                        "Balanced Corr"=as.numeric(sapply(strsplit(sum[sum$Var2=="Balanced Corr", 2], "\\:"), "[", 2)),
                                        "Unweighted Corr" =as.numeric(sapply(strsplit(sum[sum$Var2=="Unweighted Corr", 2], "\\:"), "[", 2)),check.names=FALSE)))
      temp$Difference=abs(temp[,2]-temp[,3])

      sink(paste0(home_dir, "balance/comparison values/", names(weights_models[f]), "_pre_vs_post-balance_correlations.html"))
      stargazer::stargazer(temp,type="html", digits=2, column.labels = colnames(temp),
                           summary=FALSE, rownames = TRUE, header=FALSE,
                           out=paste0(home_dir, "balance/comparison values/", names(weights_models[f]), "_pre_vs_post-balance_correlations.html"))
      sink()

      # cat("\n")
      cat("USER ALERT: see the 'balance/comparison values/' folder for a table comparing covariate relations to exposure pre and bost balancing", "\n")


      #comparing pre and post weighting standardized mean differences for binary exposure
    }else if (sum(grepl("mean",colnames(balance$balanced)))>0){
      temp=do.call(rbind, lapply(balance, as.data.frame))
      temp=temp[,3:ncol(temp)]
      temp <- tibble::rownames_to_column(temp, "rn")
      temp=as.data.frame(temp)
      status=sapply(strsplit(temp[,1], "\\."), "[", 1)
      covariate=sapply(strsplit(temp[,1], "\\."), "[", 2)
      temp=as.data.frame(cbind(status, covariate, temp[,2:ncol(temp)]))
      temp$std.mean.diff=as.numeric(unlist(temp[4]))-as.numeric(unlist(temp[3]))
      temp=as.data.frame(temp)
      sum=as.data.frame(temp)%>%dplyr::select(status, contains("mean"))%>%
        dplyr::group_by(factor(status))%>%dplyr::summarise(across(where(is.numeric), dplyr::funs(mean, sd, min, max)))
        # dplyr::group_by(factor(status))%>%dplyr::summarise_all(dplyr::funs(mean, sd, min, max))
      sum=sum[,colnames(sum)[grepl("status_", colnames(sum))==F]]
      difference=sum[1,]-sum[2,]
      sum=rbind(sum, difference)
      cat(paste0("For ", names(weights_models[f]), " the standardized group mean difference across all confounders was ", abs(sum$std.mean.diff_mean[2]) ," prior to balancing and ", abs(sum$std.mean.diff_mean[1]) ," following balancing"), "\n")
      sum=t(sum)
      colnames(sum)= c("Balanced", "Original", "Difference")

      sink(paste0(home_dir, "balance/comparison values/", names(weights_models[f]), "_pre_vs_post-balance_group_diffs.html"))
      stargazer::stargazer(sum[2:nrow(sum),],type="html", digits=2, column.labels = c("Stat", "Balanced", "Original", "Difference") ,
                           summary=FALSE, rownames = TRUE, header=TRUE,
                           out=paste0(home_dir, "balance/comparison values/", names(weights_models[f]), "_pre_vs_post-balance_group_diffs.html"))
      sink()

      # cat("\n")
      cat("USER ALERT: see the 'balance/comparison values/' folder for a table comparing covariate relations to exposure pre and bost balancing", "\n")
      set.cobalt.options(binary="std")

      }


    balance=as.data.frame(balance$balanced) #finds balanced results
    balance=as.data.frame(balance)
    balance=data.table::setDT(balance, keep.rownames = TRUE)[]

    #new balancing assessment from cobalt
    # remotes::install_github("ngreifer/cobalt")

    #ISSUE: not sure how to save this out well
    # bal.tab=cobalt::bal.tab(weights_models[f][[1]])
    # # bal.tab=data.table::as.data.table(bal.tab[1:2])
    # write.csv(bal.tab[1:2],paste0(home_dir, "balance/post-balance correlation values/", names(weights_models[f]), "all_std_mean_diffs.html") )
    # stargazer::stargazer(bal.tab,type="html", digits=2 ,
    #                      summary=FALSE, rownames = FALSE, header=TRUE,
    #                      out=paste0(home_dir, "balance/post-balance correlation values/", names(weights_models[f]), "all_std_mean_diffs.html"))

    # browser()
    #makes a plot showing std mean differences for all covariates before and after balancing
    suppressWarnings(cobalt::love.plot(weights_models[f][[1]],
                      var.order="unadjusted",
                      line=TRUE,
                      thresholds=c(balance_thresh), #save out this one plotting std m diffs for all covariates
                      drop.distance=TRUE,
                      colors = c("red", "blue"),
                      title= "Covariate Balance Summary",
                      shapes = c("triangle filled", "circle filled"),
                      sample.names = c("Unweighted", "Weighted"))+
    # position = c(.8, .2)) +
      ggplot2::theme(legend.box.background = ggplot2::element_rect(),
                     legend.position = "bottom",
                     legend.box.margin = ggplot2::margin(1, 1, 1, 1),
                     axis.text=ggplot2::element_text(size=7))+
      ggplot2::scale_x_continuous(breaks=c(seq(-1, 1, by=0.1))))
    suppressWarnings(ggplot2::ggsave(paste0(home_dir, "balance/plots/", names(weights_models[f]), "_summary_balance_plot.jpeg")))

    # cat("\n")
    cat(paste("USER ALERT: see the 'balance/plots/' folder for a summary balance plot for ", names(weights_models[f])), "\n")
    cat("\n")

    # cobalt::bal.plot(fit, "RGotFodS", which="both", mirror=T)
    # cobalt::bal.plot(fit, "WIND.6", which="both")


    all_post_balance_corrs[[paste("unbalanced_vars_", names(weights_models)[f], sep="")]] <-balance

  }




  #averages across imputed datasets to determine final variables that remain correlated with exposure at r>0.12
  sig_unbalanced_averaged_k=list()
  unbalanced_covariates_for_models=list()
  all_post_balance_correlations=list()

  all_corrs_df={}


  for (z in seq(length(outcomes))){

    for (x in 1:length(exposures)){

      significant_corrs_remaining={}

      #determining exposure type
      exposure_type=ifelse(length(na.omit(unique(data[,colnames(data)[colnames(data)==exposures[x]]])))<3, "binary", "continuous")


      for (y in 1:length(exposure_time_pts)){

        #gathering all imputed data set lists for a given exposure, outcome, and time point
        exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposures[x], "-", outcomes[z], "-", exposure_time_pts[y]),
                                                                                 names(all_post_balance_corrs))]]
        exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[sapply(strsplit(names(all_post_balance_corrs), "-"),"[",3)==exposure_time_pts[y] &
                                                                             grepl(exposures[x], sapply(strsplit(names(all_post_balance_corrs), "-"),"[",1)) &
                                                                             sapply(strsplit(names(all_post_balance_corrs), "-"),"[",2)==outcomes[z]]]
        # exposure_list=all_post_balance_corrs[names(all_post_balance_corrs)[grepl(paste0(exposures[x], "_", time_pts[time_pts!=time_pts[y]]), names(all_post_balance_corrs))==F]] #making sure other time pts not included for nested cases (e.g., 15 and 154)

        #for binary exposures, evaluate standardized mean differnce between groups (0.1 or user-determined)
        if (exposure_type=="binary"){
          exposureA=dplyr::bind_rows(exposure_list, .id = "imp")
          exposureA=exposureA%>%dplyr::select(c(contains("0"), rn, imp))%>%
            dplyr::select(c(rn, contains("std")))%>%
            dplyr::group_by(rn)%>%
            dplyr::summarise_all("mean")

          exposureB=dplyr::bind_rows(exposure_list, .id = "imp")
          exposureB=exposureB%>%dplyr::select(c(contains("1"), rn, imp))%>%
            dplyr::select(c(rn, contains("std")))%>%
            dplyr::group_by(rn)%>%
            dplyr::summarise_all("mean")

          std_mean_diff=cbind(exposureA, exposureB)
          std_mean_diff$std_mean_diff=std_mean_diff[,2]-std_mean_diff[,4]
          std_mean_diff=std_mean_diff[,colnames(std_mean_diff)[!duplicated(colnames(std_mean_diff))]]

          all_post_balance_correlations[[paste0(exposures[x],"_", exposure_time_pts[y], "-", outcomes[z])]] <- std_mean_diff

          sig=std_mean_diff[abs(std_mean_diff$std_mean_diff)>balance_thresh,]
          colnames(sig)= c("covariate", apply(expand.grid("corr_imp", 1:m), 1, paste, sep="", collapse="_"), "mean_corr")
          sig_unbalanced_averaged_k[[paste("sig_unbalanced_vars_",  exposures[x], exposure_time_pts[y], "-", outcomes[z],"_", sep="")]] <- sig

          significant_corrs_remaining=rbind(significant_corrs_remaining, sig)


        } else if (exposure_type=="continuous"){
          #for continuous exposures, evaluate correlations

          exposure_time_corrs=do.call(cbind, exposure_list) #cbinds them
          names=exposure_time_corrs[,1] #extracts variable names
          exposure_time_corrs=as.data.frame((exposure_time_corrs))
          exposure_time_corrs=exposure_time_corrs[,unlist(lapply(exposure_time_corrs, is.numeric), use.names = FALSE)] #extracts numeric columns (i.e., correlations)
          exposure_time_corrs$mean_cor=rowMeans(exposure_time_corrs) #takes row means to find average correlation
          exposure_time_corrs=cbind(names, exposure_time_corrs) #adds names back


          #compiles all post-balance correlations
          all_post_balance_correlations[[paste0(exposures[x],"_", exposure_time_pts[y], "-", outcomes[z])]] <- exposure_time_corrs

          # browser()

          #gathers only remaining mean correlations over 0.12--the ones left unbalanced
          sig=exposure_time_corrs[abs(exposure_time_corrs$mean_cor)>balance_thresh,]
          colnames(sig)= c("covariate", apply(expand.grid("corr_imp", 1:m), 1, paste, sep="", collapse="_"), "mean_corr")
          sig_unbalanced_averaged_k[[paste("sig_unbalanced_vars_",  exposures[x], exposure_time_pts[y], "-", outcomes[z],"_", sep="")]] <- sig

          significant_corrs_remaining=rbind(significant_corrs_remaining, sig)

        }
      }


      #compiles all remaining correlations over user-specified threshold for each treatment (including all time points) to print and save
      # if there are signfiicant unbalanced covariates
      significant_corrs_remaining=significant_corrs_remaining[order(significant_corrs_remaining$covariate),]
      significant_corrs_remaining=significant_corrs_remaining[!duplicated(significant_corrs_remaining$covariate),] #pulls unique covariates for inspection

      if (nrow(significant_corrs_remaining)>0){
        cat(paste0("USER ALERT: Inspect the following list of unbalanced covariates for exposure ", exposures[x], "-", outcomes[z], " across all exposure time points:"),"\n")
        cat("\n")
        cat(knitr::kable(significant_corrs_remaining), sep="\n")
        cat("\n")


      }else { #if there are no significant unbalanced covariates
        cat(paste0("USER ALERT: There are no unbalanced covariates for exposure ", exposures[x], "-", outcomes[z]),"\n")
        # write.csv(significant_corrs_remaining, file=paste(home_dir, "balance/post-balance correlation values/all_sig_post_balance_cors_", exposures[x], "-", outcomes[z],".csv", sep=""))
        cat("\n")

      }



      covariates_for_model=unique(c(as.character(unlist(significant_corrs_remaining[,1]))))
      #re-naming any factors (level appended on by CBPS)
      covariates_for_model[grepl(paste(factor_covariates, collapse="|"), paste0(covariates_for_model))]=
        substring(covariates_for_model[grepl(paste(factor_covariates, collapse="|"), paste0(covariates_for_model))], 1, nchar(covariates_for_model[grepl(paste(factor_covariates, collapse="|"), paste0(covariates_for_model))])-1)
      covariates_for_model=covariates_for_model[!grepl(exposures[x], covariates_for_model)] #removing exposure (as this will be modeled explicitly)
      covariates_for_model=covariates_for_model[! covariates_for_model %in% exposures[x]]
      covariates_for_model=paste0(covariates_for_model, sep="", collapse=" + ")
      covariates_for_model=gsub(".", "_", covariates_for_model, fixed=TRUE)

      unbalanced_covariates_for_models[[paste0(exposures[x], "-", outcomes[z])]] <-covariates_for_model

    }
  }

  # browser()

  #makes all corrs/std mean differences into data frame to sort and save out for easy viewing
  for (b in 1:length(all_post_balance_correlations)){
    f=all_post_balance_correlations[b]
    corrs=f[[1]]

    if(exposure_type=="continuous"){
      colnames(corrs)=c("Covariate", apply(expand.grid("Imp", seq(m)),1, paste0, sep="", collapse="."), "Mean Corr")
    }else{ colnames(corrs)=c("Covariate", "Mean1", "Mean2", "Std Mean Diff")}

    corrs$exposure_time=sapply(strsplit(names(f), "-"), "[", 1)
    corrs$outcome=sapply(strsplit(names(f), "-"), "[", 2)
    all_corrs_df=rbind(all_corrs_df, corrs)
  }

  if(exposure_type=="continuous"){
    all_corrs_df=all_corrs_df[order(abs(all_corrs_df$`Mean Corr`), all_corrs_df$outcome, decreasing = T),]
  }else{
    all_corrs_df=all_corrs_df[order(abs(all_corrs_df$`Std Mean Diff`), all_corrs_df$outcome, decreasing = T),]
  }

  #save out correlations/std mean differences
  sink(paste0(home_dir, "balance/comparison values/all_post-balance_correlations.html"))
  stargazer::stargazer(all_corrs_df,type="html", digits=2, column.labels = colnames(all_corrs_df),summary=FALSE, rownames = FALSE, header=FALSE,
                       out=paste0(home_dir, "balance/comparison values/all_post-balance_correlations.html"))
  sink()

  cat(paste0("USER ALERT: Check 'balance/comparison values/' folder for a table of all post-balance correlations or standardized mean differences and averaged across imputed datasets"),"\n")
  # cat("\n")


  test=as.data.frame(unbalanced_covariates_for_models)
  #save out correlations
  sink(paste0(home_dir, "balance/unbalanced covariates/all_unbalanced_covariates.html"))
  stargazer::stargazer(test,type="html", digits=2, column.labels = colnames(test),summary=FALSE, rownames = FALSE, header=FALSE, no.space=F,column.sep.width = "5pt",

                       out=paste0(home_dir, "balance/unbalanced covariates/all_unbalanced_covariates.html"))
  sink()

  # cat("\n")
  cat(paste0("USER ALERT: Check 'balance/unbalanced covariates/' folder for table of unbalanced covariates for each exposure-outcome pair based on the cutoff correlation or standardized mean difference of ", as.character(balance_thresh), " averaged across imputed datasets"),"\n")
  cat("\n")


  return(unbalanced_covariates_for_models)

}


