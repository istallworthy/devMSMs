
#code to calculate balance stats based on Jackson paper (either weighted or unweighted)
calcBalStats <-function(object, imp_wide_data, forms, f_type, exposure, outcome, k=1, weighted=0, histories=1){

  # #for testing
  # outcome=object$outcomes
  # exposure="InRatioCor"
  # imp_wide_data=wide_long_datasets
  # imp_wide_data=data_for_model_with_weights
  # forms=full_forms

  # k=1
  # weighted=1
  # long_data=imputed_datasets
  #

  # browser()
  ID=object$ID
  home_dir=object$home_dir
  exposure_time_pts=object$exposure_time_pts
  balance_thresh=object$balance_thresh
  factor_covariates=object$factor_covariates

  weights_method=ifelse(weighted==1, object$weights_method, "no weights")
  h=histories
  histories=NULL

  library(cobalt)
  # library(tidyr)

  set.seed(1234)

  form_name=f_type

  #creating initial data frames
  #data frame with all sampling weights for all exposures at all exposure time points for all histories
  all_prop_weights=data.frame(id=NA,exposure=NA, exp_time=NA, history=NA)
  colnames(all_prop_weights)[colnames(all_prop_weights)=="id"]=ID

  #gathering all balance stats over all exposure time points, covariates and lagged covariates, and exposure histories
  # bal_stats=data.frame(exposure=NA, exp_time=NA,covariate=NA,covar_time=NA,history=NA,balance_stat=NA,balanced=NA, prop_weight=NA, n=NA)
  # bal_stats=data.frame(exposure=NA, exp_time=NA,covariate=NA,covar_time=NA,balance_stat=NA,balanced=NA, prop_weight=NA, n=NA)
  all_bal_stats=data.frame()

  # browser()

  #wide data
  data= tryCatch({
    imp_wide_data[[names(imp_wide_data)[grepl(paste0("fit_", k, "_", exposure, "-", outcome), names(imp_wide_data))]]]
  }, error=function(e){
    # imp_wide_data[[names(imp_wide_data)[grepl(paste0("imp", k), names(imp_wide_data))]]]
    # imp_wide_data[[names(imp_wide_data)[grepl(k, names(imp_wide_data))]]]
    imp_wide_data[[k]]
    # mice::complete(imp_wide_data,k)
  })

  # data=wide_long_datasets[[paste("imp", 1, "_widelong", sep="")]]
  data[,ID]=as.numeric(levels(data[,ID]))
  # data_long=tidyr::complete(long_data,k)
  # data$history=NA
  # numeric_vars=numeric_vars[numeric_vars %in% colnames(data)]
  # data[,numeric_vars] <- lapply(data[,numeric_vars] , as.numeric)


  folder=ifelse(weighted==0, "pre-balance/", "balance/")


  if (weighted==1 & sum(grepl("weight", colnames(data)))==0){
    stop('You have elected the weighted option but have not provided any weights in the wide dataset')}

  exposure_type=ifelse(length(unique(data[,paste0(exposure,".", exposure_time_pts[1])]))<3, "binary", "continuous")


  for (z in 1:length(exposure_time_pts)){ #cycles through exposure time points
    exposure_time_pt=exposure_time_pts[z]

    lagged_time_pts=exposure_time_pts[exposure_time_pts<exposure_time_pt]

    #GETS COVARIATES FROM FORM FOR ASSESSING BALANCE
    full_form= forms[[names(forms)[grepl(paste0("form_", exposure, "-", outcome, "-", exposure_time_pt), names(forms))]]]
    covars=paste(deparse(full_form[[3]], width.cutoff = 500), collapse="")
    covar_time=sapply(strsplit(unlist(strsplit(as.character(covars), "\\+")), "\\."), "[", 2)
    covars=as.character(unlist(strsplit(covars, "\\+")))
    covars=gsub(" ", "", covars)

    #GETTING BALANCE STATS FOR T=1 W/ NO EXPOSURE HISTORY (ok to standardize)
    if(length(lagged_time_pts)==0){
      temp=data

      # unweighted pre-balance checking
      if(weighted==0){
        if (exposure_type=="continuous"){
          bal_stats= cobalt::col_w_cov(temp[, c(covars)], temp[, paste0(exposure, ".", exposure_time_pt)], std=T) #finding correlation
        }
        if (exposure_type=="binary"){
          bal_stats=cobalt::col_w_smd(temp[, c(covars)], temp[, paste0(exposure, ".", exposure_time_pt)], std=T) #finding smd
        }
      }

      #weighted balance checking
      if(weighted==1){
        if(exposure_type=="continuous"){
          bal_stats= cobalt::col_w_cov(temp[, c(covars)], temp[, paste0(exposure, ".", exposure_time_pt)], std=T, #finding cor
                                       weights=temp[,"weights"])
        }
        if (exposure_type=="binary")
          bal_stats=cobalt::col_w_smd(temp[, c(covars)], temp[, paste0(exposure, ".", exposure_time_pt)], std=T, #finding smd
                                      weights=temp[,"weights"])
      }

      bal_stats=as.data.frame(bal_stats)
      colnames(bal_stats)="std_bal_stats"
      bal_stats=bal_stats%>%dplyr::mutate(balanced=ifelse(abs(std_bal_stats)<balance_thresh,1,0)) #compare to balance threshold

    }



    #ASSIGNING HISTORIES FOR EXP TIME POINTS T>1
    #if there are histories, identify sampling weights based on history up until t
    if(length(lagged_time_pts)>0){
      #creating proportion weights based on proportion of individuals in a given exposure history
      prop_weights=data.frame(id=data[,ID], exposure=exposure,exp_time=exposure_time_pt,history=NA)
      colnames(prop_weights)[colnames(prop_weights)=="id"]=ID

      #finding histories up until exp time point T
      histories=apply(gtools::permutations(2, length(lagged_time_pts), c(1,0), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")
      histories=as.data.frame(histories)

      #cycle through histories to identify individuals with each history
      for (h in 1:nrow(histories)){
        his=as.data.frame(histories[h,])
        his=as.character(unlist(his[1,]))

        #get wide data variable names for exposures at each time in each history
        exps_time=apply(expand.grid(exposure, as.character(lagged_time_pts)), 1, paste, collapse=".", sep=".")
        data$flag=0
        for (t in 1:length(lagged_time_pts)){ #cycles through times in each history and flags individuals who have given history
          time=lagged_time_pts[t]
          exp=as.numeric(sapply(strsplit(as.character(his), "-"), "[", t))
          flag=t-1
          if (exp==0){ #low levels
            if (exposure_type=="continuous"){
              data$flag=ifelse(data[, exps_time[t]]<=median(data[,paste0(exposure, ".", exposure_time_pt)]) & data$flag==flag, t, NA) #findings those flagged at prev t's and meeting new criteria
            }else{
              data$flag=ifelse(data[, exps_time[t]]==0 & data$flag==flag, t, NA)
            }
          }
          if (exp==1){ #hi levels
            if (exposure_type=="continuous"){
              data$flag=ifelse(data[, exps_time[t]]>median(data[,paste0(exposure, ".", exposure_time_pt)]) & data$flag==flag, t, NA)#findings those flagged at prev t's and meeting new criteria
            }else{
              data$flag=ifelse(data[, exps_time[t]]==1 & data$flag==flag, t, NA)
            }
          }
        } #ends history time points
        ids=data[data$flag==t, ID]#finding ids with that history
        prop_weights[prop_weights[,ID] %in% ids, "history"] <- paste(his, collapse=",")
        data$flag=NULL
      } #ends history loop

      #summarize contributors to each history
      prop_sum=prop_weights%>%dplyr::group_by(as.factor(history))%>%dplyr::summarize(n=dplyr::n())

      if(h==1){
        cat(paste0("For exposure ", exposure, ", imputation ", k, " at time point ", exposure_time_pt, ", individuals were distributed within histories across lagged time point(s) ", paste(lagged_time_pts, collapse=" "), " as follows:"),
            "\n")
        prop_sum$`as.factor(history)`=as.character(prop_sum$`as.factor(history)`)
        colnames(prop_sum)=c("history","n")
        # print(prop_sum)
        cat(knitr::kable(prop_sum, caption="History Distribution for Weighting Balance Statistics", format='pipe'),  sep="\n")
        cat("\n")
      }




      #GET WEIGHTED BALANCE STATISTICS FOR T>1
      temp=merge(data, prop_weights, by=ID, all.x=T)

      #removing any histories with n=1 (cannot calc bal stats)
      if (sum(prop_sum$n==1)>0){
        cat(paste0("USER ALERT: the following history/histories, ", as.character(as.data.frame(prop_sum)[prop_sum$n==1, 1]), ", has/have been ommitted from balance checking for exposure ", exposure,
                   ",imputation ", k, ", at time point ", exposure_time_pt), "\n")
        temp=temp[temp$history!=as.character(as.data.frame(prop_sum)[prop_sum$n==1, 1]),]
      }

      #unweighted pre-balance checking
      if(weighted==0){ #no weighting
        if (exposure_type=="continuous"){
          bal_stats= sapply(sort(unique(temp$history)) , function(i){ #finding balance by history
            temp2=temp%>%dplyr::filter(history==i)
            cobalt::col_w_cov(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std=F, #finding covariance
                              subset=temp2$history[temp2$history==i]==i)
          })
          #getting weighted mean across histories (cols), weighting by proportion of those w/ that same history
          weighted_bal_stats=sapply(seq(nrow(bal_stats)), function(i){
            weighted.mean(t(bal_stats[i,])[1,],tabulate(as.factor(temp$history))/nrow(temp))})

          bal_stats=as.data.frame(cbind(bal_stats, weighted_bal_stats))

          #standardizing balance statistics after weighting by history
          bal_stats=bal_stats%>%dplyr::mutate(std_bal_stats=weighted_bal_stats/
                                                (sapply(seq(ncol(data[,covars])), function(x){
                                                  sd(as.numeric(data[,covars][,x]), na.rm=T)#unweighted covar sd
                                                })*sd(data[,paste0(exposure, ".", exposure_time_pt)], na.rm=T)))%>% #exposure SD
            dplyr::mutate(balanced=ifelse(abs(std_bal_stats)<balance_thresh,1,0)) #compare to balance threshold
          bal_stats=bal_stats%>%dplyr::select(contains(c("std", "balanced")))

        }

        if (exposure_type=="binary"){
          bal_stats=
            sapply(sort(unique(temp$history)), function(i){
              temp2=temp%>%dplyr::filter(history==i)
              cobalt::col_w_smd(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std=F, #finding covariance
                                subset=temp2$history[temp2$history==i]==i)
            })

          #getting weighted mean across histories, weighting by proportion of those w/ that same history
          weighted_bal_stats=sapply(seq(nrow(bal_stats)), function(i){
            weighted.mean(t(bal_stats[i,])[1,],tabulate(as.factor(temp$history))/nrow(temp))})

          bal_stats=as.data.frame(cbind(bal_stats, weighted_bal_stats))

          #standardizing balance statistics after finding weighted balance stats
          bal_stats=bal_stats%>%dplyr::mutate(std_bal_stats=weighted_bal_stats/
                                                sapply(seq(ncol(data[,covars])), function(x){
                                                  sqrt(mean( #dividing by pool SD estimate (unajusted)
                                                    var(as.numeric(data[paste0(exposure, ".", exposure_time_pts[1])==1,covars[x]])), #treated var
                                                    var(as.numeric(data[paste0(exposure, ".", exposure_time_pts[1])==0,covars[x]]))))}))%>% #untreated var
            dplyr::mutate(balanced=ifelse(abs(std_bal_stats)<balance_thresh,1,0)) #compare to balance threshold
          bal_stats=bal_stats%>%dplyr::select(contains(c("std", "balanced")))
        }
      }


      # browser()
      #weighted balance checking
      if (weighted==1) #if weighted, use IPTW weights from weightitmsm
        if (exposure_type=="continuous"){
          #finds balance for each covariate clustered/subset by history
          bal_stats=
            sapply(sort(unique(temp$history)), function(i){
              temp2=temp%>%dplyr::filter(history==i)
              cobalt::col_w_cov(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std=F, #finding covariance
                                subset=temp2$history[temp2$history==i]==i,
                                weights=temp2[, "weights"])
            })

          #getting weighted mean across histories, weighting by proportion of those w/ that same history
          weighted_bal_stats=sapply(seq(nrow(bal_stats)), function(i){
            weighted.mean(t(bal_stats[i,])[1,],tabulate(as.factor(temp$history))/nrow(temp))})

          bal_stats=as.data.frame(cbind(bal_stats, weighted_bal_stats))

          #standardizing balance statistics after weighting by history
          bal_stats=bal_stats%>%dplyr::mutate(std_bal_stats=weighted_bal_stats/
                                                (sapply(seq(ncol(data[,covars])), function(x){
                                                  sd(as.numeric(data[,covars][,x]), na.rm=T)#unweighted covar sd
                                                })*sd(data[,paste0(exposure, ".", exposure_time_pt)], na.rm=T)))%>% #exposure SD
            dplyr::mutate(balanced=ifelse(abs(std_bal_stats)<balance_thresh,1,0)) #compare to balance threshold

          #for a weighted_bal_stat of 0, make std stat also 0
          bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] =0
          bal_stats$balanced=ifelse(abs(bal_stats$std_bal_stats)<balance_thresh,1,0)

          bal_stats=bal_stats%>%dplyr::select(contains(c("std", "balanced")))
        }

      if (exposure_type=="binary"){
        #finds balance for each covariate clustered/subset by history
        bal_stats=
          sapply(sort(unique(temp$history)), function(i){
            temp2=temp%>%dplyr::filter(history==i)
            cobalt::col_w_smd(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std=F, #finding mean difference
                              subset=temp2$history[temp2$history==i]==i,
                              weights=temp2[,"weights"])
          })

        #getting weighted mean across histories, weighting by proportion of those w/ that same history
        weighted_bal_stats=sapply(seq(nrow(bal_stats)), function(i){
          weighted.mean(t(bal_stats[i,])[1,],tabulate(as.factor(temp$history))/nrow(temp))})

        bal_stats=as.data.frame(cbind(bal_stats, weighted_bal_stats))

        #standardizing balance statistics after finding weighted balance stats
        bal_stats=bal_stats%>%dplyr::mutate(std_bal_stats=weighted_bal_stats/
                                              sapply(seq(ncol(data[,covars])), function(x){
                                                sd(as.numeric(data[,covars][,x]), na.rm=T)*#unweighted covar sd
                                                  sd(data[,paste0(exposure, ".", exposure_time_pts[1])], na.rm=T) #getting sd of first time pt exp only
                                              }))%>% #exposure
          dplyr::mutate(balanced=ifelse(abs(std_bal_stats)<balance_thresh,1,0)) #compare to balance threshold

        #for a weighted_bal_stat of 0, make std stat also 0
        bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] =0
        bal_stats$balanced=ifelse(abs(bal_stats$std_bal_stats)<balance_thresh,1,0)

        bal_stats=bal_stats%>%dplyr::select(contains(c("std", "balanced")))

      }
      #collect proportions for histories at this time point
      all_prop_weights=rbind(all_prop_weights, prop_weights)
    } #ends lag>0 loops



    #ADDS INFO TO BAL STATS
    bal_stats=as.data.frame(bal_stats)
    bal_stats$covariate=rownames(bal_stats)
    #renames factor covariates
    bal_stats$covariate[sapply(strsplit(sapply(strsplit(bal_stats$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates] <-sapply(strsplit(bal_stats$covariate, "_"), "[", 1)[sapply(strsplit(sapply(strsplit(bal_stats$covariate, "_"), "[", 1), "\\."), "[",1) %in% factor_covariates]

    bal_stats=bal_stats%>%dplyr::mutate(
      exposure=exposure, exp_time=exposure_time_pt)%>%
      dplyr::mutate(covar_time=sapply(strsplit(covariate, "\\."), "[", 2))

    all_bal_stats=rbind(all_bal_stats, bal_stats)
    all_bal_stats$covar_time[is.na(all_bal_stats$covar_time)]=0



    x_lab=ifelse(exposure_type=="continuous", "Exposure-Covariate Correlation", "Standardized Mean Difference")
    labels=ifelse(bal_stats$balanced==0, bal_stats$covariate, "")
    min_val=ifelse(min(bal_stats$std_bal_stats)<0, min(bal_stats$std_bal_stats)-0.1, balance_thresh-0.05)
    max_val=ifelse(max(bal_stats$std_bal_stats)>0, max(bal_stats$std_bal_stats)+0.1, balance_thresh+0.05)

    # browser()

    #make love plot per exposure time point
    lp <- ggplot2::ggplot(aes(x =  bal_stats$std_bal_stats, y = bal_stats$covariate), data = bal_stats) +
      ggplot2::geom_point(aes(y = as.factor(bal_stats$covariate),
                              x = bal_stats$std_bal_stats,
                              fill = "white",
                              alpha = 1))+
      ggplot2::geom_text(aes(label=labels, hjust=-0.2, vjust=0.2), size=1.5, color="red")+
      ggplot2::geom_vline(xintercept = balance_thresh,linetype = "dashed", color = "red")+
      ggplot2::geom_vline(xintercept = -balance_thresh,linetype = "dashed", color = "red")+
      ggplot2::xlab(x_lab)+
      ggplot2::ylab("Covariate")+
      ggplot2::xlim(min_val, max_val)+
      ggplot2::ggtitle(paste0(exposure, "(t=", exposure_time_pt, ") Balance"))+
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                     axis.text.x = ggplot2::element_text(color = "black"),
                     axis.text.y = ggplot2::element_text(color = "black"),
                     axis.text = ggplot2::element_text(size=8),
                     # axis.ticks.length.y = unit(0.85, "cm"),
                     panel.border = ggplot2::element_rect(fill = NA, color = "black"),
                     plot.background = ggplot2::element_blank(),
                     plot.title=ggplot2::element_text(size=10),
                     legend.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank(),
                     legend.position="none")+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if(nrow(bal_stats)>40){ #stagger covariate labels if there are many
      lp <- lp+ ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge=2))
    }
    suppressMessages(ggplot2::ggsave(paste0(home_dir, folder, "/plots/", form_name, "_imp_",k, "_", exposure,"_", exposure_time_pt,"_",weights_method, "_summary_balance_plot.jpeg"),
                                     width=6, height=8))
    cat(paste0("A ", gsub("/", "", folder), "  summary plot for ", form_name, " ", exposure,  " imputation ", k," at time ", exposure_time_pt," for weighting method ",weights_method, " has now been saved in the '", folder, "plots/' folder."), "\n")
    cat("\n")
  } #ends exp_time_pt



  bal_summary_exp=all_bal_stats%>%dplyr::group_by(exp_time)%>%
    dplyr::summarize(balanced_n=sum(balanced==1),
                     imbalanced_n=sum(balanced==0),
                     n=dplyr::n())

  # browser()
  # write.csv(bal_summary_exp, paste0(home_dir, folder, exposure, "_",k, "_imbalanced_covars.csv"))


  write.csv(bal_summary_exp, paste0(home_dir, folder, exposure, "_",k,"_",weights_method, "_balance_stat_summary.csv"))
  cat(paste0("Balance statistics for ", form_name, " ",exposure, ", imputation ", k, ", weighting method ",weights_method," have been saved in the '", folder, "' folder"), "\n")


  write.csv(all_prop_weights, paste0(home_dir, folder, exposure, "_", k, "_",weights_method, "_history_sample_weight.csv"))
  cat(paste0("Sampling weights ", "for ", form_name, " ", exposure, ", imputation ", k,", weighting method ",weights_method,  ", weighted=", weighted, " have been saved in the '", folder, "' folder"), "\n")

  cat("\n")

  #GETS total possible COVARIATES FROM FORM FOR ASSESSING BALANCE
  all_form=as.data.frame(do.call(rbind, forms))
  tot_covars=deparse(all_form[,3], width.cutoff = 500)
  tot_covars=as.character(unlist(strsplit(tot_covars, "\\+")))[!grepl("form", as.character(unlist(strsplit(tot_covars, "\\+"))))]
  tot_covars=gsub(" ", "", tot_covars)
  tot_covars=na.omit(sapply(strsplit(tot_covars, "\\."), "[",1)[!duplicated(sapply(strsplit(tot_covars, "\\."), "[",1))])

  cat(paste0("USER ALERT: For exposure ", exposure, " imputation ", k, " using ", weights_method, " , ",
             sum(bal_summary_exp$imbalanced_n, na.rm=T), " out of ", sum(bal_summary_exp$n, na.rm=T),
             " (", round((sum(bal_summary_exp$imbalanced_n, na.rm=T)/sum(bal_summary_exp$n))*100,0), "%) covariates across time points corresponding to ",
             length(sapply(strsplit(all_bal_stats[all_bal_stats$balanced==0, "covariate"], "\\."), "[", 1)[!duplicated( sapply(strsplit(all_bal_stats[all_bal_stats$balanced==0, "covariate"], "\\."), "[", 1))]) ,
             " out of ", length(tot_covars), " unique constructs remain imbalanced with an average absolute value correlation/std mean difference of ",
             round(mean(abs(all_bal_stats[all_bal_stats$balanced==0, "std_bal_stats"]), na.rm = T),2), " (range= ",
             round(min(abs(all_bal_stats[all_bal_stats$balanced==0, "std_bal_stats"]), na.rm = T),2), "-",
             round(max(abs(all_bal_stats[all_bal_stats$balanced==0, "std_bal_stats"]), na.rm = T),2), "), based on the ", form_name, " as shown below:"), "\n")
  cat("\n")
  # print(bal_summary_exp)
  cat(knitr::kable(bal_summary_exp, caption=paste0("Imbalanced Covariates using ", weights_method, " and ", form_name), format='pipe'),  sep="\n")


  rownames(all_bal_stats)<-NULL

  return(all_bal_stats)

}
