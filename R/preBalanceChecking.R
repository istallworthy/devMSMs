
#examining initial imbalance for each exposure, each time point, each covariate for each history

preBalanceChecking <- function(object, forms, imputed_datasets, data_for_model_with_weights){

  ID=object$ID
  home_dir=object$home_dir
  exposures=object$exposures
  exposure_time_pts=object$exposure_time_pts
  outcomes=object$outcomes
  hi_cutoff=object$hi_cutoff
  lo_cutoff=object$lo_cutoff
  balance_thresh=object$balance_thresh

  #values determining presence and absence of exposure for continuous exposures
  hi_val=quantile(data_long$InRatioCor, 0.60)
  lo_val=quantile(data_long$InRatioCor, 0.40)


  library(cobalt)
  library(tidyr)

  #does this need to be done for all imputed datasets?
  data=wide_long_datasets[[paste("imp", 1, "_widelong", sep="")]]
  data[,ID]=as.numeric(levels(data[,ID]))

  data_long=complete(imputed_datasets,1)

  #data frame with all sampling weights for all exposures at all exposure time points for all histories
  all_prop_weights=data.frame(id=NA,exposure=NA, exp_time=NA, history=NA,prop_weight=NA)
  colnames(all_prop_weights)[colnames(all_prop_weights)=="id"]=ID


  for (y in 1:length(outcomes)){
    outcome=outcomes[y]

    for (x in 1:length(exposures)){
      exposure=exposures[x]
      exposure_type=ifelse(length(unique(data[,colnames(data)[colnames(data)==paste0(exposure,".", exposure_time_pts[1])]]))<3, "binary", "continuous")

      #gathering all balnace stats over all exposure time points, covariates and lagged covariates, and exposure histories
      bal_stats=data.frame(exposure=NA, exp_time=NA,covariate=NA,covar_time=NA,history=NA,balance_stat=NA,balanced=NA)

      for (z in 1:length(exposure_time_pts)){
        exposure_time_pt=exposure_time_pts[z]

        full_form= forms[[paste0("form_", exposure, "-", outcome, "-", exposure_time_pt)]]

        lagged_time_pts=exposure_time_pts[exposure_time_pts<exposure_time_pt]

        #if there are histories
        if(length(lagged_time_pts)>0){
          #creating proportion weights based on proportion of individuals in a given exposure history
          prop_weights=data.frame(id=data[,ID], exposure=exposure,exp_time=exposure_time_pt,history=NA,prop_weight=NA)
          colnames(prop_weights)[colnames(prop_weights)=="id"]=ID

          #identifying historis
          # histories=expand.grid(apply(expand.grid(lagged_time_pts, c(0,1)), 1, paste, sep="", collapse="_"),
          #                       apply(expand.grid(lagged_time_pts, c(0,1)), 1, paste, sep="", collapse="_"))

          histories=apply(gtools::permutations(2, length(lagged_time_pts), c(1,0), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")
          histories=as.data.frame(histories)
          # gtools::combinations(4, nrow(lagged_time_pts), lagged_time_pts, repeats.allowed=TRUE)

          # #removing duplicates
          # if (length(lagged_time_pts)==1){
          #   histories=as.data.frame(histories[!duplicated(histories[,1]),1])
          # }else{
          #   histories=histories[sapply(strsplit(paste0(histories[,1]), "_"), "[", 1)!=sapply(strsplit(paste0(histories[,2]), "_"), "[", 1),] #removes dupes
          #   histories=as.data.frame(histories)
          #   histories=histories[ !duplicated(apply(histories, 1, sort), MARGIN = 2), ]
          # }

          #cycle through histories
          for (h in 1:nrow(histories)){ #cyles through each history

            # times=ncol(histories)
            his=as.data.frame(histories[h,])
            his=as.character(unlist(his[1,]))

            #get wide data variable names for exposures at each time in each history
            # exps_time=apply(expand.grid(exposure, c(sapply(strsplit(his, "-"), "[",1))), 1, paste, collapse=".", sep=".")
            exps_time=apply(expand.grid(exposure, as.character(lagged_time_pts)), 1, paste, collapse=".", sep=".")

            data$flag=0

            for (t in 1:length(lagged_time_pts)){ #cycles through times in each history and flags individuals who have given history
              # time=sapply(strsplit(as.character(his[t]), "_"), "[", 1)
              time=lagged_time_pts[t]
              exp=as.numeric(sapply(strsplit(as.character(his), "-"), "[", t))
              flag=t-1
              #
              if (exp==0){ #low levels
                if (exposure_type=="continuous"){
                  data$flag=ifelse(data[, exps_time[t]]<lo_val & data$flag==flag, t, NA) #findings those flagged at prev t's and meeting new criteria
                }else{
                  data$flag=ifelse(data[, exps_time[t]]==0 & data$flag==flag, t, NA)
                }
              }
              if (exp==1){ #hi levels
                if (exposure_type=="continuous"){
                  data$flag=ifelse(data[, exps_time[t]]>hi_val & data$flag==flag, t, NA)#findings those flagged at prev t's and meeting new criteria
                }else{
                  data$flag=ifelse(data[, exps_time[t]]==1 & data$flag==flag, t, NA)
                }
              }
            } #ends history time points

            prop=sum(data$flag==t, na.rm=T)/nrow(data)
            ids=data[data$flag==t, ID]#finding ids with that history
            data$flag=NULL

            #gather proportions for each history
            prop_weights[prop_weights[,ID] %in% ids, "prop_weight"] <- prop
            prop_weights[prop_weights[,ID] %in% ids, "history"] <- paste(his, collapse=",")

          } #ends history loopp

          #summarize contributors to each history
          sum=prop_weights%>%dplyr::group_by(as.factor(history))%>%dplyr::summarize(prop=mean(prop_weight, na.rm=T),
                                                                                    n=dplyr::n())

          cat(paste0("For exposure ", exposure, " at time poit ", exposure_time_pt, ", ", round(sum(sum$prop, na.rm=T)*100, 1), " % of individuals contributed to pre-balance checking in the following exposure histories across times: ", paste(lagged_time_pts, collapse=" ")),
              "\n")
          print(na.omit(sum))

          if(sum(sum$n==1)>0){
            cat(paste0("USER ALERRT: There are is only one individual with one or more exposure histories, based on the high and low cutoffs provided for the continuous exposure, ", exposure, ". Please provide different high and low cutoff values in the msmObject and re-run."), "\n")
          }


          #create dataset w/ sampling weights for assessing balance at this exposure time point; will only have ppl who contributed to histories
          temp=merge(data, na.omit(prop_weights), by=ID, all.x=F)
          temp=na.omit(temp)

          cat(paste0(nrow(temp), " individuals contribute to pre-balance checking for exposure ", exposure, " at time point ", exposure_time_pt), "\n")

          balance=bal.tab(full_form, data=temp, stats = c("c"), thresholds = c(cor=0.1), un=T,continuous="std", binary="std",
                          which.time = .none, cluster = "history", s.weights="prop_weight")


          for (m in 1:nrow(histories)){ #gets balance stats for all histories
            results=balance$Cluster.Balance[[m]]
            bal_stats=rbind(bal_stats,
                            data.frame(exposure=exposure, exp_time=exposure_time_pt,covariate=rownames(results$Balance[1]),covar_time=sapply(strsplit(rownames(results$Balance[1]), "\\."), "[", 2),
                                       history=paste(as.character(unlist(histories[m,])), collapse=" "), balance_stat=results$Balance$Corr.Un, balanced=results$Balance$R.Threshold.Un))
          }



        }else{
          histories=NA #if it's the first time point and there are no histories

          #assess pre-balance
          balance=bal.tab(full_form, data=data, stats = c("c"), thresholds = c(cor=0.1), un=T,continuous="std", binary="std",
                          which.time = .all)

          bal_stats=rbind(bal_stats, data.frame(exposure=exposure,exp_time=exposure_time_pt,covariate=rownames(balance$Balance[1]),covar_time=NA,
                                                history=NA, balance_stat=balance$Balance$Corr.Un,balanced=balance$Balance$R.Threshold.Un))
        }
        #collect proportions for histories at this time point
        all_prop_weights=rbind(all_prop_weights, prop_weights)
      } #ends exp_time_pt



      #summarize balancing stats
      bal_summ=bal_stats%>%dplyr::group_by(exposure, exp_time, covariate, covar_time)%>%
        dplyr::summarize(mean_bal_stat=mean(balance_stat), n=dplyr::n())%>%
        dplyr::mutate(balanced=ifelse(abs(mean_bal_stat)<balance_thresh, 1, 0))

      write.csv(bal_stats, paste0(home_dir, "pre-balance/", exposure, "_balance_stats.csv"))
      cat(paste0("Pre-balance statistics for ", exposure, " have been saved in the pre-balance/ folder"), "\n")

    } #ends exposure loop
  }#ends outcome loop

  write.csv(all_prop_weights, paste0(home_dir, "pre-balance/history_sample_weight.csv"))
  cat(paste0("Pre-balance sampling weights have been saved in the pre-balance/ folder"), "\n")


}


