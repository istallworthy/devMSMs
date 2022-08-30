

#' Title
#'
#' @return
#' @export
#' @importFrom ggplot2 ggplot

#' @examples
condenseWeights <-function(ID, home_dir, exposures, time_pts){

  ####Multiplies weights for each tx across all time points so you end up with one weight per person across time for each imputed dataset.
  #Saves histogram of weights for each tx (across all time points), saves new weights, and assigns new weights to treat_k_mean_weights for each imputed dataset. Set PLOT_CO
  imp_temp={}

  #iterates over imputed datasets
  for (k in 1:5){

    #iterates through exposures
    for (y in 1:length(exposures)){

      exposure=exposures[y]
      treat_temp={}

      #Cycles through all time points
      for (x in 1:length(time_pts)){
        time=time_pts[x]

        #Collects all of the weights for each time pt for a given tx
        treat_temp=as.data.frame(cbind(treat_temp, get(paste("fit_",k, "_", exposure, "_", time, sep=""))$weight))

      }
      colnames(treat_temp)=time_pts
      #Finds the product of weights across all time pts for a given tx
      treat_temp_prod=as.data.frame(rowProds(as.matrix(treat_temp)))


      #Plots histogram of weight product for each tx
      ggplot(data=as.data.frame(treat_temp_prod), aes(x = as.numeric(unlist(treat_temp_prod)))) +
        geom_histogram(color = 'black', bins = 15)
      ggsave(paste("Hist_imp_", k, "_", exposure, "_ALL_TIMES", ".png", sep=""), path=paste0(home_dir, "combined weights", height=8, width=14))

      #Assigns tx weights per imputed dataset
      assign(paste(exposures, "_", k, "_mean_weight", sep=""),treat_temp_prod)

      # #Writes out new weights per each tx for each imputed dataset.
      write.csv(x=as.data.frame(treat_temp_prod), file=paste(home_dir, "combined weights/Imp_", k, "_", exposure, "_mean_weight.csv", sep=""))

    }

  }

  ####Averages weights across all imputed datasets to generate mean weights across time for each treatment, with ID merged
  data_for_model=read.csv(paste0(home_dir, "data_for_model2.csv"))


  for (y in 1:length(exposures)){
    exposure=exposures[y]
    temp=(rep(1, 1292))
    for (k in 1:5){ #should go thru 5 when we have all 5 imputed datasets
      temp=cbind(temp, get(paste(exposure, "_", k, "_mean_weight", sep=""))) #col binds all imputed datasets
    }
    treat_mean=RowMeans(temp) #finds mean across rows

    assign(paste(exposure, "_", "mean_weight_all_imp", sep=""),treat_mean)
    n=as.data.frame(cbind(data_for_model$s_id, treat_mean)) #adds IDs to mean weights


    write.csv(x=as.data.frame(n), file=paste(home_dir, "final weights/", exposure, "_mean_weight_all_imp.csv", sep=""))
  }

  data_for_model=read.csv(paste0(home_dir, "data_for_model2.csv"))

  #makes dataframe with all ID numbers for adding

  #ADDS WEIGHTS FOR EACH TX
  all_weights<- data_for_model %>%
    dplyr::select(s_id)

  for (k in 1:length(exposures)){
    exposure=exposures[k]

    #Make a dataframe with original treatments with missing using any wave
    mean_weight=read_csv(paste(home_dir, "final weights/", exposure, "_mean_weight_all_imp.csv", sep="")) #read in the mean weights for each treatment
    mean_weight$X1=NULL
    colnames(mean_weight)[colnames(mean_weight)=="V1"] <- ID
    colnames(mean_weight)[colnames(mean_weight)=="treat_mean"] <- paste(exposure, "_weight", sep="")

    all_weights=merge(all_weights, mean_weight, by=ID)

  }


  data_for_model=merge(data_for_model, all_weights, by="s_id", all.x=T)


}
