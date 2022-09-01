#code to condense weights by multiplying across time point and averaging across imputed dataset

#' condenseWeights
#' @param ID person-level identifier in your dataset
#' @param home_dir path to home directory for the project
#' @param m number of imputed datasets from Amelia
#' @param weights_models output from createWeights
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param time_pts list of time points along your developmental path of interest for which you have at least one measurement
#' @return data_for_model_with_weights
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom matrixStats rowProds

#' @examples
condenseWeights <-function(ID, home_dir, m, weights_models, exposures, time_pts){

  ####Multiplies weights for each tx across all time points so you end up with one weight per person across time for each imputed dataset.
  #Saves histogram of weights for each tx (across all time points), saves new weights, and assigns new weights to treat_k_mean_weights for each imputed dataset. Set PLOT_CO
  weights_exposure_k=list()

  #iterates over imputed datasets
  for (k in 1:m){

    #iterates through exposures
    for (y in 1:length(exposures)){

      exposure=exposures[y]

      treat_temp={}
      #Cycles through all time points
      for (x in 1:length(time_pts)){
        time=time_pts[x]

        #Collects all of the weights for each time pt for a given tx
        # treat_temp=as.data.frame(cbind(treat_temp, get(paste("fit_",k, "_", exposure, "_", time, sep=""))$weight))
        weights=weights_models[[paste("fit_", k, "_", exposure, "_", time, sep="")]]
        weights=weights$weight
        treat_temp=cbind(treat_temp, weights)
      }


      colnames(treat_temp)=time_pts
      #Finds the product of weights across all time pts for a given tx
      treat_temp_prod=as.data.frame(matrixStats::rowProds(as.matrix(treat_temp)))


      #Plots histogram of weight product for each tx and imputed dataset
      ggplot2::ggplot(data=as.data.frame(treat_temp_prod), aes(x = as.numeric(unlist(treat_temp_prod)))) +
        geom_histogram(color = 'black', bins = 15)
      ggplot2::ggsave(paste("Hist_imp_", k, "_", exposure, "_ALL_TIMES", ".png", sep=""), path=paste0(home_dir, "combined weights"), height=8, width=14)
      print(paste0("A histograms of weights for imputation ", k, " and exposure ", exposure, " has been saved to the 'combined weights' folder"))

      #Assigns tx weights per imputed dataset
      weights_exposure_k[[paste(exposure, "_", k, "_mean_weight", sep="")]] <-treat_temp_prod
      # assign(paste(exposures, "_", k, "_mean_weight", sep=""),treat_temp_prod)

      # #Writes out new weights per each tx for each imputed dataset.
      write.csv(x=as.data.frame(treat_temp_prod), file=paste(home_dir, "combined weights/Imp_", k, "_", exposure, "_mean_weight.csv", sep=""))
      print(paste0("Weights  for imputation ", k, " and exposure ", exposure, " have been saved as csv files in 'combined weights'"))

    }

  }


  ####Averages weights across all imputed datasets to generate mean weights across time for each treatment, with ID merged
  data_for_model=read.csv(paste0(home_dir, "data_for_final_model.csv"))
  weights_exposure=list()
  all_weights<- as.data.frame(data_for_model[,colnames(data_for_model)==ID])
  colnames(all_weights)=ID

  for (y in 1:length(exposures)){
    exposure=exposures[y]
    temp=(rep(1, nrow(data_for_model)))

    for (k in 1:m){
      temp=cbind(temp, weights_exposure_k[[paste(exposure, "_", k, "_mean_weight", sep="")]]) #col binds all imputed datasets
    }
    treat_mean=rowMeans(temp) #finds mean across rows

    weights_exposure[[paste(exposure, "_", "mean_weight_all_imp", sep="")]] <- treat_mean
    mean_weight=as.data.frame(cbind(data_for_model[,colnames(data_for_model)==ID], treat_mean)) #adds IDs to mean weights
    colnames(mean_weight)=c(ID, paste(exposure, "_weight", sep=""))

    all_weights=merge(all_weights, mean_weight, by=ID)

    write.csv(x=as.data.frame(n), file=paste(home_dir, "final weights/", exposure, "_mean_weight_all_imp.csv", sep=""))
    print(paste0("Weights for exposure ", exposure, " are now saved as a csv file in the 'final weights' folder"))
  }

  data_for_model_with_weights=merge(data_for_model, all_weights, by=ID, all.x=T)
  write.csv(data_for_model_with_weights, paste0(home_dir, "final weights/data_for_model_with_weights.csv"))
  print("USER ALERT: Final dataset including weights for each treatment is now saved in the 'final weights' folder as a csv file")

  return(data_for_model_with_weights)

}
