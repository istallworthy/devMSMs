#' Condense weights across time points and imputed datasets
#'
#' Condenses weights by multiplying across time point and averaging across imputed dataset to return a single weight per exposure per person in data_for_model_with_weights
#'
#' @param object msm object that contains all relevant user inputs
#' @param weights_models output from createWeights
#' @return data_for_model_with_weights
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom matrixStats rowProds
#' @seealso [createWeights()] for more on weights_models input
#' @examples condenseWeights(object, weights_models)
#'
condenseWeights <-function(object, data, weights_models){

  ID=object$ID
  home_dir=object$home_dir
  m=object$m
  exposures=object$exposures
  outcomes=object$outcomes
  exposure_time_pts=object$exposure_time_pts

  ####Multiplies weights for each tx across all time points so you end up with one weight per person across time for each imputed dataset.
  weights_exposure_k=list()


  #iterates over imputed datasets
  for (k in 1:m){

    for (z in seq(length(outcomes))){
      outcome=outcomes[z]

      #iterates through exposures
      for (y in 1:length(exposures)){

        exposure=exposures[y]
        treat_temp={}

        #determining exposure type
        exposure_type=ifelse(length(unique(data[,colnames(data)[colnames(data)==exposure]]))<3, "binary", "continuous")

        if (exposure_type=="continuous"){

          #Cycles through all time points
          for (x in 1:length(exposure_time_pts)){
            time=exposure_time_pts[x]

            #Collects all of the weights for each time pt for a given tx
            weights=weights_models[[paste("fit_", k, "_", exposure, "-", outcome, "_", time, sep="")]]
            weights=weights$weight
            treat_temp=cbind(treat_temp, weights)
          }

          colnames(treat_temp)=exposure_time_pts
          #Finds the product of weights across all time pts for a given tx
          treat_temp_prod=as.data.frame(matrixStats::rowProds(as.matrix(treat_temp)))

          #Plots histogram of weight product for each tx and imputed dataset
          ggplot2::ggplot(data=as.data.frame(treat_temp_prod), aes(x = as.numeric(unlist(treat_temp_prod)))) +
            ggplot2::geom_histogram(color = 'black', bins = 15)
          ggplot2::ggsave(paste("Hist_imp_", k, "_", exposure, "-", outcome, "_ALL_TIMES", ".png", sep=""), path=paste0(home_dir, "combined weights/histograms/"), height=8, width=14)
          cat(paste0("A histogram of weights for imputation ", k, " and exposure ", exposure, " on ", outcome, " has been saved to the 'combined weights/histograms/' folder"),"\n")

          #Assigns tx weights per imputed dataset
          weights_exposure_k[[paste(exposure, "-", outcome, "_", k, "_mean_weight", sep="")]] <-treat_temp_prod
          # assign(paste(exposures, "_", k, "_mean_weight", sep=""),treat_temp_prod)

          # #Writes out new weights per each tx for each imputed dataset.
          write.csv(x=as.data.frame(treat_temp_prod), file=paste(home_dir, "combined weights/values/Imp_", k, "_", exposure, "-", outcome, "_mean_weight.csv", sep=""))
          cat(paste0("Weights  for imputation ", k, " and exposure ", exposure, " on ", outcome, " have been saved as csv files in 'combined weights/values/'"),"\n")


        } else if (exposure_type=="binary"){

          weights_exposure_k=weights_models #original weights already by imputation

          cat(paste0("No additional weights have been saved for the binary exposure ", exposure, " given that the original weights encompassed all exposure time points"), "/n")

        }
      }
    }
  }


  ####Averages weights across all imputed datasets to generate mean weights across time for each treatment, with ID merged
  data_for_model=read.csv(paste0(home_dir, "data_for_final_model.csv"))
  weights_exposure=list()
  all_weights<- as.data.frame(data_for_model[,colnames(data_for_model)==ID])
  colnames(all_weights)=ID

  for (z in seq(outcomes)){
    outcome=outcomes[z]

    for (y in 1:length(exposures)){
      exposure=exposures[y]
      temp=(rep(NA, nrow(data_for_model)))

      for (k in 1:m){
        temp=cbind(temp, weights_exposure_k[[paste(exposure, "-", outcome, "_", k, "_mean_weight", sep="")]]) #col binds all imputed datasets
      }
      treat_mean=rowMeans(temp[,2:ncol(temp)]) #finds mean across rows

      weights_exposure[[paste(exposure,  "-", outcome,"_", "mean_weight_all_imp", sep="")]] <- treat_mean
      mean_weight=as.data.frame(cbind(data_for_model[,colnames(data_for_model)==ID], treat_mean)) #adds IDs to mean weights
      colnames(mean_weight)=c(ID, paste(exposure,  "-", outcome, "_weight", sep=""))

      all_weights=merge(all_weights, mean_weight, by=ID)

      write.csv(x=as.data.frame(all_weights), file=paste(home_dir, "final weights/values/", exposure, "-", outcome, "_mean_weight_all_imp.csv", sep=""))
      cat(paste0("Weights for exposure ", exposure, " on ", outcome, " are now saved as a csv file in the 'final weights/values/' folder"),"\n")
    }
  }

  data_for_model_with_weights=merge(data_for_model, all_weights, by=ID, all.x=T)
  write.csv(data_for_model_with_weights, paste0(home_dir, "final weights/data_for_model_with_weights.csv"))
  cat("USER ALERT: Final dataset including weights for each treatment is now saved in the 'final weights' folder as a csv file","\n")

  return(data_for_model_with_weights)

}
