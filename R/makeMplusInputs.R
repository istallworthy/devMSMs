
#' Create MPlus inputs for running weighted growth curve models
#' Code for assessing effects of exposure histories with repeated outcomes
#' @param ID person-level identifier in your dataset
#' @param home_dir path to home directory for the project
#' @param missing missing data marker in your dataset
#' @param exposures list of variables that represent your exposures/treatments of interest
#' @param exposure_epochs data frame containing epochs and correponding time points for exposure histories
#' @param covariates any covariates you wish to include in the model
#' @param reference reference history to compare comparison histories
#' @param outcomes list of variables that represent your outcomes of interest
#' @param hi_cutoff integer for percentile considered "high" for exposure
#' @param lo_cutoff integer for percentile considered "low" for exposure
#' @importFrom dplyr %>%
#' @return none
#'
#' @examples makeMplusInputs(ID, home_dir, missing, exposures, exposure_epochs, covariates, reference, outcomes, hi_cutoff=.75, lo_cutoff=.25)

makeMplusInputs <- function(ID, home_dir,data_file, missing, time_var, exposures, exposure_epochs, covariates=NULL, reference="LLL", outcomes, hi_cutoff=.75, lo_cutoff=.25) {

  #makes MPlus dir if necessary
  if(dir.exists(paste0(home_dir, "for Mplus/"))==F){dir.create(paste0(home_dir, "for Mplus/"))}

  # #making a long dataset with weights
  # weights=read.csv(paste0(home_dir,"final weights/data_for_model_with_weights_cutoff_0.9.csv"))
  # weights=weights%>%
  #   dplyr::select(ID, contains("weight"))
  # data=read.csv(data_path)
  # colnames(data)[colnames(data)==time_var] <- "WAVE"
  # write.csv(data, paste0(home_dir, "for Mplus/data_long_weights.csv"))
  # data_file=paste0(home_dir, "for Mplus/data_long_weights.csv")
  #
  # data[data == missing] <- NA
  #
  # data=data%>%dplyr::select(ID, "WAVE", contains(c(exposures, outcomes)), unlist(covariates$covariates))
  # data=merge(data, weights, by=ID, all.x=T)
  #
  #
  data=read.csv(data_file)
  data=data%>%dplyr::select(ID, "WAVE", contains("treatmn"), contains( c(unlist(covariates$outcomes), exposures)), unlist(covariates$covariates))
  #


  epochs=exposure_epochs$epochs
  #
  # if (length(covariates>0)){
  #   if ("transform" %in% colnames(covariates)){
  #     for (i in seq(nrow(covariates))){
  #       new_var=paste0(covariates$transform[i], covariates$outcome[i])
  #       if (covariates$transform[i]=="ln")
  #       temp=log(data[,colnames(data)==covariates$outcome[i]], 10)
  #     }
  #   }
  # }


  #cycles through outcomes
  for (y in 1:unlist(covariates$outcomes)){
    outcome=outcomes[y]

    #cycling through exposures
    for (e in 1:length(exposures)){

      exposure=exposures[e]

      #calculates the mean value for each exposure for each exposure epoch
      for (e in 1:nrow(exposure_epochs)){
        epoch=exposure_epochs[e,1]
        temp=data.frame(row.names=1:nrow(data))
        new_var=paste0(exposure, "_", epoch)

        #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
        for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
          level=as.numeric(unlist(exposure_epochs[e,2]))[l]
          z=data[,which(grepl(paste0(exposure, "_", level), names(data)))]
          temp=cbind(temp, z)
        }
        #adds a new variable of the exposure averaged within epoch
        data=data%>%dplyr::mutate(!!new_var :=rowMeans(temp, na.rm=T))
      }

      #from truncateWweights
      weights=paste0(exposure, "_weight_cutoff")

      #all possible exposure histories
      exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")

      #gathering epoch information for each exposure for deriving betas
      epoch_info=as.data.frame(rep(exposure, length(epochs)))
      epoch_info$time=epochs
      epoch_info$low=NA
      epoch_info$high=NA
      #cycling through epochs to find hi and lo values of the exposure for each epoch
      for (t in 1:length(epochs)){
        var_name=paste(exposure, epochs[t], sep="_")
        epoch_info$low[t]=as.numeric(quantile(data[,var_name],probs= lo_cutoff, na.rm=T))
        epoch_info$high[t]=as.numeric(quantile(data[,var_name],probs= hi_cutoff, na.rm=T))
      }


      #df specifying mean values for each epoch (main effects)
      test=data.frame(names=paste(substr(exposure_epochs$epochs,1,1), "es", sep="_"))
      for (v in 1:nrow(exposure_epochs)){
        vals=as.character(apply(expand.grid(exposure, as.character(as.numeric(unlist(exposure_epochs$values[v])))), 1, paste0, sep="", collapse=""))
        vals=noquote(paste0("mean(", paste(vals, sep=" ", collapse=" "), ");"))
        vals=paste(test[v,1], "=", vals)
        test[v,2]=vals
      }

      main_effects_terms=test$names

      #gathering all variables --ID, time-varying exposures, wave, covariates, outcomes --what else goes into this list?
      variables=c(apply(expand.grid(exposures, as.character(as.numeric(unlist(exposure_epochs$values)))), 1, paste0, sep="", collapse=""), outcomes)
      variables=as.data.frame(c(ID, variables[order(variables)], covariates, "WAVE")) #add others here?
      variables[nrow(variables),1]=paste0(variables[nrow(variables),1], ";")
      colnames(variables)= "V1"

      #listing out just mean effects
      epoch_info2=as.data.frame(test$V2)
      colnames(epoch_info2)= "V1"

      #specifying 2-way interactions
      temp=as.data.frame(cbind(apply(substr(combn(test$names, 2),1,1), 2, paste, sep="", collapse="_"),
                               apply(combn(test$names, 2), 2, paste, sep="", collapse="*"),
                               rep(";")))
      two_way_ints=as.data.frame(paste0(temp[,1]," = ", temp[,2], temp[,3]))
      colnames(two_way_ints)= "V1"
      two_way_int_terms=temp$V1

      #specifying 3-way interaction; changed to "e_m_l" for consistency
      three_way_ints=paste0(paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_"), "=", paste(test$names, sep="", collapse="*"), ";") #adds interaction
      three_way_int_terms=paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_")

      #variables to use: main effect terms, 2-way ints, 3-way ints, wave, outcomes, covariates --anything else?
      use_var=as.data.frame(c(test$names, paste0(temp[,1]),paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_"), outcome, covariates, "WAVE")) #not quite sure all that goes in this field
      colnames(use_var)= "V1"

      #between-subject variables
      between=as.data.frame(c(test$names, paste0(temp[,1]), paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_"), covariates)) #what else goes here?
      colnames(between)= "V1"

      #specifying outcome model (between-person)
      outcome_on=as.data.frame(paste(as.character(between[1:(nrow(between)-1),1])))
      outcome_on[,2]=paste0("(", apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""), ")")
      outcome_on[,3]=paste(outcome_on[,1], outcome_on[,2])
      outcome_on_final=as.data.frame(c((outcome_on[,3]), unlist(covariates)))
      colnames(outcome_on_final)= "V1"

      #specifying model for the slopes (between-person)
      s_on=as.data.frame(paste(as.character(between[1:(nrow(between)-1),1])))
      s_on[,2]=paste0("(", apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""),"a",")")
      s_on[,3]=paste(s_on[,1], s_on[,2])
      s_on_final=as.data.frame(c((s_on[,3])))
      colnames(s_on_final)= "V1"

      #specifying model for the quadratic effect (between-person)
      qd_on=as.data.frame(paste(as.character(between[1:(nrow(between)-1),1])))
      qd_on[,2]=paste0("(", apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""),"b",")")
      qd_on[,3]=paste(qd_on[,1], qd_on[,2])
      qd_on_final=as.data.frame(c((qd_on[,3])))
      colnames(qd_on_final)= "V1"

      #listing out constraints for outcome, slope, quadratic, and a slope:quadratic interaction? does the order matter?
      histories=toupper(exposure_levels)
      histories=noquote(gsub("-", "", histories))
      constraints=noquote(paste0(paste("new ("), paste0(histories, sep="", collapse=" "), " ",
                                 paste0(noquote(apply(expand.grid(c("s", "q", "sd", "qd"), histories), 1, paste0, sep="", collapse="_")[order(apply(expand.grid(c("s", "q", "sd", "qd"), histories), 1, paste0, sep="", collapse="_"))]), sep="", collapse=" "),
                                 ");", sep=" ", collapse=" "))

      #making contrasts
      #finds reference beta values and formula based on parameters of best-fitting model
      ref_parameters=setNames(data.frame(matrix(ncol = length(histories)+2, nrow = 7)), c("sequence", "parameter", histories))
      ref_parameters$parameter= as.data.frame(unlist((apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""))))
      ref_parameters$sequence=outcome_on[,1]
      ref_parameters=as.data.frame(ref_parameters)

      for (x in 1:nrow(ref_parameters)){
        for (y in 1:length(histories)){
          seq=colnames(ref_parameters[y+2])
          if (grepl("es", ref_parameters[x,1])){ #finding main effects
            epoch_info_row=which(substr(epoch_info$time,1,1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1))
            epoch_info_col=which(substr(colnames(epoch_info), 1,1)==tolower(substr(colnames(ref_parameters[y+2]),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)))))
            beta=epoch_info[epoch_info_row, epoch_info_col]

            ref_parameters[x,(y+2)]=beta
          }else{ #find interaction
            if(length(strsplit(ref_parameters$sequence[x], "_")[[1]])==3){
              a=epoch_info[which(substr(epoch_info$time,1,1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[", 1)),
                           which(substr(colnames(epoch_info), 1,1)==tolower(substr(colnames(ref_parameters[y+2]),which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)))))]
              b=epoch_info[which(substr(epoch_info$time,1,1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[", 2)),
                           which(substr(colnames(epoch_info), 1,1)==tolower(substr(colnames(ref_parameters[y+2]),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)))))]
              c=epoch_info[which(substr(epoch_info$time,1,1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[", 3)),
                           which(substr(colnames(epoch_info), 1,1)==tolower(substr(colnames(ref_parameters[y+2]),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)))))]
              beta=a*b*c
              ref_parameters[x,(y+2)]=beta

            }else{
              a=epoch_info[which(substr(epoch_info$time,1,1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[", 1)),
                           which(substr(colnames(epoch_info), 1,1)==tolower(substr(colnames(ref_parameters[y+2]),which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)))))]
              b=epoch_info[which(substr(epoch_info$time,1,1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[", 2)),
                           which(substr(colnames(epoch_info), 1,1)==tolower(substr(colnames(ref_parameters[y+2]),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)),
                                                                                   which(sapply(strsplit(main_effects_terms,"_"), "[", 1)==sapply(strsplit(ref_parameters$sequence[x], "_"), "[",1)))))]
              beta=a*b
              ref_parameters[x,(y+2)]=beta
            }
          }
        }
      }

      #creating forms of parameters and beta weights for each history
      forms=setNames(data.frame(matrix(ncol = 6, nrow = length(histories))), c("history", "form"))
      for (x in 1:length(histories)){
        #outcome model
        form=noquote(paste0(histories[x], " = p1+", (paste0("(", as.character(unlist(ref_parameters[,2])), "*",
                                                            as.character(unlist(ref_parameters[,colnames(ref_parameters)[colnames(ref_parameters)==histories[x]]])), ")",
                                                            sep="", collapse=" + ")), ";", sep="", collapse=""))
        forms[x,1]=histories[x]
        forms[x,2]=form

        #slope model --not sure why this only includes p1-p4? will this always be the case?
        forms[x,3]=paste0("s_", histories[x])
        form_s=noquote(paste0(paste0("s_", histories[x]), " = p1a+", (paste0("(", as.character(unlist(ref_parameters[1:3,2])), "a*",
                                                                             as.character(unlist(ref_parameters[,colnames(ref_parameters)[colnames(ref_parameters)==histories[x]]][1:3])), ")",
                                                                             sep="", collapse=" + ")), ";", sep="", collapse=""))
        forms[x,4]=form_s

        #quadratic model--not sure why this only includes p1-p4?
        forms[x,5]=paste0("q_", histories[x])
        form_q=noquote(paste0(paste0("q_", histories[x]), " = p1b+", (paste0("(", as.character(unlist(ref_parameters[1:3,2])), "b*",
                                                                             as.character(unlist(ref_parameters[,colnames(ref_parameters)[colnames(ref_parameters)==histories[x]]][1:3])), ")",
                                                                             sep="", collapse=" + ")), ";", sep="", collapse=""))
        forms[x,6]=form_q
      }
      outcome_forms=as.data.frame(forms[,2])
      colnames(outcome_forms)="V1"
      s_forms=as.data.frame(forms[,4])
      colnames(s_forms)="V1"
      q_forms=as.data.frame(forms[,6])
      colnames(q_forms)="V1"


      #contrasts using ref event
      d_contrasts=as.data.frame(paste0(paste0(as.character(unlist(strsplit(paste0("d", "_", tolower(histories[!histories %in% reference]), sep=" ", collapse=""), " "))), " = ",
                                              apply(expand.grid(tolower(histories[!histories %in% reference]), tolower(reference)),1, paste0, sep="", collapse=" - ")), ";"))
      colnames(d_contrasts)="V1"

      sd_contrasts=as.data.frame(paste0(paste0(as.character(unlist(strsplit(paste0("sd", "_", tolower(histories[!histories %in% reference]), sep=" ", collapse=""), " "))), " = ",
                                               apply(expand.grid(paste0("s_",tolower(histories[!histories %in% reference])), paste0("s_",tolower(reference))),1, paste0, sep="", collapse=" - ")), ";"))
      colnames(sd_contrasts)="V1"

      qd_contrasts=as.data.frame(paste0(paste0(as.character(unlist(strsplit(paste0("qd", "_", tolower(histories[!histories %in% reference]), sep=" ", collapse=""), " "))), " = ",
                                               apply(expand.grid(paste0("q_",tolower(histories[!histories %in% reference])), paste0("q_",tolower(reference))),1, paste0, sep="", collapse=" - ")), ";"))
      colnames(qd_contrasts)="V1"


      empty_line=data.frame(matrix(nrow = 1, ncol = 1))
      colnames(empty_line)="V1"


      #creating the .inp text for the Mplus file
      model_input <- as.data.frame(NULL)
      model_input[1,1] <- paste0("TITLE: ", exposure, " ", outcome)

      model_input[3,1] <- paste(" ")

      model_input[4,1] <- paste0("   DATA: file is ", data_file)
      model_input=rbind(model_input, empty_line)

      model_input[6,1] <- paste("   define:")
      model_input[7,1] <- paste("   Q = TIME*TIME;")
      model_input[8,1] <- paste("   C = TIME*TIME*TIME;")
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, epoch_info2)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, two_way_ints)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, three_way_ints)
      model_input=rbind(model_input, empty_line)

      model_input= rbind(model_input, paste("   VARIABLE:   NAMES ="))
      model_input= rbind(model_input, variables)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste0("   MISSING ARE ALL (", missing, ");"))
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste0("  USEVAR ="))
      model_input=rbind(model_input, use_var)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, ";")
      model_input=rbind(model_input, paste0("   within= TIME Q;")) #will time_var always be the same?
      model_input=rbind(model_input, paste("   between ="))
      model_input=rbind(model_input, c(main_effects, two_way_int_terms, three_way_int_terms, covariates))
      model_input=rbind(model_input, between, ";")
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste0("  cluster = ", ID, ";"))
      model_input=rbind(model_input, paste0("  bweight = ", weights, ";")) #truncated weights
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste("   ANALYSIS:  TYPE= twolevel random;"))
      model_input=rbind(model_input, paste("estimator = mlf;"))
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste("  MODEL: "))
      model_input=rbind(model_input, paste("  %within%"))
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, paste0("S |", outcome, " on TIME;" )) #will time_var always be the same?
      model_input=rbind(model_input, paste0("QD |", outcome, " on q;")) #what is q?
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, "  %between%")
      model_input=rbind(model_input, paste0(outcome, ";"))
      model_input=rbind(model_input, paste0("S;")) #always S?
      model_input=rbind(model_input, paste0("QD;"))
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, paste0(outcome, " with S QD;"))
      model_input=rbind(model_input, paste0("S with QD;")) #always?
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, paste0(outcome, " on"))
      model_input=rbind(model_input, outcome_on_final)
      model_input=rbind(model_input, ";")
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, paste("s on"))
      model_input=rbind(model_input, s_on_final)
      model_input=rbind(model_input, ";")
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, paste("qd on"))
      model_input=rbind(model_input, qd_on_final)
      model_input=rbind(model_input, ";")
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, paste0("[", outcome, "] (p1);"))
      model_input=rbind(model_input, paste0("[s] (p1a);"))
      model_input=rbind(model_input, paste0("[qd] (p1b);"))
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste("model constraint:"))
      model_input=rbind(model_input, constraints)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, outcome_forms)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, s_forms)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, q_forms)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input,d_contrasts)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input,sd_contrasts)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input,qd_contrasts)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste("    OUTPUT: SAMPSTAT residual STDYX modindices (3.94) ;"))
      model_input=rbind(model_input, paste("    tech1 tech8 tech11 tech14;"))
      model_input=rbind(model_input, paste("    plot:"))
      model_input=rbind(model_input, paste("    type is plot1;"))
      model_input=rbind(model_input, empty_line)


      #write out mplus input file
      write.inp.file(model_input,fixPath(file.path(paste0(home_dir, "for Mplus/", exposure, "_", outcome, ".inp",sep=""))))

      # write.inp.file(round2input,fixPath(file.path(dir,"round2calibration.inp",sep="")))

    }
  }



  print("Please carefully inspect all input files in the 'for Mplus' folder prior to running them")

}
