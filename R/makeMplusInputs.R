
#' Create MPlus inputs for running weighted growth curve models
#' Code for assessing effects of exposure histories with repeated outcomes
#' @param object msm object that contains all relevant user inputs
#' @param covariates any covariates you wish to include in the model
#' @param reference reference history to compare comparison histories
#' @param outcomes list of variables that represent your outcomes of interest
#' @param hi_cutoff integer for percentile considered "high" for exposure
#' @param lo_cutoff integer for percentile considered "low" for exposure
#' @importFrom dplyr %>%
#' @return none
#'
#' @examples makeMplusInputs(object, mplusObject)

makeMplusInputs <- function(object, mplusObject) {


  ID=object$ID
  home_dir=object$home_dir
  missing=object$missing
  exposures=object$exposures
  exposure_epochs=object$exposure_epochs
  outcomes=covariates$outcomes #gets outcomes from covariate data frame
  covars=covariates$covariates

  data_file=mplusObject$data_file
  reference=mplusObject$reference
  covariates=mplusObject$covariates
  hi_cutoff=mplusObject$hi_cutoff
  lo_cutoff=mplusObject$lo_cutoff


  #makes MPlus dir if necessary
  if(dir.exists(paste0(home_dir, "for Mplus/"))==F){dir.create(paste0(home_dir, "for Mplus/"))}

  data=read.csv(data_file)
  data=as.data.frame(data)
  data[data==missing] <- NA
  #


  epochs=exposure_epochs$epochs


  #cycles through outcomes
  for (a in seq(unlist(covariates$outcomes))){
    outcome=outcomes[a]

    add_covars=c(as.character(unlist(covariates$covariates[covariates$outcomes==outcome])))

    #cycling through exposures
    for (h in 1:length(exposures)){

      exposure=exposures[h]

      #calculates the mean value for each exposure for each exposure epoch and apends to data
      for (e in 1:nrow(exposure_epochs)){
        epoch=exposure_epochs[e,1]
        temp=data.frame(row.names=1:nrow(data))
        new_var=paste0(exposure, "_", epoch)
        #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
        for (l in 1:length(as.numeric(unlist(exposure_epochs[e,2])))){
          level=as.numeric(unlist(exposure_epochs[e,2]))[l]
          # z=data[,which(grepl(paste0(exposure, "_", as.character(level)), names(data)) & grepl())]
          z=data[,which(grepl(exposure, sapply(strsplit(names(data), "_"), "[", 1)) & (level == suppressWarnings(as.numeric(sapply(strsplit(names(data), "_"), "[", 2)))))]
          temp=cbind(temp, z)
        }
        #adds a new variable of the exposure averaged within epoch
        data=data%>%dplyr::mutate(!!new_var :=rowMeans(temp, na.rm=T))
      }

      #from truncateWweights
      # weights=paste0(exposure, "-", outcome, "_weight_cutoff") #uncomment after testing!!!
      weights=paste0(exposure, "_", outcome)

      #finding all possible exposure histories
      exposure_levels=apply(gtools::permutations(2, nrow(exposure_epochs), c("l", "h"), repeats.allowed=TRUE), 1, paste, sep="", collapse="-")

      #gathering epoch information for each exposure --finds high and low values for each epoch based on user-specified quantile values
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


      #df parameterizing mean values for each epoch (main effects) --i technically calculated these already and appended to dataset
      test=data.frame(names=paste(substr(exposure_epochs$epochs,1,1), "es", sep="_"))
      for (v in 1:nrow(exposure_epochs)){
        vals=as.character(apply(expand.grid(exposure, as.character(as.numeric(unlist(exposure_epochs$values[v])))), 1, paste0, sep="", collapse="_"))
        vals=noquote(paste0("mean(", paste(vals, sep=" ", collapse=" "), ");"))
        vals=paste(test[v,1], "=", vals)
        test[v,2]=vals
      }

      #listing out main effects of each epoch
      main_effects_terms=test$names

      #gathering all variables --ID, time-varying exposures, wave, covariates, outcomes --should these just be all variables in dataset?
      # variables=c(apply(expand.grid(exposures, as.character(as.numeric(unlist(exposure_epochs$values)))), 1, paste0, sep="", collapse=""), outcomes)
      # variables=as.data.frame(c(ID, variables[order(variables)], covariates, "WAVE")) #add others here?
      # variables[nrow(variables),1]=paste0(variables[nrow(variables),1], ";")
      variables=as.data.frame(colnames(data))
      variables=rbind(variables, ";")
      colnames(variables)= "V1"

      #listing out just mean effects in one string
      epoch_info2=as.data.frame(test$V2)
      colnames(epoch_info2)= "V1"

      #specifying all 2-way interactions between epoch main effects
      temp=as.data.frame(cbind(apply(substr(combn(test$names, 2),1,1), 2, paste, sep="", collapse="_"),
                               apply(combn(test$names, 2), 2, paste, sep="", collapse="*"),
                               rep(";")))
      two_way_ints=as.data.frame(paste0(temp[,1]," = ", temp[,2], temp[,3]))
      colnames(two_way_ints)= "V1"
      two_way_int_terms=temp$V1

      #specifying 3-way interaction for main effects; changed label to "e_m_l" for consistency
      three_way_ints=paste0(paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_"), " = ", paste(test$names, sep="", collapse="*"), ";") #adds interaction
      three_way_int_terms=paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_")

      #variables to use: main effect terms, 2-way ints, 3-way ints, wave, outcomes, covariates --anything else?
      use_var=as.data.frame(c(outcome, unlist(as.character(add_covars)), "WAVE",
                              test$names, "q", paste0(temp[,1]),paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_"))) #not quite sure all that goes in this field
      colnames(use_var)= "V1"

      #between-subject variables
      between=as.data.frame(c(test$names, paste0(temp[,1]), paste(substr(exposure_epochs$epochs,1,1), sep="", collapse="_"), add_covars)) #what else goes here?
      colnames(between)= "V1"

      #specifying outcome model (between-person)
      outcome_on=as.data.frame(paste(as.character(between[1:(nrow(between)-length(add_covars)),1])))
      outcome_on[,2]=paste0("(", apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""), ")")
      outcome_on[,3]=paste(outcome_on[,1], outcome_on[,2])
      outcome_on_final=as.data.frame(c((outcome_on[,3]), unlist(add_covars)))
      colnames(outcome_on_final)= "V1"

      #specifying model for the slopes (between-person)
      s_on=as.data.frame(paste(as.character(between[1:(nrow(between)-length(add_covars)),1])))
      s_on[,2]=paste0("(", apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""),"a",")")
      s_on[,3]=paste(s_on[,1], s_on[,2])
      s_on_final=as.data.frame(c((s_on[,3])))
      colnames(s_on_final)= "V1"

      #specifying model for the quadratic effect (between-person)
      qd_on=as.data.frame(paste(as.character(between[1:(nrow(between)-length(add_covars)),1])))
      qd_on[,2]=paste0("(", apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""),"b",")")
      qd_on[,3]=paste(qd_on[,1], qd_on[,2])
      qd_on_final=as.data.frame(c((qd_on[,3])))
      colnames(qd_on_final)= "V1"

      #listing out constraints for outcome, slope, quadratic
      histories=toupper(exposure_levels)
      histories=noquote(gsub("-", "", histories))
      constraints=noquote(paste0(paste("new ("), paste0(histories, sep="", collapse=" "), " ",
                                 paste0(noquote(apply(expand.grid(c("s", "q"), histories), 1, paste0, sep="", collapse="_")[order(apply(expand.grid(c("s", "q"), histories), 1, paste0, sep="", collapse="_"))]),
                                        sep="", collapse=" "), " ",
                                 paste0(noquote(apply(expand.grid(c("d", "sd", "qd"), histories[! histories %in% reference]) #those with no reference
                                                      , 1, paste0, sep="", collapse="_")[order(apply(expand.grid(c("d","sd", "qd"), histories[! histories %in% reference]),
                                                                                                     1, paste0, sep="", collapse="_"))]), sep="", collapse=" "),
                                 ");", sep=" ", collapse=" "))



      #making contrasts
      #finds reference beta values and formula based on parameters of best-fitting model
      ref_parameters=setNames(data.frame(matrix(ncol = length(histories)+2, nrow = 7)), c("sequence", "parameter", histories))
      ref_parameters$parameter= as.data.frame(unlist((apply(expand.grid("p", 2:(2*nrow(epoch_info)+2)),1, paste0, sep="", collapse=""))))
      ref_parameters$sequence=outcome_on[,1]
      ref_parameters=as.data.frame(ref_parameters)

      # browser()

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
        form=noquote(paste0(histories[x], " = p1+", (paste0("(", as.character(unlist(ref_parameters[,2])), "*(",
                                                            as.character(unlist(ref_parameters[,colnames(ref_parameters)[colnames(ref_parameters)==histories[x]]])), "))",
                                                            sep="", collapse=" + ")), ";", sep="", collapse=""))
        forms[x,1]=histories[x]
        forms[x,2]=form

        #slope model --not sure why this only includes p1-p4? will this always be the case?
        forms[x,3]=paste0("s_", histories[x])
        form_s=noquote(paste0(paste0("s_", histories[x]), " = p1a+", (paste0("(", as.character(unlist(ref_parameters[1:3,2])), "a*(",
                                                                             as.character(unlist(ref_parameters[,colnames(ref_parameters)[colnames(ref_parameters)==histories[x]]][1:3])), "))",
                                                                             sep="", collapse=" + ")), ";", sep="", collapse=""))
        forms[x,4]=form_s

        #quadratic model--not sure why this only includes p1-p4?
        forms[x,5]=paste0("q_", histories[x])
        form_q=noquote(paste0(paste0("q_", histories[x]), " = p1b+", (paste0("(", as.character(unlist(ref_parameters[1:3,2])), "b*(",
                                                                             as.character(unlist(ref_parameters[,colnames(ref_parameters)[colnames(ref_parameters)==histories[x]]][1:3])), "))",
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

      #write out new data file with no column names
      data2=data
      data2[is.na(data2)] <- missing #switching back to old missing value for mplus

      write.table(as.data.frame(data2), col.names = F, row.names=F, sep=",", paste0(home_dir, "for MPlus/", exposure, "_", outcome, "_for_Mplus.csv"))


      #creating the .inp text for the Mplus file
      model_input <- as.data.frame(NULL)
      model_input[1,1] <- paste0("TITLE: ", exposure, " ", outcome)

      model_input[3,1] <- paste(" ")

      model_input[4,1] <- paste0("   DATA: file is ", paste0(home_dir, "for MPlus/", exposure, "_", outcome, "_for_Mplus.csv;"))
      model_input=rbind(model_input, empty_line)

      model_input[6,1] <- paste("   define:")
      model_input[7,1] <- paste("   Q = WAVE*WAVE;")
      model_input[8,1] <- paste("   C = WAVE*WAVE*WAVE;")
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, epoch_info2)
      model_input=rbind(model_input, empty_line)
      model_input=rbind(model_input, two_way_ints)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, three_way_ints)
      model_input=rbind(model_input, empty_line)

      model_input= rbind(model_input, paste("   VARIABLE:     NAMES ="))
      model_input= rbind(model_input, variables)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste0("   MISSING ARE ALL (-9999);"))
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, paste0("  USEVAR ="))
      model_input=rbind(model_input, use_var)
      model_input=rbind(model_input, empty_line)

      model_input=rbind(model_input, ";")
      model_input=rbind(model_input, paste0("   within= WAVE Q;")) #will time_var always be the same?
      model_input=rbind(model_input, paste("   between ="))
      # model_input=rbind(model_input, c(main_effects_terms, two_way_int_terms, three_way_int_terms, add_covars))
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
      model_input=rbind(model_input, paste0("S |", outcome, " on WAVE;" )) #will time_var always be the same?
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
