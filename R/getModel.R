
#' Fits outcome model
#'
#' @param d wide data frame
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param epochs data frame of exposure epoch labels and values
#' @param exp_epochs
#' @param int_order integer specification of highest order exposure main effects interaction for interaction models
#' @param model
#' @param fam function specification for svyglm model
#' @param covariates list of characters reflecting variable names of covariates for covariate models
#' @param interactions
#' @param verbose TRUE or FALSE indicator for user output
#'
#' @return list of fitted model(s)
#' @export
#'
#' @examples
getModel <- function(d, outcome, epochs, exp_epochs, int_order, model, fam, covariates, interactions, verbose){


  if (sum(duplicated(d$"ID")) > 0){
    stop("Please provide wide data with a single row per ID.", call. = FALSE)
  }

  # exposure epochs
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }
  else { #add epochs by averaging exposure time points
    #adds exposure epochs
    #calculates the mean value for each exposure for each exposure epoch
    for (e in seq_len(nrow(epochs))){
      epoch <- epochs[e, 1]
      temp <- data.frame(row.names = seq_len(nrow(d)))
      new_var <- paste0(exposure, "_", epoch)
      if (! new_var %in% colnames(d)){
        #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
        for (l in seq_len(length(as.numeric(unlist(epochs[e, 2]))))){
          level <- as.numeric(unlist(epochs[e, 2]))[l]
          z  <- d[,which(grepl(paste0(exposure, ".", level), names(d)))]
          temp <- cbind(temp, z)
        }
        #adds a new variable of the exposure averaged within epoch
        d  <- d %>% dplyr::mutate(!!new_var := rowMeans(temp, na.rm=T))
        d[,new_var] <- as.numeric(d[,new_var])
      }
    }
  }


  # Covariate models checking
  if (model == "m1" | model == "m3" | model == "covs") {
    if (sum(covariates %in% colnames(d)) < length(covariates)){
      stop("Please only include covariates that correspond to variables in the wide dataset.", call. = FALSE)
    }
    covariate_list <- paste(as.character(covariates), sep = "", collapse = " + ")
  } else{
    covariate_list <- NULL
  }


  # interaction model checking
  if (model == "m2" | model == "m3"){
    if (int_order > nrow(epochs)){
      stop("Please provide an interaction order equal to or less than the total number of epochs/time points.", call. = FALSE)
    }
    interactions <- paste(
      lapply(2:int_order, function(z) {
        paste(apply(combn(exp_epochs, z), 2, paste, sep = "", collapse = ":"), sep = "", collapse = " + ")
      }),
      collapse = " + "
    )
    #create interactions in data
    for (x in seq_len(length(unlist(strsplit(interactions, "\\+"))))) {
      name <- gsub(" ", "", unlist(strsplit(interactions, "\\+"))[x])
      if (! name %in% colnames(d)){
        temp <- d[, c(gsub(" ", "", as.character(unlist(strsplit(unlist(strsplit(interactions, "\\+"))[x], "\\:"))))) ]
        d  <- d %>% dplyr::mutate(!! name := matrixStats::rowProds(as.matrix(temp), na.rm=T))
      }
    }
  }
  else {
    interactions <- NULL
  }


  s <- survey::svydesign(
    id = ~1,
    data = d,
    weights = ~ weights
  )


  #Null models
  # Fitting intercept-only model
  if (model == "int"){
    fi <- paste(outcome, "~ 1")
    mi <- survey::svyglm(as.formula(fi), family = fam,  design = s) # List of model fitted to all imputed datasets
    return(mi)
  }
  else if (model == "covs"){
    fc <- paste(outcome, "~", covariate_list)
    mc <- survey::svyglm(as.formula(fc), family = fam,  design = s) # List of model fitted to all imputed datasets
    return(mc)
  }


  # Fitting baseline model w/ main effects only (m0) for all models
  f0 <- paste(outcome, "~", paste0(exp_epochs, sep = "", collapse = " + "))
  m0 <- survey::svyglm(as.formula(f0), family = fam,  design = s) # List of model fitted to all imputed datasets

  if (model == "m0") {
    return(m0) # Save model
  }
  else {
    # Baseline + sig covar model OR baseline + sig covar + int model
    if (model == "m1" | model == "m3") {
      f1 <- paste(f0, "+", covariate_list) # Baseline + covariate model
      m1 <- survey::svyglm(as.formula(f1), family = fam, design = s)

      # Baseline + imbalanced covars
      if (model == "m1") {
        return(m1)
      }
    }

    # Baseline + interaction OR baseline + covars + interactions
    if (model == "m2" | model == "m3") {
      f2 <- paste(f0, "+", paste(interactions, sep = "", collapse = " + "))

      # Baseline + interactions
      if (model == "m2") {
        # Fitting m2
        m2 <- survey::svyglm(as.formula(f2), family = fam, design = s)
        # Baseline + interactions
        return(m2)
      }

      # Baseline + covars + interactions
      if (model == "m3") {
        # Fitting m3
        f3 <- paste0(f1, "+", paste(interactions, sep = "", collapse = " + "))
        f3 <- as.formula(f3)
        m3 <- survey::svyglm(f3, family = fam, design = s)
        # Baseline + covars+ interactions

        return(m3)
      }
    }
  }
}
