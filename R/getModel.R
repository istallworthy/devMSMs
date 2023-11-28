
#' Fits outcome model
#'
#' @param d wide data frame
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param epochs data frame of exposure epoch labels and values
#' @param exp_epochs list of exposure epochs
#' @param int_order integer specification of highest order exposure main effects
#'   interaction for interaction models
#' @param model character indicating one of the following outcome models:
#'  * "m0" (exposure main effects)
#'  * "m1" (exposure main effects & covariates)
#'  * "m2" (exposure main effects & their interactions)
#'  * "m3" (exposure main effects, their interactions, & covariates)
#'  * "covs" (covariate only model)
#'  * "int" (intercept only model)
#' @param fam function specification for svyglm model
#' @param covariates list of characters reflecting variable names of covariates
#'   for covariate models
#' @param verbose TRUE or FALSE indicator for user output (default is TRUE)
#' @return list of fitted model(s)
#' @export
#' @examples
#' test <- data.frame(ID = 1:50,
#'                    A.1 = rnorm(n = 50),
#'                    A.2 = rnorm(n = 50),
#'                    A.3 = rnorm(n = 50),
#'                    B.1 = rnorm(n = 50),
#'                    B.2 = rnorm(n = 50),
#'                    B.3 = rnorm(n = 50),
#'                    C = rnorm(n = 50),
#'                    D.3 = rnorm(n = 50))
#' test[, c("A.1", "A.2", "A.3")] <- lapply(test[, c("A.1", "A.2", "A.3")], as.numeric)
#'
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
#'
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
#'                     values = I(list(c(1, 2), c(3))))
#' e <- apply(expand.grid("A", as.character(epochs[, 1])), 1,
#'            paste, sep = "", collapse = ".")
#' test$weights <- w[[1]]$weights
#'
#' g <- getModel(d = test,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               epochs = epochs,
#'               exp_epochs = e,
#'               fam = gaussian,
#'               model = "m0")

getModel <- function(d, exposure, exposure_time_pts, outcome, exp_epochs, 
                     int_order, model, fam, covariates, verbose, epochs = NULL) {
  
  if (any(duplicated(d[["ID"]]))) {
    stop("Please provide wide data with a single row per ID.",
         call. = FALSE)
  }
  
  # exposure epochs
  
  if (is.null(epochs)) { #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }
  else { #add epochs by averaging exposure time points
    
    #checks epoch time points
    if (!is.numeric(unlist(epochs$values))|| 
        any(!as.numeric(unlist(epochs$values)) %in% exposure_time_pts)) {
      stop("Please supply epochs with numeric values that are included in the exposure time points.",
           call. = FALSE)
    }
    
    
    #adds exposure epochs
    #calculates the mean value for each exposure epoch
    
    for (e in seq_len(nrow(epochs))) {
      epoch <- epochs[e, 1]
      temp <- data.frame(row.names = seq_len(nrow(d)))
      new_var <- paste0(exposure, ".", epoch)
      
      if (!new_var %in% colnames(d)) {
        
        #finds data from each time point in each epoch, horizontally aligns all exposure values within the epoch for averaging
        
        for (l in seq_len(length(as.numeric(unlist(epochs[e, 2]))))) {
          level <- as.numeric(unlist(epochs[e, 2]))[l]
          z <- d[, names(d)[grepl(exposure, names(d))]] #finds exposure vars
          cols <- colnames(z)[as.logical(sapply(strsplit(names(z), "\\."), "[", 2) 
                                         == as.character(level))]
          cols <- cols[!is.na(cols)]
          z <- as.numeric(as.character(unlist(z[, cols])))
          temp <- cbind(temp, z)
        }
        
        #adds a new variable of the exposure averaged within epoch
        
        x <- as.data.frame(rowMeans(temp, na.rm = TRUE))
        colnames(x) <- c(new_var)
        d <- cbind(d, x)
        d[, new_var] <- as.numeric(d[, new_var])
      }
    }
  }
  
  
  #adding in any covars interactions to d --unfortunately, even though svyglm doesn't need this per se,
  # i think we have to do this in order for these interaction variables to show up in the data field 
  # of the output from this function (which we will need in compareHistories, specifically in avg_predictions()) --otherwise, int variable names don't seem to show up 
  #split factors
  factor_covariates <- names(d)[sapply(d, is.factor)]
  factor_covariates <- setdiff(factor_covariates, "ID")
  
  if (length(factor_covariates) > 0) {
    d <- cobalt::splitfactor(d, factor_covariates, drop.first = "if2")
    
    factors_split <- names(d)[sapply(strsplit(names(d), "\\_"), "[", 1) 
                              %in% factor_covariates]
  }
  
  # if (!missing(covariates)) {
  #   if (any(grepl("\\:", covariates))) {
  #     ints <- covariates[grepl("\\:", covariates)]
  #     
  #     #making interactions w/ split factors 
  #     
  #     for (x in seq_len(length(ints))) {
  #       vars <- as.character(unlist(strsplit(ints[x], "\\:")))
  #       num_comp <- length(vars)
  #       
  #       f_vars <- NULL
  #       if (any(vars %in% factor_covariates)) {
  #         vars <- do.call(c, lapply(vars, function(y) {
  #           if (y %in% factor_covariates) {
  #             f_vars <- factors_split[sapply(strsplit(factors_split, "\\_"), "[", 1) %in% y]
  #             y <- f_vars } 
  #           y
  #         }))
  #       }
  #       
  #       if (any(as.logical(unlist(lapply(vars, function(x) {
  #         any(!x %in% names(d))}))))) {
  #         stop("Please only include covariate interactions between variables in your data",
  #              call. = FALSE)
  #       }
  #       
  #       ints2 <- combn(vars, num_comp)
  #       ints2 <- as.data.frame(ints2[, sapply(strsplit(ints2[1, ], "\\_"), "[", 1) != 
  #                                      sapply(strsplit(ints2[2, ], "\\_"), "[", 1)])
  #       ints2 <- unlist(lapply(1:ncol(ints2), 
  #                              function(y) {paste(ints2[, y], collapse = ":")} ))
  #       ints2 <- ints2[!duplicated(ints2)]
  #       
  #       prods <- lapply(ints2, function(z) {
  #         v <- as.character(unlist(strsplit(z, "\\:")))
  #         temp <- as.data.frame(d[, v])
  #         prod <- apply(as.matrix(temp), 1, prod)
  #         prod
  #       })
  #       prods <- do.call(rbind.data.frame, prods)
  #       prods <- as.data.frame(t(prods))
  #       names(prods) <- ints2
  #       
  #       #make factor class if both components are factors
  #       for (f in seq_len(length(ints2))) {
  #         vars <- as.character(unlist(strsplit(ints2[f], "\\:")))
  #         if (all(vars %in% factor_covariates)) {
  #           prods[, names(prods)[any(as.logical(unlist(lapply(names(prods), function(k) { 
  #             as.character(unlist(strsplit(k, "\\:"))) %in% f_vars}))))]] <- 
  #             as.data.frame(lapply(prods[, names(prods)[any(as.logical(unlist(lapply(names(prods),function(l) {
  #               as.character(unlist(strsplit(l, "\\:"))) %in% f_vars}))))]], 
  #               as.factor))
  #         }
  #       }
  #       #adding to dataset
  #       
  #       d <- cbind(d, prods)
  #     }
  #   }
  #   
  #   covariates <- c(covariates[!grepl("\\:", covariates)], 
  #                   names(d)[grepl("\\:", names(d))])
  # }
  
  
  # Covariate models checking
  
  if (model %in% c("m1", "m3", "covs")) {
    
    if (any(grepl("\\.", covariates))) {
      cov <- as.character(unlist(strsplit(covariates, "\\:")))
      tv_cov <- cov[grepl("\\.", cov)]
      if (any(as.numeric(gsub("_.*", "", sub(".*\\.(.)", "\\1", 
                                             as.character(unlist(strsplit(tv_cov, "\\:")))))) > 
              exposure_time_pts[1])) {
        warning("Please only include covariates (including interaction instituents) that are time invariant or measured only at the first exposure time point.",
                call. = FALSE)
        cat("\n")
      }
    }
    
    # if (!all(covariates[!grepl("\\:", covariates)] %in% colnames(d))) {
    #   stop("Please only include covariates that correspond to variables in the wide dataset.",
    #        call. = FALSE)
    # }
    
    covariate_list <- paste(c(as.character(covariates)), sep = "", 
                            collapse = " + ")
    
  }
  else {
    covariate_list <- NULL
  }
  
  
  # interaction model checking
  
  if (model %in% c("m2", "m3")) {
    
    if (int_order > nrow(epochs)) {
      stop("Please provide an interaction order equal to or less than the total number of epochs/time points.",
           call. = FALSE)
    }
    
    interactions <- paste(
      lapply(2:int_order, function(z) {
        paste(apply(combn(exp_epochs, z), 2, paste, sep = "", collapse = ":"),
              sep = "", collapse = " + ")
      }),
      collapse = " + "
    )
    
    #create exposure main effect interactions in data
    
    for (x in seq_along(unlist(strsplit(interactions, "\\+")))) {
      name <- gsub(" ", "", unlist(strsplit(interactions, "\\+"))[x])
      
      if (!name %in% colnames(d)) {
        temp <- d[, c(gsub(" ", "", 
                           as.character(unlist(strsplit(unlist(strsplit(interactions, 
                                                                        "\\+"))[x], 
                                                        "\\:"))))) ]
        new <- apply(as.matrix(temp), 1, prod, na.rm = TRUE)
        names(new) <- name
        d <- cbind(d, new)
      }
    }
  }
  else {
    interactions <- NULL
  }
  
  s <- survey::svydesign(
    id = ~ 1,
    data = d,
    weights = ~ weights
  )
  
  #Null models for omnibus testing
  # Fitting intercept-only model
  
  if (model == "int") {
    fi <- paste(outcome, "~ 1")
    mi <- survey::svyglm(as.formula(fi),
                         family = fam,
                         design = s)
    return(mi)
  }
  
  if (model == "covs") {
    fc <- paste(outcome, "~", covariate_list)
    mc <- survey::svyglm(as.formula(fc),
                         family = fam,
                         design = s)
    return(mc)
  }
  
  
  # Fitting baseline model w/ main effects only (m0) for all models
  
  f <- paste(outcome, "~", paste0(exp_epochs, sep = "", collapse = " + "))
  
  # Baseline + sig covar model OR baseline + sig covar + int model
  
  if (model == "m1") {
    f <- paste(f, "+", covariate_list) # Baseline + covariate model
  }
  else if (model == "m2") {
    
    # Baseline + interactions
    
    f <- paste(f, "+", paste(interactions, sep = "", collapse = " + "))
  }
  else if (model == "m3") {
    
    # Baseline + covars + interactions
    
    f <- paste(f, "+", covariate_list) # Baseline + covariate model
    f <- paste(f, "+", paste(interactions, sep = "", collapse = " + "))
  }
  
  survey::svyglm(as.formula(f),
                 family = fam,
                 design = s)
}
