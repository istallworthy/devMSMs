
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

getModel <- function(d, exposure, exposure_time_pts, outcome, epochs, exp_epochs, 
                     int_order, model, fam, covariates, verbose) {

  if (sum(duplicated(d$"ID")) > 0) {
    stop ("Please provide wide data with a single row per ID.",
         call. = FALSE)
  }

  # exposure epochs

  if ( is.null(epochs)) { #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }

  else { #add epochs by averaging exposure time points

    #adds exposure epochs
    #calculates the mean value for each exposure epoch

    for (e in seq_len(nrow(epochs))) {
      epoch <- epochs[e, 1]
      temp <- data.frame(row.names = seq_len(nrow(d)))
      new_var <- paste0(exposure, ".", epoch)

      if (! new_var %in% colnames(d)) {

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


  # Covariate models checking

  if (model == "m1" | model == "m3" | model == "covs") {

    if (sum(covariates %in% colnames(d)) < length(covariates)) {
      stop ("Please only include covariates that correspond to variables in the wide dataset.",
           call. = FALSE)
    }
    covariate_list <- paste(as.character(covariates), sep = "", 
                            collapse = " + ")

  } else {
    covariate_list <- NULL
  }


  # interaction model checking
  
  if (model == "m2" | model == "m3") {

    if (int_order > nrow(epochs)) {
      stop ("Please provide an interaction order equal to or less than the total number of epochs/time points.",
           call. = FALSE)
    }

    interactions <- paste(
      lapply(2:int_order, function(z) {
        paste(apply(combn(exp_epochs, z), 2, paste, sep = "", collapse = ":"),
              sep = "", collapse = " + ")
      }),
      collapse = " + "
    )

    #create interactions in data

    for (x in seq_len(length(unlist(strsplit(interactions, "\\+"))))) {
      name <- gsub(" ", "", unlist(strsplit(interactions, "\\+"))[x])

      if (! name %in% colnames(d)) {
        temp <- d[, c(gsub(" ", "", as.character(unlist(strsplit(unlist(strsplit(interactions, "\\+"))[x], "\\:"))))) ]
        new <- matrixStats::rowProds(as.matrix(temp), na.rm = TRUE)
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
  else if (model == "covs") {
    fc <- paste(outcome, "~", covariate_list)
    mc <- survey::svyglm(as.formula(fc),
                         family = fam,
                         design = s)
    return(mc)
  }


  # Fitting baseline model w/ main effects only (m0) for all models

  f0 <- paste(outcome, "~", paste0(exp_epochs, sep = "", collapse = " + "))
  m0 <- survey::svyglm(as.formula(f0),
                       family = fam,
                       design = s)

  if (model == "m0") {
    return(m0) # Save model
  }
  else {

    # Baseline + sig covar model OR baseline + sig covar + int model

    if (model == "m1" | model == "m3") {
      f1 <- paste(f0, "+", covariate_list) # Baseline + covariate model
      m1 <- survey::svyglm(as.formula(f1),
                           family = fam,
                           design = s)

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

        m2 <- survey::svyglm(as.formula(f2),
                             family = fam,
                             design = s)

        # Baseline + interactions

        return(m2)
      }

      # Baseline + covars + interactions

      if (model == "m3") {

        # Fitting m3
        
        f3 <- paste0(f1, "+", paste(interactions, sep = "", collapse = " + "))
        f3 <- as.formula(f3)
        m3 <- survey::svyglm(f3,
                             family = fam,
                             design = s)

        # Baseline + covars+ interactions

        return(m3)
      }
    }
  }
}
