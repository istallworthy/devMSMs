#' Calculate balance stats based on Jackson paper
#'
#' Calculate weighted or unweighted standardized balance statistics for a given
#' exposure time point, using all relevant confounders. Draws on Jackson, 2016
#' approaches to assessing balance for time-varying exposures by weighting
#' statistics based on sample distribution in exposure histories.
#'
#' @param home_dir (optional) path to home directory (required if save.out =
#'   TRUE)
#' @param data data in wide format as: a data frame, path to folder of imputed
#'   .csv files, or mids object
#' @param formulas list of balancing formulas at each time point output from
#'   createFormulas()
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param balance_thresh (optional) one or two numbers between 0 and 1
#'   indicating a single balancingn threshold or thresholds for more and less
#'   important confounders, respectively
#' @param k (optional) imputation number
#' @param weights (optional) list of IPTW weights output from createWeights
#' @param imp_conf (optional) list of variable names reflecting important
#'   confounders (required if two balance thresholds are provided)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally
#' @return data frame of balance statistics
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
#' c <- calcBalStats(data = test,
#'                   formulas = f,
#'                   exposure = "A",
#'                   exposure_time_pts = c(1, 2, 3),
#'                   outcome = "D.3",
#'                   balance_thresh = 0.1,
#'                   save.out = FALSE)
#' c <- calcBalStats(data = test,
#'                   formulas = f,
#'                   exposure = "A",
#'                   exposure_time_pts = c(1, 2, 3),
#'                   outcome = "D.3",
#'                   balance_thresh = c(0.05, 0.1),
#'                   imp_conf = "B2",
#'                   save.out = FALSE)
#' c <- calcBalStats(data = test,
#'                   formulas = f,
#'                   exposure = "A",
#'                   exposure_time_pts = c(1, 2, 3),
#'                   outcome = "D.3",
#'                   balance_thresh = 0.1,
#'                   weights = w[[1]],
#'                   save.out = FALSE)

calcBalStats <- function(data, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = 0, weights = NULL,
                         imp_conf = NULL, home_dir = NULL, verbose = TRUE, save.out = TRUE) {
  
  if (!is.list(formulas) | is.data.frame(formulas)) {
    stop("Please provide a list of formulas for each exposure time point",
         call. = FALSE)
  }
  if (!is.null(weights) && !inherits(weights, "weightitMSM")) {
    stop("Please supply a list of weights output from the createWeights function (via WeightIt::WeightItMSM).",
         call. = FALSE)
  }
  
  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.",
           call. = FALSE)
    }
    if (!is.character(home_dir)) {
      stop("Please provide a valid home directory path as a string if you wish to save output locally.",
           call. = FALSE)
    }
    if (!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.",
           call. = FALSE)
    }
  }
  
  if (!is.logical(verbose)) {
    stop("Please set verbose to either TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(verbose) != 1) {
    stop("Please provide a single TRUE or FALSE value to verbose.",
         call. = FALSE)
  }
  
  if (verbose) {
    rlang::check_installed("knitr")
  }
  
  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[", 1)
  
  exposure_name1 <- paste0(exposure, ".", exposure_time_pts[1])
  
  
  #checking exposure type
  
  if (any(as.logical(unlist(lapply(data[, paste0(exposure, '.', exposure_time_pts)], function(x){
    !inherits(x, "numeric") && !inherits(x, "integer")
  }))))) {
    stop("Please provide an exposure in numeric or integer form.",
         call. = FALSE)
  }
  else if (any(as.logical(unlist(lapply(data[, paste0(exposure, '.', exposure_time_pts)], function(x) {
    inherits(x, "integer") && unique(x) != c(1, 0) } ))))) {
    stop("Please make sure your exposure levels are 1s and 0s for integer exposures.",
         call. = FALSE)
  }
  
  exposure_type <- if (is.numeric(data[[exposure_name1]])) "continuous" else "binary"
  
  
  weighted <- !is.null(weights)
  
  factor_covariates <- names(data)[sapply(data, is.factor)]
  factor_covariates <- setdiff(factor_covariates, "ID")
  
  if (weighted) {
    weights_method <- weights$method
    data$weights <- as.numeric(weights$weights) #IPTW weights
  }
  else {
    weights_method <- "no weights"
  }
  
  folder <- if (weighted) "weighted/" else "prebalance/"
  
  data_type <- if (k == 0) "single" else "imputed"
  
  if (data_type == "imputed" && verbose) {
    cat(paste0("**Imputation ", k, "**"), "\n")
  }
  
  #checking for char vars
  
  if (any(sapply(data, is.character))) {
    stop(sprintf("The following variables are characters. Please convert them to factor variables: %s",
                 paste(names(data[sapply(data, is.character)]), collapse = ", ")),
         call. = FALSE)
  }
  
  #split factors
  
  if (length(factor_covariates) > 0) {
    data2 <- cobalt::splitfactor(data, factor_covariates, drop.first = "if2")
  }
  else {
    data2 <- data
  }
  
  factors_split <- names(data2)[sapply(strsplit(names(data2), "\\_"), "[", 1) 
                                %in% factor_covariates]
  
  
  
  #adding in any covars interactions to data2
  
  all_covars <- lapply(formulas, function(x) {
    covars <- deparse1(x[[3]], collapse = "") # gets covariates
    covar_time <- sapply(strsplit(unlist(strsplit(as.character(covars), "\\+")), 
                                  "\\."), tail, 1)
    covars <- as.character(unlist(strsplit(covars, "\\+")))
    covars <- gsub(" ", "", covars)
  })
  all_covars <- unname(do.call(c, all_covars))
  all_covars <- all_covars[!duplicated(all_covars)]
  
  if (any(grepl("\\:", all_covars))) {
    ints <- all_covars[grepl("\\:", all_covars)]
    
    #making interactions w/ split factors 
    
    for (x in seq_len(length(ints))) {
      vars <- as.character(unlist(strsplit(ints[x], "\\:")))
      num_comp <- length(vars)
      if (any(vars %in% factor_covariates)) {
        vars <- do.call(c, lapply(vars, function(y) {
          if (y %in% factor_covariates) {
            f_vars <- factors_split[sapply(strsplit(factors_split, "\\_"), "[", 1) %in% y]
            y <- f_vars } 
          y
        }))
      }
      
      ints2 <- combn(vars, num_comp)
      ints2 <- as.data.frame(ints2[, sapply(strsplit(ints2[1, ], "\\_"), "[", 1) != 
                                     sapply(strsplit(ints2[2, ], "\\_"), "[", 1)])
      ints2 <- unlist(lapply(1:ncol(ints2), 
                             function(y) {paste(ints2[, y], collapse = ":")} ))
      ints2 <- ints2[!duplicated(ints2)]
      
      prods <- lapply(ints2, function(z) {
        v <- as.character(unlist(strsplit(z, "\\:")))
        temp <- as.data.frame(data2[, v])
        prod <- apply(as.matrix(temp), 1, prod)
        prod
      })
      prods <- do.call(rbind.data.frame, prods)
      prods <- as.data.frame(t(prods))
      names(prods) <- ints2
      
      #make factor class if both components are factors
      for (f in seq_len(length(ints2))) {
        vars <- as.character(unlist(strsplit(ints2[f], "\\:")))
        if (all(vars %in% factor_covariates)) {
          prods[, names(prods)[any(as.logical(unlist(lapply(names(prods), function(k) { 
            as.character(unlist(strsplit(k, "\\:"))) %in% f_vars}))))]] <- 
            as.data.frame(lapply(prods[, names(prods)[any(as.logical(unlist(lapply(names(prods),function(l) {
              as.character(unlist(strsplit(l, "\\:"))) %in% f_vars}))))]], 
              as.factor))
        }
      }
      #adding to dataset
      
      data2 <- cbind(data2, prods)
    }
    
  }
  
  
  
  #creating initial data frames
  #data frame with all sampling weights for all exposures at all exposure time points for all histories
  
  all_prop_weights <- data.frame("ID" = NA,
                                 exposure = NA,
                                 exp_time = NA,
                                 history = NA)
  all_bal_stats <- data.frame()
  
  
  #cycles through user-specified exposure time points
  
  for (z in seq_along(exposure_time_pts)) {
    exposure_time_pt <- exposure_time_pts[z]
    lagged_time_pts <- exposure_time_pts[exposure_time_pts < exposure_time_pt]
    
    # GETS COVARIATES FROM FORM FOR ASSESSING BALANCE
    
    full_form <- formulas[[names(formulas)[as.numeric(sapply(strsplit(names(formulas), 
                                                                      "-"), "[", 2)) == exposure_time_pt]]]
    covars <- deparse1(full_form[[3]], collapse = "") # gets covariates
    covar_time <- sapply(strsplit(unlist(strsplit(as.character(covars), "\\+")), 
                                  "\\."), "[", 2)
    covars <- as.character(unlist(strsplit(covars, "\\+")))
    covars <- gsub(" ", "", covars)
    
    
    #checking covariates from formulas 
    
    if (any(duplicated(covars))) {
      stop(sprintf("The following variable(s) are duplicated in the formula for exposure at time point %s: %s.",
                   exposure_time_pt, paste(covars[duplicated(covars)], 
                                           collapse = ", ")),
           call. = FALSE)
    }
    
    if (!all(covars[!grepl("\\:", covars)] %in% c(names(data)))) {
      stop(sprintf("The following variable(s) included in the balancing formulas are not present in your data: %s,",
                   paste(covars[!grepl("\\:", covars)][! covars %in% names(data)], collapse = ", ")),
           call. = FALSE
      )
    }
    
    
    exposure_name <- paste0(exposure, ".", exposure_time_pt)
    
    
    
    #adding appropriate interaction terms for factors to covar list
    
    if (any(grepl("\\:", covars))) {
      ints <- covars[grepl("\\:", covars)]
      comps <- lapply(ints, function (x) {
        as.character(unlist(strsplit(x, "\\:")))})
      comps <- do.call(c, comps)
      
      if (!all(comps %in% c(names(data)))) {
        stop(sprintf("The following variable(s) included in the interaction terms in the balancing formulas are not present in your data: %s,",
                     paste(comps[! comps %in% names(data)], collapse = ", ")),
             call. = FALSE
        )
      }
      
      factor_ints <- ints[as.logical(unlist(lapply(ints, function(x) {
        vars <- as.character(unlist(strsplit(x, "\\:")))
        any(vars %in% factor_covariates)
      })))]
      
      temp <- names(data2)[grepl("\\:", names(data2))]
      f_ints_full <- lapply(factor_ints, function(x){
        f_vars <- temp[gsub("_.", "", temp) == x]
        f_vars
      })
      f_ints_full <- do.call(c, f_ints_full)
      covars <- c(covars[!covars %in% factor_ints], f_ints_full)
    }
    
    #making factor covars split variables in covar list (ignoring ints)
    
    if (length(factor_covariates) > 0) {
      
      data_cov <- data[, covars[!grepl("\\:", covars)]]
      data_cov <- cobalt::splitfactor(data_cov,
                                      names(data_cov)[sapply(data_cov, class ) == "factor"],
                                      drop.first = "if2" )
      covars <- c(covars[grepl("\\:", covars)], colnames(data_cov))
    }
    
    # GETTING BALANCE STATS FOR T=1 W/ NO EXPOSURE HISTORY (ok to standardize immediately)
    
    if (length(lagged_time_pts) == 0) {
      temp <- data2
      
      # Unweighted pre-balance checking
      
      if (!weighted) {
        if (exposure_type == "continuous") {
          
          bal_stats <- cobalt::col_w_cov(temp[covars], 
                                         temp[[exposure_name]], 
                                         std = TRUE) # finding correlation
        }
        else if (exposure_type == "binary") {
          bal_stats <- cobalt::col_w_smd(temp[covars], 
                                         temp[[exposure_name]], 
                                         std = TRUE) # finding smd
          
        }
      }
      
      # Weighted balance checking
      
      else {
        if (exposure_type == "continuous") {
          
          # finding cor
          
          bal_stats <- cobalt::col_w_cov(temp[covars], 
                                         temp[[exposure_name]], 
                                         std = TRUE,
                                         weights = temp[["weights"]]) #IPTW weights
        }
        else if (exposure_type == "binary") {
          
          # finding smd
          
          bal_stats <- cobalt::col_w_smd(temp[covars], 
                                         temp[[exposure_name]], 
                                         std = TRUE,
                                         weights = temp[["weights"]]) #IPTW weights
          
        }
      }
      bal_stats <- as.data.frame(bal_stats)
      
      
      names(bal_stats) <- "std_bal_stats"
    } #ends lag=0
    
    
    # ASSIGNING HISTORIES FOR EXP TIME POINTS T>1
    
    else {
      
      # creating proportion weights based on proportion of individuals in a given exposure history
      
      prop_weights <- data.frame(ID = data[["ID"]],
                                 exposure = exposure,
                                 exp_time = exposure_time_pt,
                                 history = NA)
      
      # finding histories up until exp time point T
      
      histories <- apply(perm2(length(lagged_time_pts), c(1, 0)),
                         1, paste, sep = "", collapse = "-")
      histories <- as.data.frame(histories)
      
      # cycle through histories to identify individuals with each history
      
      for (h in seq_len(nrow(histories))) {
        his <- as.data.frame(histories[h, ])
        his <- as.character(unlist(his[1, ]))
        
        # get wide data variable names for exposures at each time in each history
        
        exps_time <- apply(expand.grid(exposure, as.character(lagged_time_pts)), 1,
                           paste, collapse = ".", sep = ".")
        
        # initiate flag marking whether each individual falls into the given history at 0 (t(1)-1)
        
        data$flag <- 0
        
        # cycles through times in each history and sequentially flags individuals who meet each criterion
        # (e.g., hi at 6, lo at 15, etc.) for that history
        
        for (t in seq_along(lagged_time_pts)) {
          time <- lagged_time_pts[t]
          exp <- as.numeric(sapply(strsplit(as.character(his), "-"), "[", t)) # hi/lo indicator
          flag <- t - 1
          
          if (exp == 0) { # low levels/absent
            if (exposure_type == "continuous") {
              data$flag <- ifelse(data[[exps_time[t]]] <= median(data[[exposure_name]])
                                  & data$flag == flag, t, NA) # finding those w/ vals <= median exp @ time pt & flagged at prev t's
            }
            else { # for binary exp
              data$flag <- ifelse(data[[exps_time[t]]] == 0 & data$flag == flag, 
                                  t, NA) # if exposure is absent & flagged at prev t's
            }
            
          }
          else { # hi levels/present
            if (exposure_type == "continuous") {
              data$flag <- ifelse(data[[exps_time[t]]] > median(data[[exposure_name]])
                                  & data$flag == flag, t, NA) # finding those w/ vals > median exp @ time pt & flagged at prev t's
            }
            else { # binary exp
              data$flag <- ifelse(data[[exps_time[t]]] == 1 & data$flag == flag, 
                                  t, NA) # if exposure is present & flagged at prev t's
              
            }
            
          }
          
        } # ends history's constituent time pts time points (e.g., 6, 15, 24)
        
        
        # finding ids who met criteria for that history
        
        ids <- data[["ID"]][!is.na(data$flag) & data$flag == t]
        
        # labels those ids w/ that history
        
        prop_weights[["history"]][prop_weights[["ID"]] %in% ids] <- paste(his, collapse = ",")
        data$flag <- NULL # resets flag
      } # ends history loop (e.g., "l-l-l")
      
      prop_sum <- aggregate(exposure ~ as.factor(history), data = prop_weights,
                            FUN = length)
      
      
      # GET BALANCE STATISTICS FOR T>1 (when there is a history to weight on)
      
      if (length(lagged_time_pts) > 0) {
        
        # Merge data with prop_weights by ID
        
        temp <- merge(data2, prop_weights, by = "ID", all.x = TRUE) #data with factors split up
        
        
        # Removing any histories that only have 1 or 0 person contributing (cannot calc bal stats)
        
        if (any(prop_sum$exposure == 1) || any(prop_sum$exposure == 0)) {
          
          omitted_histories <- as.character(as.data.frame(prop_sum)[[1]][prop_sum$exposure == 1 | 
                                                                           prop_sum$exposure == 0])
          
          
          if (data_type == "imputed"){
            
            cat(sprintf("USER ALERT: the following history/histories, %s has/have been omitted from
                        balance checking for exposure %s imputation %s at time point %s due to insufficient counts:",
                        omitted_histories, exposure, k, exposure_time_pt))
            cat("\n")
          }
          
          else {
            
            cat(sprintf("USER ALERT: the following history/histories, %s has/have been omitted from
                        balance checking for exposure %s at time point %s due to insufficient counts:",
                        omitted_histories, exposure, exposure_time_pt))
            cat("\n")
          }
          
          temp <- temp[!temp$history %in% omitted_histories, , drop = FALSE]
        } #ends hist exc
        
        
        # Unweighted pre-balance checking
        
        if (!weighted) { # no IPTW weighting but weighting on history
          
          if (exposure_type == "continuous") {
            bal_stats <- sapply(sort(unique(temp$history)), function(i) { # finding balance by history
              
              temp2 <- temp[temp$history == i, , drop = FALSE]
              
              cobalt::col_w_cov(temp2[covars], 
                                temp2[[exposure_name]], 
                                std = FALSE)
              
            })
            
            # getting weighted mean across histories (cols), weighting by proportion of those w/ that same history
            
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ],
                            tabulate(as.factor(temp$history)) / nrow(temp))
            })
            
            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))
            
            bal_stats$std_bal_stats <- weighted_bal_stats /
              (sapply(rownames(bal_stats), function(x) { 
                sd(as.numeric(data2[[x]]), na.rm = TRUE)
              }) *# unweighted covar sd
                sd(data[[exposure_name]], na.rm = TRUE))  # exposure SD at that time pt
            
            
            bal_stats <- bal_stats[startsWith(names(bal_stats), "std")]
            
          } #ends continuous
          
          else if (exposure_type == "binary") {
            bal_stats <- sapply(sort(unique(temp$history)), function(i) {
              
              temp2 <- temp[temp$history == i, , drop = FALSE ]
              
              cobalt::col_w_smd(temp2[covars], 
                                temp2[[exposure_name]], 
                                std = FALSE) # finding mean difference
              
            })
            
            # getting weighted mean across histories, weighting by proportion of those w/ that same history
            
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ],
                            tabulate(as.factor(temp$history)) / nrow(temp))
            })
            
            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))
            # 
            
            # standardizing balance statistics after finding weighted balance stats
            
            bal_stats$std_bal_stats <- weighted_bal_stats/
              sapply(covars, function(x) {
                sqrt(mean( #dividing by pool SD estimate (unadjusted)
                  var(as.numeric(data[[x]][data[[exposure_name1]] == 1])), #treated var
                  var(as.numeric(data[[x]][data[[exposure_name1]] == 0])) #untreated var
                ))
              })
            
            
            # For a weighted_bal_stat of 0, make std stat also 0 so as not to throw an error
            
            bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] <- 0
            
            bal_stats <- bal_stats[startsWith(names(bal_stats), "std")]
            
          } #ends binary
        } #ends weighted=0
        
        
        # Weighted balance checking
        
        else { # if weighted, use IPTW weights from weightitmsm and weight by history
          
          if (exposure_type == "continuous") {
            
            # finds balance for each covariate clustered/subset by history
            
            bal_stats <- sapply(sort(unique(temp$history)), function(i) {
              
              temp2 <- temp[temp$history == i,, drop = FALSE]
              
              cobalt::col_w_cov(temp2[covars], 
                                temp2[[exposure_name]], 
                                std = FALSE, # finding covariance
                                weights = temp2[["weights"]]) # adding IPTW weights
            })
            
            # getting weighted mean across histories, weighting by proportion of those w/ that same history
            
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ],
                            tabulate(as.factor(temp$history)) / nrow(temp))
            })
            
            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))
            
            # standardizing balance statistics after weighting by history
            
            bal_stats$std_bal_stats <- weighted_bal_stats /
              (sapply(rownames(bal_stats), function(x) { #issue: looking in data for unweighted vals but factors have additional vars
                sd(as.numeric(data2[[x]]), na.rm = TRUE) }) *# unweighted covar sd
                 sd(data[[exposure_name]], na.rm = TRUE))  # exposure SD at that time pt
            
            
            # For a weighted_bal_stat of 0, make std stat also 0 so as not to throw an error
            
            bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] <- 0
            
            bal_stats <- bal_stats[startsWith(names(bal_stats), "std")]
            
          } #ends continuous
          
          else if (exposure_type == "binary") {
            
            # finds balance for each covariate clustered/subset by history
            
            bal_stats <- sapply(sort(unique(temp$history)), function(i) {
              temp2 <- temp[temp$history == i,, drop = FALSE]
              
              cobalt::col_w_smd(temp2[covars], 
                                temp2[[exposure_name]], 
                                std = FALSE, # finding mean difference
                                weights = temp2[["weights"]]) # adding IPTW weights
            })
            
            # getting weighted mean across histories, weighting by proportion of those w/ that same history
            
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ],
                            tabulate(as.factor(temp$history)) / nrow(temp))
            })
            
            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))
            
            # 
            # standardizing balance statistics after finding weighted balance stats
            
            bal_stats$std_bal_stats <- weighted_bal_stats/
              sapply(covars, function(x) {
                sqrt(mean( #dividing by pool SD estimate (unadjusted)
                  var(as.numeric(data[[x]][data[[exposure_name1]] == 1])), #treated var
                  var(as.numeric(data[[x]][data[[exposure_name1]] == 0])) #untreated var
                ))
              })
            
            # For a weighted_bal_stat of 0, make std stat also 0 so as not to throw an error
            
            bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] <- 0
            
            bal_stats <- bal_stats[startsWith(names(bal_stats), "std")]
            
          } #ends binary
        } #ends weighted
        
        # collect proportions for histories at this time point
        
        all_prop_weights <- rbind(all_prop_weights, prop_weights)
      }
    } # ends lag>0 loops
    
    
    
    # ADDS INFO TO BAL STATS
    
    bal_stats <- as.data.frame(bal_stats)
    
    bal_stats$covariate <- rownames(bal_stats)
    
    
    #averages across factor levels to create one bal stat per factor variable?
    covariate <- NULL #to avoid check() binding problem
    data$ID <- as.numeric(data$ID)
    f_vars <- colnames(data)[sapply(data, is.factor)]
    
    if (length(f_vars) > 0) {
      f_vars <- rownames(bal_stats)[sapply(strsplit(rownames(bal_stats), 
                                                    "_"), "[", 1) %in% f_vars]
      f_vars <- f_vars[!grepl("\\:", f_vars)] #excludes ints
      f_stats <- subset(bal_stats, covariate %in% f_vars)
      f_stats$name <- sapply(strsplit(rownames(f_stats), "_"), "[", 1)
      test <- aggregate(std_bal_stats ~ name, data = f_stats,
                        FUN = mean)
      
      colnames(test) <- c("covariate", "std_bal_stats")
      new <- data.frame(std_bal_stats = test$std_bal_stats,
                        covariate = test$covariate)
      rownames(new) <- new$covariate
      bal_stats <- rbind(as.data.frame(subset(bal_stats, !covariate %in% f_vars)), new)
    }
    
    #adds custom bal thresh info
    
    if (!is.null(imp_conf)) {
      bal_stats$bal_thresh <- ifelse(bal_stats$covariate %in% imp_conf,
                                     balance_thresh[1], balance_thresh[2])
      bal_stats$balanced <- ifelse(abs(bal_stats$std_bal_stats) < 
                                     bal_stats$bal_thresh, 1, 0)
    }
    else {
      bal_stats$bal_thresh <- balance_thresh
      bal_stats$balanced <- ifelse(abs(bal_stats$std_bal_stats) < 
                                     bal_stats$bal_thresh, 1, 0)
      
    }
    
    bal_stats$exposure <- exposure
    bal_stats$exp_time <- exposure_time_pt
    bal_stats$covar_time <- NA
    bal_stats$covar_time[grepl("\\.", bal_stats$covariate)] <- 
      as.numeric(sapply(strsplit(bal_stats$covariate[grepl("\\.", bal_stats$covariate)], "\\."), tail, 1))
    
    all_bal_stats <- rbind(all_bal_stats, bal_stats)
    all_bal_stats$covar_time[is.na(all_bal_stats$covar_time)] <- 0
    
    # Make love plot per exposure time point
    
    make_love_plot(home_dir = home_dir, folder = folder, exposure = exposure,
                   exposure_time_pt = exposure_time_pt,
                   exposure_type = exposure_type, 
                   k = k, form_name = form_name, balance_stats = bal_stats,
                   data_type = data_type, balance_thresh = balance_thresh, 
                   weights_method = weights_method, imp_conf = imp_conf,
                   save.out = save.out, verbose = verbose)
    
  }     # Ends exp_time_pt loop
  
  
  if (verbose & save.out) {
    if (data_type == "imputed") {
      
      cat(sprintf("For each time point and imputation, %s summary plots for  %s
                 formulas weighting method %s have now been saved in the %s plots/' folder.\n",
                  gsub("/", "", folder), form_name, weights_method, folder))
      cat("\n")
      
    }
    else {
      
      cat(sprintf("For each time point, %s summary plots for  %s
                 formulas weighting method %s have now been saved in the %s plots/' folder.\n",
                  gsub("/", "", folder), form_name, weights_method, folder))
      cat("\n")
      
    }
  }
  
  
  # Summarizing balance
  
  bal_summary_exp <- as.data.frame(aggregate(balanced ~ exp_time,
                                             data = all_bal_stats,
                                             FUN = function(x) c(balanced_n = sum(x == 1),
                                                                 imbalanced_n = sum(x == 0),
                                                                 n = length(x))))
  bal_summary_exp <- do.call(data.frame, bal_summary_exp) #? IS: b/c aggregate() makes new variables into a single weird matrix thing
  names(bal_summary_exp) <- c("exp_time", "balanced_n", "imbalanced_n", "n")
  
  
  if (save.out) {
    write.csv(bal_summary_exp,
              file.path(home_dir, "balance", folder, 
                        sprintf("%s_%s_%s_%s_balance_stat_summary.csv",
                                form_name, exposure, k, weights_method)))
    
    if (verbose) {
      cat("\n")
      if (data_type == "imputed") {
        
        cat(sprintf("Balance statistics using %s formulas for %s imputation %s, using %s have been saved in the 'balance/%s' folder. \n",
                    form_name, exposure, k, weights_method, folder))
        cat("\n")
        
        cat(sprintf("Sampling weights using the %s for %s imputation %s have been saved in the 'balance/%s' folder., \n",
                    form_name, exposure, k, folder))
        cat("\n")
      }
      else {
        
        cat(sprintf("Balance statistics using %s formulas for %s using %s have been saved in the 'balance/%s' folder. \n",
                    form_name, exposure, weights_method, folder))
        cat("\n")
        
        cat(sprintf("Sampling weights using the %s for %s have been saved in the 'balance/%s' folder., \n",
                    form_name, exposure, folder))
      }
      cat("\n")
    }
  }
  
  
  # tallies total possible COVARIATES FROM FORM FOR ASSESSING BALANCE
  
  all_form <- as.data.frame(do.call(rbind, formulas))
  tot_covars <- deparse1(all_form[, 3])
  tot_covars <- as.character(unlist(strsplit(tot_covars, "\\+")))[
    !grepl("form", as.character(unlist(strsplit(tot_covars, "\\+"))))]
  tot_covars <- gsub(" ", "", tot_covars)
  tot_covars <- na.omit(sapply(strsplit(tot_covars, "\\."), head, 1)[
    !duplicated(sapply(strsplit(tot_covars, "\\."), "[", 1))])
  
  imbalanced_covars <- sum(bal_summary_exp$imbalanced_n, na.rm = TRUE)
  total_covars <- sum(bal_summary_exp$n, na.rm = TRUE)
  total_domains <- length(tot_covars)
  
  if (imbalanced_covars > 0) {
    percentage_imbalanced <- round((imbalanced_covars / total_covars) * 100, 0)
    
    remaining_imbalanced_domains <- length(sapply(strsplit(all_bal_stats[all_bal_stats$balanced == 0, "covariate"], "\\."),
                                                  "[", 1)[!duplicated(sapply(strsplit(all_bal_stats[all_bal_stats$balanced == 0, "covariate"],
                                                                                      "\\."), "[", 1))])
    
    # remaining_avg_abs_corr <- round(mean(abs(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"]), na.rm = TRUE), 2)
    remaining_avg_abs_corr <- round(median(abs(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"]), 
                                           na.rm = TRUE), 2)
    remaining_corr_range <- paste0(round(min(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"], 
                                             na.rm = TRUE), 2),
                                   "-", round(max(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"], 
                                                  na.rm = TRUE), 2))
  }
  
  
  if (verbose) {
    if (data_type == "imputed") {
      
      cat(sprintf("USER ALERT: For exposure %s imputation %s using %s and %s formulas: \n",
                  exposure, k, weights_method, form_name))
      
      cat(sprintf("The median absolute value relation between exposure and confounder is %s (range = %s - %s).\n",
                  round(median(abs(all_bal_stats$std_bal_stats)), 2),
                  round(min(all_bal_stats$std_bal_stats), 2),
                  round(max(all_bal_stats$std_bal_stats), 2)))
      cat("\n")
      
      if (imbalanced_covars > 0) {
        cat("\n")
        cat(sprintf("As shown below, %s out of %s (%s%%) covariates across time points, corresponding to %s out of %s domains, remain imbalanced with a remaining median absolute value correlation/std mean difference of %s (range= %s):\n",
                    imbalanced_covars,
                    total_covars,
                    percentage_imbalanced,
                    remaining_imbalanced_domains,
                    total_domains,
                    remaining_avg_abs_corr,
                    remaining_corr_range))
        
        cat(knitr::kable(bal_summary_exp,
                         caption = sprintf("Imbalanced Covariates for imputation %s using %s and %s formulas",
                                           k, weights_method, form_name),
                         format = 'pipe'),
            sep = "\n")
      }
      else {
        cat(sprintf("No covariates remain imbalanced for imputation %s using %s and %s formulas. \n",
                    k, weights_method, form_name))
      }
      
      cat("\n")
      cat("\n")
      
    }
    else {
      
      if (imbalanced_covars > 0) {
        cat(sprintf("As shown below, %s out of %s (%s%%) covariates across time points, corresponding to %s out of %s domains, remain imbalanced with a remaining median absolute value correlation/std mean difference of %s (range= %s):\n",
                    imbalanced_covars,
                    total_covars,
                    percentage_imbalanced,
                    remaining_imbalanced_domains,
                    total_domains,
                    remaining_avg_abs_corr,
                    remaining_corr_range))
        
        cat("\n")
        cat(knitr::kable(bal_summary_exp,
                         caption = sprintf("Imbalanced covariates using %s and %s formulas", 
                                           weights_method, form_name),
                         format = 'pipe'),
            sep = "\n")
      }
      else {
        cat(sprintf("No covariates remain imbalanced using %s and %s formulas. \n",
                    weights_method, form_name))
      }
      cat("\n\n")
    }
  }
  
  rownames(all_bal_stats) <- NULL
  
  all_bal_stats
}

