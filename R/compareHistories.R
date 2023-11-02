#' Estimate, compare, and visualize exposure histories
#'
#' Takes fitted model output to created predicted values for user-specified
#' histories (pooling for imputed data), before conducting contrast comparisons
#' (pooling for imputed data), correcting for multiple comparisons, and then
#' plotting results.
#' @seealso {[marginaleffects::avg_predictions()],
#'   <https://cran.r-project.org/web/packages/marginaleffects/marginaleffects.pdf>}
#' @seealso {[marginaleffects::hypotheses()],
#'   <https://cran.r-project.org/web/packages/marginaleffects/marginaleffects.pdf>}
#' @seealso {[stats::p.adjust()],
#'   <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust>}
#' @param home_dir path to home directory (required if 'save.out' = TRUE)
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param model list of model outputs from fitModel()
#' @param epochs (optional) data frame of exposure epoch labels and values
#' @param hi_lo_cut (optional) list of two numbers indicating quantile values
#'   that reflect high and low values, respectively, for continuous exposure
#'   (default is median split)
#' @param reference (optional) list sof one or more strings of "-"-separated "l"
#'   and "h" values indicative of a reference exposure history to which to
#'   compare comparison, required if comparison is supplied
#' @param comparison (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is supplied
#' @param mc_comp_method (optional) character abbreviation for multiple
#'   comparison correction method for stats::p.adjust, default is
#'   Benjamini-Hochburg ("BH")
#' @param dose_level (optional) "l" or "h" indicating whether low or high doses
#'   should be tallied in tables and plots (default is high "h")
#' @param exp_lab (optional) character label for exposure variable in plots
#'   (default is variable name)
#' @param out_lab (optional) character label for outcome variable in plots
#'   (default is variable name)
#' @param colors (optional) character specifying Brewer palette or list of
#'   colors (n(epochs)+1) for plotting (default is "Dark2" palette)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return data frame of history comparisons
#' @export
#' @examples
#' f <- createFormulas(exposure = "A",
#'                     exposure_time_pts = c(1, 2, 3),
#'                     outcome = "D.3",
#'                     tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                     ti_confounders = "C",
#'                     type = "full",
#'                     save.out = FALSE)
#'
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
#' w <- createWeights(data = test,
#'                    exposure = "A",
#'                    outcome = "D.3",
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' m <- fitModel(data = test,
#'               weights = w,
#'               exposure = "A",
#'               exposure_time_pts = c(1, 2, 3),
#'               outcome = "D.3",
#'               model = "m0",
#'               save.out = FALSE)
#'
#' r <- compareHistories(exposure = "A",
#'                       exposure_time_pts = c(1, 2, 3),
#'                       outcome = "D.3",
#'                       model = m,
#'                       save.out = FALSE)
#' r <- compareHistories(exposure = "A",
#'                       exposure_time_pts = c(1, 2, 3),
#'                       outcome = "D.3",
#'                       model = m,
#'                       reference = "l-l-l",
#'                       comparison = "h-h-h",
#'                       save.out = FALSE)
#' r <- compareHistories(exposure = "A",
#'                       exposure_time_pts = c(1, 2, 3),
#'                       outcome = "D.3",
#'                       model = m,
#'                       reference = "l-l-l",
#'                       comparison = c("h-h-h", "h-l-l"),
#'                       save.out = FALSE)
#' r <- compareHistories(exposure = "A",
#'                       exposure_time_pts = c(1, 2, 3),
#'                       outcome = "D.3",
#'                       model = m,
#'                       reference = c("l-l-l", "l-h-h"),
#'                       comparison = c("h-h-h"),
#'                       save.out = FALSE)
#' r <- compareHistories(exposure = "A",
#'                       exposure_time_pts = c(1, 2, 3),
#'                       outcome = "D.3",
#'                       model = m,
#'                       reference = c("l-l-l", "l-h-h"),
#'                       comparison = c("h-h-h", "l-l-h"),
#'                       save.out = FALSE)
#' r <- compareHistories(exposure = "A",
#'                       exposure_time_pts = c(1, 2, 3),
#'                       outcome = "D.3",
#'                       model = m,
#'                       reference = c("l-l-l", "l-h-h"),
#'                       comparison = c("h-h-h", "l-l-h"),
#'                       hi_lo_cut = c(0.60, 0.30),
#'                       mc_comp_method = "BH",
#'                       dose_level = "l",
#'                       exp_lab = "Hello",
#'                       out_lab = "Goodbye",
#'                       colors = "Set1",
#'                       save.out = FALSE)
#'                       
compareHistories <- function(home_dir, exposure, exposure_time_pts, outcome, model, epochs = NULL, hi_lo_cut = NULL,
                             reference = NULL, comparison = NULL, mc_comp_method = NA, dose_level = NA, exp_lab = NA, out_lab = NA,
                             colors = NULL, verbose = TRUE, save.out = TRUE ) {
  
  if (missing(exposure)) {
    stop("Please supply a single exposure.",
          call. = FALSE)
  }
  if (!is.character(exposure) || length(exposure) != 1) {
    stop("Please supply a single exposure as a character.",
          call. = FALSE)
  }
  else if (grepl("\\.", exposure)) {
    stop ("Please supply an exposure without the '.time' suffix or any '.' special characters. Note that the exposure variables in your dataset should be labeled with the '.time' suffix.",
          call. = FALSE)
  }
  
  if (missing(outcome)) {
    stop("Please supply a single outcome.",
          call. = FALSE)
  }
  if (!is.character(outcome) || length(outcome) != 1) {
    stop("Please supply a single outcome as a character.",
          call. = FALSE)
  }
  else if (!grepl("\\.", outcome)) {
    stop ("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
          call. = FALSE)
  }
  else if (as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) != 
           exposure_time_pts[length(exposure_time_pts)] && 
           !as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) > 
           exposure_time_pts[length(exposure_time_pts)] ) {
    stop ("Please supply an outcome variable and time point that is equal to or greater than the last exposure time point.",
          call. = FALSE)
  }
  
  if (missing(exposure_time_pts)) {
    stop("Please supply the exposure time points at which you wish to create weights.",
          call. = FALSE)
  }
  if (!is.numeric(exposure_time_pts)) {
    stop("Please supply a list of exposure time points as integers.",
          call. = FALSE)
  }
  else if (!length(exposure_time_pts) > 1) {
    stop ("Please supply at least two exposure time points.",
          call. = FALSE)
  }
  
  if (missing(model)) {
    stop("Please supply a list of model output",
          call. = FALSE)
  }
  if (!is.list(model) || is.data.frame(model)) {
    stop("Please provide a list of model output from the fitModel function.",
          call. = FALSE)
  }
  if (!all(sapply(model, inherits, "svyglm"))) {
    stop("Please supply a model as a list from the createWeights function.",
          call. = FALSE)
  }
  
  if ((!is.character(mc_comp_method) && !is.na(mc_comp_method)) 
      || length(mc_comp_method) != 1) {
    stop("Please provide a single valid character string abbreviation for a multiple comparison method compatible with the stats::p.adjust() function.",
          call. = FALSE)
  }
  if (is.na(mc_comp_method)) {
    mc_comp_method <- "BH"
  }
  
  if ((!is.character(dose_level) && !is.na(dose_level)) || 
      length(dose_level) != 1 ) {
    stop("Please provide a single valid character string of either 'l' or 'h' to indicate whether you with to tally doses of low or high levels of exposure, respectively",
          call. = FALSE)
  }
  if (!dose_level %in% c('h', 'l') && !is.na(dose_level)) {
    stop("Please provide a single valid character string of either 'l' or 'h' to indicate whether you with to tally doses of low or high levels of exposure, respectively",
          call. = FALSE)
  }
  if (is.na(dose_level)) { 
    dose_level <- "h"
  }
  
  if (!is.character(colors) && !is.null(colors)) {
    stop("To specify plotting colors, as character strings, please provide either a list of valid colors equal to the number of exposure main effects + 1 or a valid Brewer palette.",
          .call = FALSE)
  }
  if (is.null(colors)) {
    colors <- "Dark2"
  }
  
  if ((!is.character(exp_lab) && !is.na(exp_lab)) || length(exp_lab) > 1) {
    stop("To specify an alternate exposure label for plotting, please provide a single name as a character string.",
          .call = FALSE)
  }
  
  if ((!is.character(out_lab) && !is.na(out_lab)) || length(out_lab) > 1) {
    stop("To specify an alternate outcome label for plotting, please provide a single name as a character string.",
          .call = FALSE)
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
  
  if (!is.logical(save.out)) {
    stop("Please set save.out to either TRUE or FALSE.",
          call. = FALSE)
  }
  if (length(save.out) != 1) {
    stop("Please provide a single TRUE or FALSE value to save.out.",
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
           .call = FALSE)
    }
    
    rlang::check_installed("stargazer")
    
    histories_dir <- file.path(home_dir, "histories")
    if (!dir.exists(histories_dir)) {
      dir.create(histories_dir)
    }
    plot_dir <- file.path(home_dir, "plots")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
  }
  
  exposure_type <- if (inherits(model[[1]]$data[, paste0(exposure, '.', exposure_time_pts[1])], "numeric"))
    "continuous" else "binary"
  
  
  ints <- gsub(" ", "", 
               as.character(unlist(strsplit(as.character(unlist(model[[1]]$terms)), 
                                            "\\+"))))
  ints <- ifelse(sum(grepl(":", ints)) > 0, 1, 0)
  
  
  #history inspection --done early bc function deals with epochs, etc.
  #print history sample distribution
  
  if (verbose) {
    eval_hist(data = model[[1]]$data, exposure = exposure, epochs = epochs,
              time_pts = exposure_time_pts, hi_lo_cut = hi_lo_cut, 
              ref = reference, comps = comparison)
  }
  
  
  #getting data to use for determining hi/lo values: should have any epochs created/ used in the model
  
  if (is.null(names(model))) { #imputed data
    
    # stacking imputed data to compute hi/lo vals across all imps
    
    data <- lapply(model, `[[`, "data")
    data <- do.call("rbind", data)
  }
  
  else if (!is.null(names(model))) { #single df (not imputed)
    data <- model[[1]]$data
    
  }
  
  if (is.null(epochs)) { #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }
  else {
    if ( !is.data.frame(epochs) || ncol(epochs) != 2 || 
         !all(colnames(epochs) == c("epochs", "values"))) {
      stop("If you supply epochs, please provide a dataframe with two columns of epochs and values.",
            call. = FALSE)
    }
    if (anyNA(epochs$values)) {
      stop("Please provide one or a list of several values for each epoch.",
            call. = FALSE)
    }
  }
  
  terms  <- model[[1]]$formula
  terms <- gsub(" ", "", as.character(unlist(strsplit(as.character(unlist(terms[3])), 
                                                      "\\+"))))
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, 
                      paste, sep = "", collapse = ".")
  
  if (!all(exp_epochs %in% terms)) {
    stop("The use of exposure epochs must be consistent across fitModel and compareHistories, given that the fitted model output is used as the basis of histories.
         If you specified exposure epochs in fitModel, please specify the same ones at this step.
         If you did not, please do not specify them at this step.",
          call. = FALSE)
  }
  
  
  # creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  
  exposure_levels <- apply(perm2(nrow(epochs), c("l", "h")), 1, 
                           paste, sep = "", collapse = "-")
  
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, 
                      paste, sep = "", collapse = ".")
  
  # gathering epoch information for each exposure for deriving betas
  
  eps <- epochs$epochs
  epoch_info <- as.data.frame(rep(exposure, length(eps)))
  epoch_info$time <- eps
  is.na(epoch_info$low) <- TRUE
  is.na(epoch_info$high) <- TRUE
  
  if (exposure_type == "continuous") {
    
    #finds high and low values for each epoch
    
    for (t in seq_len(length(eps))) {
      var_name <- paste(exposure, eps[t], sep = ".")
      
      if (is.null(hi_lo_cut)) { #default is median split
        epoch_info$low[t] <- as.numeric(median(data[, var_name] - 0.001, 
                                               na.rm = T)) #padded with -.001 bc otherwise breaks hypotheses()
        epoch_info$high[t] <- as.numeric(median(data[, var_name] + 0.001 , 
                                                na.rm = T))
        
      }
      else if (!is.numeric(hi_lo_cut)) {
        stop("Please provide one or two numbers to hi_lo_cut between 0-1",
              call. = FALSE)
        
      }
      else {
        if (length(hi_lo_cut) > 2) {
          stop("Please provide either one or two numeric values between 0-1 to hi_lo_cut.",
                call. = FALSE)
          
        }
        if (length(hi_lo_cut) == 2) { #if two cutoff values supplied
          hi_cutoff <- hi_lo_cut[1]
          lo_cutoff <- hi_lo_cut[2]
          
          if (hi_cutoff > 1 || hi_cutoff < 0) {
            stop('Please supply a high cutoff value to hi_lo_cut between 0 and 1',
                  call. = FALSE)
          }
          if (lo_cutoff > 1 || lo_cutoff < 0) {
            stop('Please supply a low cutoff value to hi_lo_cut between 0 and 1',
                  call. = FALSE)
          }
        }
        else { #if only one cutoff value supplied
          if (hi_lo_cut > 1 || hi_lo_cut < 0) {
            stop('Please select a hi_lo cutoff value between 0 and 1',
                  call. = FALSE)
          }
          
          hi_cutoff <- hi_lo_cut
          lo_cutoff <- hi_lo_cut
        }
        
        #find user-specfied quantile values
        
        epoch_info$low[t] <- as.numeric(quantile(data[, var_name], 
                                                 probs = lo_cutoff, na.rm = T))
        epoch_info$high[t] <- as.numeric(quantile(data[, var_name], 
                                                  probs = hi_cutoff, na.rm = T))
      }
    }
  }
  
  
  else if (exposure_type == "binary") {
    for (t in seq_len(length(eps))) {
      var_name <- paste(exposure, eps[t], sep = ".")
      epoch_info$low[t] <- 0
      epoch_info$high[t] <- 1
    }
  }
  
  # gather high and low values for each exposure epoch based on quantile values
  
  d <- data.frame(e = paste(epoch_info[[1]], epoch_info[[2]], sep = "."),
                  l = epoch_info$low,
                  h = epoch_info$high)
  d$v <- paste(d$l, d$h, sep=",")
  d$z <- lapply(seq_len(nrow(d)), function (x) {
    c(as.numeric(unlist(strsplit(unlist(strsplit(d$v[x], " ")), "\\,")))) })
  args <- d$z # creating vector of each epoch and each corresponding h/l value
  names(args) <- c(d$e)
  
  
  
  # STEPS
  
  # STEP 1: Estimated marginal predictions for each history
  # Gets estimated marginal predictions
  
  preds <- lapply(seq_len(length(model)), function (y) { # Goes through dsifferent fitted model
    final_model <- model[[y]]
    p <- marginaleffects::avg_predictions(final_model,
                                          variables = args,
                                          wts = "(weights)")
    class(p) <- c("pred_custom", class(p))
    p
  })
  
  
  #STEP 2: CONDUCT HISTORY COMPARISONS
  
  if (!is.null(reference)) {
    if (!is.character(reference)){
      stop("Please provide one or more character strings each comprising a valid reference history made up of combinations of 'h' and 'l'",
            call. = FALSE)
    }
    if (length(reference) >= length(exposure_levels)) {
      stop("Please provide at least 1 fewer reference events.",
            call. = FALSE)
    }
    
    if (sum(exposure_levels %in% reference) != length(reference)) {
      stop(paste0('If you wish to conduct custom comparisons, please select a valid reference history from the following list:
                  ', paste(apply(perm2(nrow(epochs), c("l", "h")), 1,
                                 paste, sep = ",", collapse = "-"), sep = ", ", 
                           collapse = ", ")),
            call. = FALSE)
    }
    if (is.null(comparison)) {
      stop(paste0("If you wish to conduct custom comparisons, please specify at least one valid comparison history from the following list:
                  ", paste0(apply(perm2(nrow(epochs), c("l", "h")), 1,
                                  paste, sep = ",", collapse = "-"), sep = " ", 
                            collapse = " "),
                   " otherwise, do not specify reference and comparison events to conduct all comparisons."),
            call. = FALSE)
    }
    if (any(reference %in% comparison)) {
      stop("If you wish to make custom comparisons, please provide distinct reference and comparison events.",
            call. = FALSE)
    }
  }
  
  if (!is.null(comparison)) {
    if (!is.character(comparison)) {
      stop("Please provide as a character a valid comparison history/histories made up of combinations of 'h' and 'l'",
            call. = FALSE)
    }
    if (length(comparison) >= length(exposure_levels)) {
      stop("Please provide at least 1 fewer comparison events.",
            call. = FALSE)
    }
    if (sum(exposure_levels %in% comparison) != length(comparison)) {
      stop(paste0('If you wish to specify comparison(s), please provide at least one comparison history from the following list ',
                   paste0(apply(perm2(nrow(epochs), c("l", "h")), 1,
                                paste, sep = "", collapse = "-"), sep = ", ", 
                          collapse = " "),
                   " otherwise, do not specify reference or comparison events to conduct all comparisons."))
    }
    comp_histories <- exposure_levels[exposure_levels %in% comparison]
  }
  else {
    comp_histories <- NULL
  }
  
  #conduct all pairwise comparisons if no ref/comparisons were specified by user
  
  if (is.null(reference) && is.null(comp_histories)) {
    
    # Pairwise comparisons; don't need to use custom class
    
    comps <- lapply(preds, function(y) {
      marginaleffects::hypotheses(y, "pairwise")
    })
    
    
  }
  else {
    comps <- create_custom_contrasts(d, reference, comp_histories, 
                                     exposure, preds)
  }
  
  
  
  #IMPUTED DATA
  #pooling predicted values and contrasts for imputed data
  
  if (is.null(names(model))) { #imputed data
    rlang::check_installed("mice")
    
    # STEP 1b: Create custom tidy method; not called directly but used in mice::pool()
    
    ## ALERT ## <-- how to send this function to mice::pool() without assigning 
    tidy.pred_custom <<- function(x, ...) {
      out <- NextMethod("tidy", x)
      out$term <- do.call(sprintf, c(paste0(paste(names(args), 
                                                  collapse = " = %s, "), " = %s"),
                                     as.list(unclass(out[, names(args)]))))
      out
    }
    
    # STEP 1c: Pooling predicted estimates
    # Pool results
    
    preds_pool <- summary(mice::pool(mice::as.mira(preds)), digits = 3)
    
    preds_pool <- add_histories(preds_pool, d)
    
    preds_pool <- add_dose(preds_pool, dose_level)
    
    # If the user specified reference and comparison groups, subset pred_pool for inspection and plotting
    
    if (!is.null(reference) && !is.null(comp_histories)) {
      preds_pool <- preds_pool[preds_pool$history %in% c(reference, 
                                                         comp_histories), , 
                               drop = FALSE ]
    }
    
    
    if (save.out) {
      if (is.null(hi_lo_cut)) {
        hi_lo_cut <- "median+-.1"
      }
      
      # Makes table of pooled average estimates and saves out
      
      outfile <- file.path(home_dir, "histories", 
                           sprintf("%s-%s_pooled_estimated_means_hi_lo=%s.html",
                                   exposure, outcome, 
                                   paste(hi_lo_cut, collapse = "_") ))
      sink(outfile)
      stargazer::stargazer(
        as.data.frame(preds_pool),
        type = "html",
        digits = 4,
        column.labels = colnames(preds_pool),
        summary = FALSE,
        rownames = FALSE,
        header = FALSE,
        out = outfile
      )
      sink()
    }
    
    if (verbose) {
      cat("\n")
      cat("Below are the pooled average predictions by user-specified history:") #
      cat(knitr::kable(preds_pool,
                       format = 'pipe',
                       digits = 4),
          sep = "\n")
      cat("\n")
    }
    
    
    
    #STEP 2b: pool comparison values
    #pool summary
    
    # putting data from all ref events into  same entry
    
    comps <- lapply(comps, function(x) {
      do.call(rbind.data.frame, x)
    })
    
    comps_pool <- summary(mice::pool(comps))

    comps_pool <- add_histories(comps_pool, d)
    
    comps_pool <- add_dose(comps_pool, dose_level)
    
    
    #STEP 3: conduct multiple comparison correction
    
    comps_pool <- perform_multiple_comparison_correction(comps_pool, 
                                                         reference, 
                                                         comp_histories, 
                                                         mc_comp_method, 
                                                         verbose)
    
    #rounding term values
    
    comps_pool$term <- unlist(lapply(1:length(comps_pool$term), function(x) {
      x = as.character(comps_pool$term[x])
      a = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                                                               sapply(strsplit(x, "\\ - "), "[", 1)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      b = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                                                               sapply(strsplit(x, "\\ - "), "[", 2)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                               sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      new = paste(a, b, sep = " - ")
      new
    }))
    
    
    if (save.out) {
      if (is.null(hi_lo_cut)) {
        hi_lo_cut <- "median+-.1"
      }
      
      # Makes table of pooled comparisons and saves out
      
      outfile <- file.path(home_dir, "histories", 
                           sprintf("%s-%s_pooled_comparisons_hi_lo=%s.html",
                                   exposure, outcome, 
                                   paste(hi_lo_cut, collapse = "_") ))
      sink(outfile)
      stargazer::stargazer(
        as.data.frame(comps_pool),
        type = "html",
        digits = 4,
        column.labels = colnames(comps_pool),
        summary = FALSE,
        rownames = FALSE,
        header = FALSE,
        out = outfile
      )
      sink()
    }
    if (verbose) {
      cat("\n")
      cat(paste0("USER ALERT: please inspect the following pooled comparisons :"), "\n")
      cat(knitr::kable(comps_pool,
                       format = 'pipe',
                       digits = 4),
          sep = "\n")
      cat("\n")
      cat("\n")
    }
    
    
    # for plotting
    
    comps <- comps_pool
    preds <- preds_pool
    
    
  }
  else { # NON-IMPUTED DATA
    
    #add history and dose labels
    
    preds <- add_histories(preds, d)
    
    preds <- add_dose(preds, dose_level)
    
    # If the user specified reference and comparison groups, subset preds for inspection and plotting
    
    if (!is.null(reference) && !is.null(comp_histories)) {
      preds <- preds[preds$history %in% c(reference, comp_histories), , 
                     drop = FALSE ]
    }
    
    
    if (save.out) {
      if (is.null(hi_lo_cut)) {
        hi_lo_cut <- "median+-.1"
      }
      
      # Makes table of average estimates a
      
      # lapply(seq_len(length(preds)), function(x) {
      #   y <- preds[[x]]
      
      outfile <- file.path(home_dir, "histories", 
                           sprintf("%s-%s_estimated_means_hi_lo=%s.html",
                                   exposure, outcome, 
                                   paste(hi_lo_cut, collapse = "_") ))
      
      sink(outfile)
      stargazer::stargazer(
        as.data.frame(preds),
        type = "html",
        digits = 4,
        column.labels = colnames(preds),
        summary = FALSE,
        rownames = FALSE,
        header = FALSE,
        out = outfile
      )
      sink()
      
      # })
    }
    
    if (verbose) {
      cat("\n")
      cat("Below are the average predictions by user-specified history:", "\n") # Not sure if we need to print this?
      cat(knitr::kable(preds,
                       format = 'pipe',
                       digits = 4),
          sep = "\n")
      cat("\n")
    }
    
    
    
    #STEP 3: conduct multiple comparison correction
    
    # if (length(reference) > 1) {
    comps <- do.call(rbind.data.frame, comps)
    # }
    
    comps <- add_histories(comps, d)
    # comps <- lapply(comps[[1]], function(x){
    #   add_histories(x, d)
    # })
    
    comps <- add_dose(comps, dose_level)
    # comps <- lapply(comps, function(x){
    #   add_dose(x, dose_level)
    # })

    
    
    comps <- perform_multiple_comparison_correction(comps, 
                                                    reference, 
                                                    comp_histories, 
                                                    mc_comp_method, 
                                                    verbose)
    
    
    #rounding term values
    
    comps$term <- unlist(lapply(1:length(comps$term), function(x) {
      x = comps$term[x]
      a = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                                                               sapply(strsplit(x, "\\ - "), "[", 1)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                               sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      b = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                                                               sapply(strsplit(x, "\\ - "), "[", 2)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", 
                                               sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      new = paste(a, b, sep = " - ")
      new
    }))
    
    if (save.out) {
      if (is.null(hi_lo_cut)) {
        hi_lo_cut <- "median+-.1"
      }
      
      outfile <- file.path(home_dir, "histories", 
                           sprintf("%s-%s_comparisons_hi_lo=%s.html",
                                   exposure, outcome, paste(hi_lo_cut, 
                                                            collapse = "_") ))
      
      
      # Makes table of comparisons and saves out
      
      sink(outfile)
      stargazer::stargazer(
        as.data.frame(comps),
        type = "html",
        digits = 4,
        column.labels = colnames(comps),
        summary = FALSE,
        rownames = FALSE,
        header = FALSE,
        out = outfile
      )
      sink()
      
    }
    
    if (verbose) {
      cat("\n")
      cat(paste0("USER ALERT: please inspect the following comparisons:"), "\n")
      cat(knitr::kable(comps,
                       format = 'pipe',
                       digits = 2),
          sep = "\n")
      cat("\n")
      cat("\n")
    }
    
  }
  
  
  
  #STEP 4: Plot results
  
  if (length(colors) > 1 && length(colors) != nrow(epochs) + 1) {
    stop(sprintf("Please provide either %s different colors, a Brewer color palette, or leave this entry blank. \n",
                 nrow(epochs) + 1),
          call. = FALSE)
  }
  
  
  if (is.na(exp_lab)) {
    exp_lab <- exposure
  }
  if (is.na(out_lab)) {
    out_lab <- outcome
  }
  
  comparisons <- data.frame(preds)
  comparisons$low_ci <- comparisons$estimate - (1.96 * comparisons$std.error)
  comparisons$high_ci <- comparisons$estimate + (1.96 * comparisons$std.error)
  comparisons$history <- as.factor(comparisons$history)
  comparisons$dose <- as.factor(comparisons$dose)
  
  comparisons <- comparisons[order(comparisons$dose), , drop = FALSE]
  
  if (length(colors) > 1) { # If user input a list of colors
    p <- ggplot2::ggplot(data = comparisons, ggplot2::aes(x = .data$estimate, 
                                                          y = .data$history, 
                                                          color = .data$dose)) +
      ggplot2::geom_point(size = 5) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_y_discrete(limits = c(as.character(comparisons$history)), 
                                expand = c(0, 0.2)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$low_ci, 
                                           xmax = .data$high_ci), 
                              height = 0.6) +
      ggplot2::xlab(paste0("Predicted ", out_lab, " Value")) +
      ggplot2::ylab(paste0(exp_lab, " Exposure History")) +
      ggplot2::xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci),
                    max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
      ggplot2::theme(text = ggplot2::element_text(size = 14)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), 
                     axis.line = ggplot2::element_line(colour = "black"))
    
    if (save.out) {
      ggplot2::ggsave(file.path(home_dir, "plots", 
                                sprintf("%s-%s.jpeg",
                                        exposure, outcome)),
                      plot = p)
    }
    if (verbose) {
      print(p)
    }
    
  }
  else { # If user lists a palette (default)
    p <- ggplot2::ggplot(data = comparisons, ggplot2::aes(x = .data$estimate, 
                                                          y = .data$history, 
                                                          color = .data$dose)) +
      ggplot2::geom_point(size = 5) +
      ggplot2::scale_colour_brewer(palette = colors) +
      ggplot2::scale_y_discrete(limits = c(as.character(comparisons$history)), 
                                expand = c(0, 0.2)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$low_ci, 
                                           xmax = .data$high_ci), 
                              height = 0.6) +
      ggplot2::xlab(paste0("Predicted ", out_lab, " Value")) +
      ggplot2::ylab(paste0(exp_lab, " Exposure History")) +
      ggplot2::xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci),
                    max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
      ggplot2::theme(text = ggplot2::element_text(size = 14)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), 
                     axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Dosage"))
    
    if (save.out) {
      ggplot2::ggsave(file.path(home_dir, "plots", 
                                sprintf("%s-%s.jpeg",
                                        exposure, outcome
                                )),
                      plot = p)
    }
    if (verbose) {
      print(p)
    }
  }
  
  if (save.out && verbose) {
    cat("\n")
    cat("See the '/plots/' folder for graphical representations of results.")
  }
  
  comps
}




