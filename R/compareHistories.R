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
#' @param reference (optional) string of "-"-separated "l" and "h" values
#'   indicative of a reference exposure history to which to compare comparison,
#'   required if comparison is supplied
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


compareHistories <- function(home_dir, exposure, exposure_time_pts, outcome, model, epochs = NULL, hi_lo_cut = NULL,
                             reference = NA, comparison = NULL, mc_comp_method = "BH", dose_level = "h", exp_lab = NA, out_lab = NA,
                             colors = "Dark2", verbose = TRUE, save.out = TRUE ) {

  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if(!is.character(home_dir)){
      stop("Please provide a valid home directory path as a string if you wish to save output locally.", call. = FALSE)
    }
    else if(!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.", call. = FALSE)
    }
  }

  if (missing(exposure)){
    stop("Please supply a single exposure.", call. = FALSE)
  }
  else if(!is.character(exposure) || length(exposure) != 1){
    stop("Please supply a single exposure as a character.", call. = FALSE)
  }

  if (missing(outcome)){
    stop("Please supply a single outcome.", call. = FALSE)
  }
  else if(!is.character(outcome) || length(outcome) != 1){
    stop("Please supply a single outcome as a character.", call. = FALSE)
  }

  if (missing(exposure_time_pts)){
    stop("Please supply the exposure time points at which you wish to create weights.", call. = FALSE)
  }
  else if(!is.numeric(exposure_time_pts)){
    stop("Please supply a list of exposure time points as integers.", call. = FALSE)
  }

  if (missing(model)){
    stop("Please supply a list of model output", call. = FALSE)
  }
  else if(!is.list(model) || is.data.frame(model)){
    stop("Please provide a list of model output from the fitModel function.", call. = FALSE)
  }
  else if (sum(sapply(model, function(x) {
    inherits(x, "svyglm")})) != length(model)){
    stop("Please supply a model as a list from the createWeights function.", call. = FALSE)
  }

  if(!is.logical(verbose)){
    stop("Please set verbose to either TRUE or FALSE.", call. = FALSE)
  }
  else if(length(verbose) != 1){
    stop("Please provide a single TRUE or FALSE value to verbose.", call. = FALSE)
  }

  if(!is.logical(save.out)){
    stop("Please set save.out to either TRUE or FALSE.", call. = FALSE)
  }
  else if(length(save.out) != 1){
    stop("Please provide a single TRUE or FALSE value to save.out.", call. = FALSE)
  }



  if(save.out){
    histories_dir <- file.path(home_dir, "histories")
    if (!dir.exists(histories_dir)) {
      dir.create(histories_dir)
    }
    plot_dir <- file.path(home_dir, "plots")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
  }

  exposure_type <- if(inherits(model[[1]]$data[, paste0(exposure, '.', exposure_time_pts[1])], "numeric"))
    "continuous" else "binary"


  ints <- gsub(" ", "", as.character(unlist(strsplit(as.character(unlist(model[[1]]$terms)), "\\+"))))
  ints <- ifelse(sum(grepl(":", ints)) > 0, 1, 0)


  #history inspection --done early bc function deals with epochs, etc.
  #print history sample distribution
  if (verbose){
    eval_hist(data = model[[1]]$data, exposure, epochs,
              exposure_time_pts, hi_lo_cut, reference, comparison, verbose)
  }


  #getting data to use for determining hi/lo values: should have any epochs created/ used in the model
  if (is.null(names(model))){ #imputed data
    # stacking imputed data to compute hi/lo vals across all imps
    data <- lapply(model, function(x) { x$data })
    data <- do.call("rbind", data)
  }
  else if(!is.null(names(model))){ #single df (not imputed)
    data <- model[[1]]$data
  }

  if(is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }
  else{

    if( !is.data.frame(epochs) || ncol(epochs) != 2 || sum(colnames(epochs) == c("epochs", "values")) != ncol(epochs)){
      stop("If you supply epochs, please provide a dataframe with two columns of epochs and values.",
           call. = FALSE)
    }
    if(sum(is.na(epochs$values)) > 0){
      stop("Please provide one or a list of several values for each epoch.", call. = FALSE)
    }
  }

  terms  <- model[[1]]$formula
  terms <- gsub(" ", "", as.character(unlist(strsplit(as.character(unlist(terms[3])), "\\+"))))
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, paste, sep = "", collapse = ".")

  if (sum(exp_epochs %in% terms) != length(exp_epochs)){
    stop("The use of exposure epochs must be consistent across fitModel and compareHistories, given that the fitted model output is used as the basis of histories.
         If you specified exposure epochs in fitModel, please specify the same ones at this step.
         If you did not, please do not specify them at this step.",
         call. = FALSE)
  }


  # creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  exposure_levels <- apply(gtools::permutations(2, nrow(epochs), c("l", "h"),
                                                repeats.allowed = TRUE), 1, paste, sep = "", collapse = "-")
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, paste, sep = "", collapse = ".")

  eps <- epochs$epochs
  # gathering epoch information for each exposure for deriving betas
  epoch_info <- as.data.frame(rep(exposure, length(eps)))
  epoch_info$time <- eps
  is.na(epoch_info$low) <- TRUE
  is.na(epoch_info$high) <- TRUE

  if (exposure_type == "continuous") {

    #finds high and low values for each epoch
    for (t in seq_len(length(eps))) {
      var_name <- paste(exposure, eps[t], sep = ".")

      if(is.null(hi_lo_cut)){ #default is median split
        epoch_info$low[t] <- as.numeric(median(data[, var_name] - 0.001, na.rm = T)) #padded with -.001 bc otherwise breaks hypotheses()
        epoch_info$high[t] <- as.numeric(median(data[, var_name] + 0.001 , na.rm = T))

      } else if (!inherits(hi_lo_cut, "numeric")){
        stop("Please provide one or two numbers to hi_lo_cut between 0-1", call. = FALSE)

      } else if (is.numeric(hi_lo_cut)){
        if (length(hi_lo_cut) > 2){
          stop("Please provide either one or two numeric values between 0-1 to hi_lo_cut.", call. = FALSE)

        }else if (length(hi_lo_cut) == 2){ #if two cutoff values supplied
          hi_cutoff <- hi_lo_cut[1]
          lo_cutoff <- hi_lo_cut[2]

          if (hi_cutoff > 1 || hi_cutoff < 0) {
            stop('Please supply a high cutoff value to hi_lo_cut between 0 and 1', call. = FALSE)
          }
          if (lo_cutoff > 1 || lo_cutoff < 0) {
            stop('Please supply a low cutoff value to hi_lo_cut between 0 and 1', call. = FALSE)
          }
        } else{ #if only one cutoff value supplied
          if (hi_lo_cut > 1 || hi_lo_cut < 0) {
            stop('Please select a hi_lo cutoff value between 0 and 1', call. = FALSE)
          }

          hi_cutoff <- hi_lo_cut
          lo_cutoff <- hi_lo_cut
        }

        #find user-specfied quantile values
        epoch_info$low[t] <- as.numeric(quantile(data[, var_name], probs = lo_cutoff, na.rm = T))
        epoch_info$high[t] <- as.numeric(quantile(data[, var_name], probs = hi_cutoff, na.rm = T))
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
  d$z <- lapply(seq_len(nrow(d)), function(x) {
    c(as.numeric(unlist(strsplit(unlist(strsplit(d$v[x], " ")), "\\,")))) })
  args <- d$z # creating vector of each epoch and each corresponding h/l value
  names(args) <- (c(d$e))



  # STEPS

  # STEP 1: Estimated marginal predictions for each history
  # Gets estimated marginal predictions
  preds <- lapply(seq_len(length(model)), function(y) { # Goes through dsifferent fitted model
    final_model <- model[[y]]
    p <- marginaleffects::avg_predictions(final_model,
                                          variables = args,
                                          wts = "(weights)")
    class(p) <- c("pred_custom", class(p))
    p
  })


  #STEP 2: CONDUCT HISTORY COMPARISONS

  if (!is.na(reference)) {
    if (!inherits(reference, "character")){
      stop("Please provide as a character a valid reference history made up of combinations of 'h' and 'l'",
           call. = FALSE)
    }
    if (sum(exposure_levels %in% reference) == 0) {
      stop(paste0('If you wish to conduct custom comparisons, please select a valid reference history from the following list:
                  ',
                  paste(apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                              paste, sep = ",", collapse = "-"), sep = ", ", collapse = ", ")),
           call. = FALSE)
    }
    if (is.null(comparison)) {
      stop(paste0("If you wish to conduct custom comparisons, please specify at least one valid comparison history from the following list:
                  ",
                  paste0(apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                               paste, sep = ",", collapse = "-"), sep = " ", collapse = " "),
                  " otherwise, do not specify reference and comparison events to conduct all comparisons."),
           call. = FALSE)
    }
    else if (reference %in% comparison) {
      stop("If you wish to make a custom comparison, please provide unique reference and comparison events.",
           call. = FALSE)
    }
  }

  if (!is.null(comparison)) {
    if (!inherits(comparison, "character")){
      stop("Please provide as a character a valid comparison history/histories made up of combinations of 'h' and 'l'",
           call. = FALSE)
    }
    if (sum(exposure_levels %in% comparison) == 0) {
      stop(paste0('If you wish to specify comparison(s), please provide at least one comparison history from the following list ',
                  paste0(apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                               paste, sep = "", collapse = "-"), sep = ", ", collapse = " "),
                  " otherwise, do not specify reference or comparison events to conduct all comparisons."))
    }
    comp_histories <- exposure_levels[exposure_levels %in% comparison]
  }
  else {
    comp_histories <- NULL
  }

  #conduct all pairwise comparisons if no ref/comparisons were specified by user
  if (is.na(reference) && is.null(comp_histories)){
    # Pairwise comparisons; don't need to use custom class
    comps <- lapply(preds, function(y){
      y |> marginaleffects::hypotheses("pairwise")
    })

  } else{
    comps <- create_custom_contrasts(d, reference, comp_histories, exposure, preds)
  }


  #IMPUTED DATA
  #pooling predicted values and contrasts for imputed data
  if (is.null(names(model))){ #imputed data

    # STEP 1b: Create custom tidy method; not called directly but used in mice::pool()
    tidy.pred_custom <<- function(x, ...) {
      out <- NextMethod("tidy", x)
      out$term <- do.call(sprintf, c(paste0(paste(names(args), collapse = " = %s, "), " = %s"),
                                     as.list(unclass(out[, names(args)]))))
      out
    }

    # STEP 1c: Pooling predicted estimates
    # Pool results
    preds_pool <- mice::pool(mice::as.mira(preds)) |> summary(digits = 3)

    preds_pool <- add_histories(preds_pool, d)

    preds_pool <- add_dose(preds_pool, dose_level)

    # If the user specified reference and comparison groups, subset pred_pool for inspection and plotting
    if (!is.na(reference) && !is.null(comp_histories)) {
      # preds_pool <- preds_pool %>%
      #   dplyr::filter(history %in% c(reference, comp_histories))
      preds_pool <- preds_pool[preds_pool$history %in% c(reference, comp_histories), , drop = FALSE ]
    }


    if (save.out){
      if(is.null(hi_lo_cut)){
        hi_lo_cut <- "median+-.1"
      }

      # Makes table of pooled average estimates and saves out
      outfile <- sprintf("%s/histories/%s-%s_pooled_estimated_means_hi_lo=%s.html",
                         home_dir, exposure, outcome, paste(hi_lo_cut, collapse = "_") )
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

    if (verbose){
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
    comps_pool <- mice::pool(comps) |> summary()

    comps_pool <- add_histories(comps_pool, d)

    comps_pool <- add_dose(comps_pool, dose_level)


    #STEP 3: conduct multiple comparison correction
    comps_pool <- perform_multiple_comparison_correction(comps_pool, reference, comp_histories, mc_comp_method, verbose)

    #rounding term values
    comps_pool$term <- unlist(lapply(1:length(comps_pool$term), function(x){
      x = as.character(comps_pool$term[x])
      a = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      b = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 2)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      new = paste(a, b, sep = " - ")
      new
    }))


    if (save.out){
      if(is.null(hi_lo_cut)){
        hi_lo_cut <- "median+-.1"
      }

      # Makes table of pooled comparisons and saves out
      outfile <- sprintf("%s/histories/%s-%s_pooled_comparisons_hi_lo=%s.html",
                         home_dir, exposure, outcome, paste(hi_lo_cut, collapse = "_") )
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
    if (verbose){
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


  } else{ # NON-IMPUTED DATA

    #add history and dose labels
    preds <- add_histories(preds, d)

    preds <- add_dose(preds, dose_level)

    # If the user specified reference and comparison groups, subset preds for inspection and plotting
    if (!is.na(reference) && !is.null(comp_histories)) {
      # preds <- preds %>%
      #   dplyr::filter(history %in% c(reference, comp_histories))
      preds <- preds[preds$history %in% c(reference, comp_histories), , drop = FALSE ]
    }


    if (save.out){
      if(is.null(hi_lo_cut)){
        hi_lo_cut <- "median+-.1"
      }

      # Makes table of average estimates a
      lapply(seq_len(length(preds)), function(x) {
        y <- preds[[x]]

        outfile <- sprintf("%s/histories/%s-%s_estimated_means_hi_lo=%s.html",
                           home_dir, exposure, outcome, paste(hi_lo_cut, collapse = "_") )

        sink(outfile)
        stargazer::stargazer(
          as.data.frame(y),
          type = "html",
          digits = 4,
          column.labels = colnames(y),
          summary = FALSE,
          rownames = FALSE,
          header = FALSE,
          out = outfile
        )
        sink()

      })
    }

    if (verbose){
      cat("\n")
      cat("Below are the average predictions by user-specified history:", "\n") # Not sure if we need to print this?
      cat(knitr::kable(preds,
                       format = 'pipe',
                       digits = 4),
          sep = "\n")
      cat("\n")
    }



    #STEP 3: conduct multiple comparison correction
    comps <- comps[[1]]

    comps <- add_histories(comps, d)

    comps <- add_dose(comps, dose_level)

    comps <- perform_multiple_comparison_correction(comps, reference, comp_histories, mc_comp_method, verbose)

    #rounding term values
    comps$term <- unlist(lapply(1:length(comps$term), function(x){
      x = comps$term[x]
      a = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      b = paste0("(", paste(round(as.numeric(as.character(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 2)), " "))[
        !is.na(as.numeric(unlist(strsplit(gsub("[^0-9.-]", " ", sapply(strsplit(x, "\\ - "), "[", 1)), " "))))])), 4),
        collapse = ", ", sep = ", "), ")")
      new = paste(a, b, sep = " - ")
      new
    }))

    if (save.out){
      if(is.null(hi_lo_cut)){
        hi_lo_cut <- "median+-.1"
      }

      outfile <- sprintf("%s/histories/%s-%s_comparisons_hi_lo=%s.html",
                         home_dir, exposure, outcome, paste(hi_lo_cut, collapse = "_") )

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
  if(!inherits(colors, "character")){
    stop("Please provide a character string of a Brewer palette name or list of colors for plotting.", call. = FALSE)
  }
  else if (length(colors) > 1 && length(colors) != nrow(epochs) + 1) {
    stop(
      # paste0('Please provide either: ', nrow(epochs) + 1,
      #           ' different colors, a color palette, or leave this entry blank.'),
      sprint("Please provide either %s different colors, a Brewer color palette, or leave this entry blank. \n",
             nrow(epochs) + 1),
      call. = FALSE)
  }

  if (is.na(exp_lab)){
    exp_lab <- exposure
  }
  if(is.na(out_lab)){
    out_lab <- outcome
  }

  comparisons <- data.frame(preds)
  comparisons$low_ci <- comparisons$estimate - (1.96 * comparisons$std.error)
  comparisons$high_ci <- comparisons$estimate + (1.96 * comparisons$std.error)
  comparisons$history <- as.factor(comparisons$history)
  comparisons$dose <- as.factor(comparisons$dose)

  # comparisons <- comparisons %>%
  #   dplyr::arrange(dose) # Order by dose
  comparisons <- comparisons[order(comparisons$dose), , drop = FALSE]

  if (length(colors) > 1) { # If user input a list of colors
    p <- ggplot2::ggplot(data = comparisons, ggplot2::aes(x = estimate, y = history, color = dose)) +
      ggplot2::geom_point(size = 5) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_y_discrete(limits = c(as.character(comparisons$history)), expand = c(0, 0.2)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = low_ci, xmax = high_ci), height = 0.6) +
      ggplot2::xlab(paste0("Predicted ", out_lab, " Value")) +
      ggplot2::ylab(paste0(exp_lab, " Exposure History")) +
      ggplot2::xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci),
                    max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
      ggplot2::theme(text = ggplot2::element_text(size = 14)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))

    if (save.out){
      ggplot2::ggsave(
        # paste0(home_dir, "/plots/", exposure, "-", outcome, ".jpeg"),
        sprintf("%s/plots/%s-%s.jpeg",
                home_dir, exposure, outcome),
        plot = p)
    }
    if(verbose){
      print(p)
    }

  } else { # If user lists a palette (default)
    p <- ggplot2::ggplot(data = comparisons, ggplot2::aes(x = estimate, y = history, color = dose)) +
      ggplot2::geom_point(size = 5) +
      ggplot2::scale_colour_brewer(palette = colors) +
      ggplot2::scale_y_discrete(limits = c(as.character(comparisons$history)), expand = c(0, 0.2)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = low_ci, xmax = high_ci), height = 0.6) +
      ggplot2::xlab(paste0("Predicted ", out_lab, " Value")) +
      ggplot2::ylab(paste0(exp_lab, " Exposure History")) +
      ggplot2::xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci),
                    max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
      ggplot2::theme(text = ggplot2::element_text(size = 14)) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::guides(fill = ggplot2::guide_legend(title="Dosage"))

    if (save.out){
      ggplot2::ggsave(
        # paste0(home_dir, "/plots/", exposure, "-", outcome, ".jpeg"),
        sprintf("%s/plots/%s-%s.jpeg",
                home_dir, exposure, outcome),
        plot = p)
    }
    if(verbose){
      print(p)
    }
  }

  if (save.out && verbose){
    cat("\n")
    cat("See the '/plots/' folder for graphical representations of results.")
  }

  comps
}




