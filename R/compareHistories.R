#' Compare exposure histories
#' This code uses the best-fitting model for each exposure-outcome pair to compare the effects of user-specified reference and comparison histories of exposure on outcome using linear hypothesis testing
#' @param object msm object that contains all relevant user inputs
#' @param data_for_model_with_weights_cutoff imputed datasets with truncated weights
#' @param all_models fitted models for each imputed dataset
#' @param reference optional reference event for custom comparison
#' @param compare optional comparison event(s) for custom comparison
#' @importFrom gtools permutations
#' @importFrom marginaleffects avg_predictions
#' @importFrom marginaleffects hypotheses
#' @importFrom stringr str_count
#' @importFrom dplyr filter
#' @importFrom stargazer stargazer
#' @importFrom mice pool
#' @importFrom stats p.adjust
#' @importFrom knitr kable
#' @return preds_pool
#' @examples compareHistories(object, data_for_model_with_weights_cutoff, all_models, reference=NA, compare=NA)

compareHistories <- function(home_dir, exposure, outcome, tv_confounders, model, epochs = NULL, hi_lo_cut = NA, reference = NA, comparison = NULL, mc_comp_method = "BH", dose_level = "h", exp_lab = NA, out_lab = NA, colors = "Dark2" ) {

  #error checking
  if (!dir.exists(home_dir)) {
    stop("Please provide a valid home directory path.")
  }

  histories_dir <- file.path(home_dir, "histories")
  if (!dir.exists(histories_dir)) {
    dir.create(histories_dir)
  }

  plot_dir <- file.path(home_dir, "plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }

  exposure_type <- ifelse(class(model[[1]]$data[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

  exposure_time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }

  if (length(colors) > 1 & length(colors) != nrow(epochs) + 1) {
    stop(paste0('Please provide either: ', nrow(epochs) + 1,
                ' different colors, a color palette, or leave this entry blank.'))
  }


  #history inspection
  #print history sample distribution --check this
  eval_hist(data = model[[1]]$data, exposure, tv_confounders, epochs,
            time_pts= as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2)),
            hi_lo_cut, reference, comparison)

  # gathering epoch info
  eps <- epochs$epochs
  # creates permutations of high ("h") and low ("l") levels of exposure for each exposure epoch
  exposure_levels <- apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1, paste, sep = "", collapse = "-")
  exp_epochs <- apply(expand.grid(exposure, as.character(epochs[, 1])), 1, paste, sep = "", collapse = "_")


  if (!is.na(reference)) {
    if (sum(exposure_levels %in% reference) == 0) {
      stop(paste0('If you wish to conduct custom comparisons, please select a valid reference history from the following list ',
                  paste(apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                              paste, sep = ",", collapse = "-"), sep = ", ", collapse = ", ")))
    }
    if (is.null(comparison)) {
      stop(paste0("If you wish to conduct custom comparisons, please specify at least one valid comparison history from the following list ",
                  paste0(apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                               paste, sep = ",", collapse = "-"), sep = " ", collapse = " "),
                  " otherwise, do not specify reference and comparison events to conduct all comparisons."))
    } else if (reference %in% comparison) {
      stop("If you wish to make a custom comparison, please provide unique reference and comparison events.")
    }
  }

  if (!is.null(comparison)) {
    if (sum(exposure_levels %in% comparison) == 0) {
      stop(paste0('If you wish to specify comparison(s), please provide at least one comparison history from the following list ',
                  paste0(apply(gtools::permutations(2, nrow(epochs), c("l", "h"), repeats.allowed = TRUE), 1,
                               paste, sep = "", collapse = "-"), sep = ", ", collapse = " "),
                  " otherwise, do not specify reference and comparison events to conduct all comparisons."))
    }
    comp_histories <- exposure_levels[exposure_levels %in% comparison]
  } else {
    comp_histories <- NA
  }

  ints <- gsub(" ", "", as.character(unlist(strsplit(as.character(unlist(model[[1]]$terms)), "\\+"))))
  ints <- ifelse(sum(grepl(":", ints)) > 0, 1, 0)
  # gathering epoch information for each exposure for deriving betas
  epoch_info <- as.data.frame(rep(exposure, length(eps)))
  epoch_info$time <- eps
  epoch_info$low <- NA
  epoch_info$high <- NA


  # cycling through eps to find hi and lo values of the exposure for each epoch based on user-specified values
  if (exposure_type == "continuous") {

    if (is.null(names(model))){ #imputed data
      # stacking imputed data to compute hi/lo vals across all
      data <- lapply(model, function(x) { x$data })[[1]]
      data <- do.call("rbind", data)
    }

    if(names(model) == "0"){ #single df (not imputed)
      data <- model[[1]]$data
    }

    for (t in 1:length(eps)) {
      var_name <- paste(exposure, eps[t], sep = "_")

      if(is.na(hi_lo_cut)){
        epoch_info$low[t] <- as.numeric(median(data[, var_name], na.rm = T))
        epoch_info$high[t] <- as.numeric(median(data[, var_name], na.rm = T))

      } else{

        hi_cutoff <- hi_lo_cut[1]
        lo_cutoff <- hi_lo_cut[2]

        if (hi_cutoff > 1 || hi_cutoff < 0) {
          stop('Please select a high cutoff value between 0 and 1')
        }
        if (lo_cutoff > 1 || lo_cutoff < 0) {
          stop('Please select low cutoff value between 0 and 1')
        }

        epoch_info$low[t] <- as.numeric(quantile(data[, var_name], probs = lo_cutoff, na.rm = T))
        epoch_info$high[t] <- as.numeric(quantile(data[, var_name], probs = hi_cutoff, na.rm = T))
      }
    }
  }

  if (exposure_type == "binary") {
    for (t in 1:length(eps)) {
      var_name <- paste(exposure, eps[t], sep = "_")
      epoch_info$low[t] <- 0
      epoch_info$high[t] <- 1
    }
  }


  # gather high and low values for each exposure epoch based on quantile values
  d <- data.frame(e = paste(epoch_info[[1]], epoch_info[[2]], sep = "_"),
                  l = epoch_info$low,
                  h = epoch_info$high)
  d$v <- paste(d$l, d$h, sep=",")
  d$z <- lapply(1:nrow(d), function(x) { c(as.numeric(unlist(strsplit(unlist(strsplit(d$v[x], " ")), "\\,")))) })
  args <- d$z # creating vector of each epoch and each corresponding h/l value
  names(args) <- (c(d$e))


  # STEP 1: Estimated marginal predictions for each history
  # Gets estimated marginal predictions
  preds <- lapply(1:length(model), function(y) { # Goes through different fitted model
    final_model <- model[[y]]
    p <- marginaleffects::avg_predictions(final_model,
                                          variables = args,
                                          wts = "(weights)")
    class(p) <- c("pred_custom", class(p))
    p
  })


  #STEP 2: CONDUCT HISTORY COMPARISONS
  #conduct all pairwise comparisons if no ref/comparisons were specified by user
  if (sum(is.na(reference) & is.na(comp_histories))==1){
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
    preds_pool <- mice::pool(mice::as.mira(y)) |> summary(digits = 3)

    preds_pool <- add_histories(preds_pool, d)

    preds_pool <- add_dose(preds_pool, dose_level)

    # If the user specified reference and comparison groups, subset pred_pool for inspection and plotting
    if (!is.na(reference) & sum(is.na(comp_histories)) == 0) {
      preds_pool <- preds_pool %>% dplyr::filter(history %in% c(reference, comp_histories))

    }

    cat("\n")
    cat("Below are the pooled average predictions by user-specified history:") # Not sure if we need to print this?
    print(preds_pool)

    # Makes table of average estimates and saves out
    sink(paste0(home_dir, "/histories/", exposure, "-", outcome, "_estimated_means_",
                "_hi=", hi_cutoff, "_lo=", lo_cutoff, ".html"))
    stargazer::stargazer(
      as.data.frame(y),
      type = "html",
      digits = 4,
      column.labels = colnames(y),
      summary = FALSE,
      rownames = FALSE,
      header = FALSE,
      out = paste0(home_dir, "/histories/", exposure, "-", outcome, "_estimated_means_",
                   "_hi=", hi_cutoff, "_lo=", lo_cutoff, ".html"
      )
    )
    sink()


    #STEP 2b: pool comparison values
    #pool summary
    comps_pool <- mice::pool(comps) |> summary()

    comps_pool <- add_histories(comps_pool, d)

    comps_pool <- add_dose(comps_pool, dose_level)


    #STEP 3: conduct multiple comparison correction
    comps_pool <- perform_multiple_comparison_correction(comps_pool, method)

    cat(paste0("USER ALERT: please inspect the following comparisons :"), "\n")
    cat(knitr::kable(comps_pool, format = 'pipe', digits = 2), sep = "\n")
    cat("\n")
    cat("\n")


    # for plotting
    comps <- comps_pool
    preds <- preds_pool


  } else{ # NON-IMPUTED DATA

    #add history and dose labels
    preds <- add_histories(preds, d)

    preds <- add_dose(preds, dose_level)

    # If the user specified reference and comparison groups, subset preds for inspection and plotting
    if (!is.na(reference) & sum(is.na(comp_histories)) == 0) {
      preds <- lapply(preds, function(y) {
        as.data.frame(y) %>% dplyr::filter(history %in% c(reference, comp_histories))
      })
    }

    cat("\n")
    cat("Below are the average predictions by user-specified history:") # Not sure if we need to print this?
    print(preds)

    # Makes table of average estimates and saves out per cutoff value
    lapply(1:length(preds), function(x) {
      y <- preds[[x]]
      sink(paste0(home_dir, "/histories/", exposure, "-", outcome, "_estimated_means_",
                  names(preds)[x], "_hi=", hi_cutoff, "_lo=", lo_cutoff, ".html"))
      stargazer::stargazer(
        as.data.frame(y),
        type = "html",
        digits = 4,
        column.labels = colnames(y),
        summary = FALSE,
        rownames = FALSE,
        header = FALSE,
        out = paste0(home_dir, "/histories/", exposure, "-", outcome, "_estimated_means_",
                     names(preds)[x], "_hi=", hi_cutoff, "_lo=", lo_cutoff, ".html"
        )
      )
      sink()
    })


    #STEP 3: conduct multiple comparison correction

    comps <- add_histories(comps, d)

    comps <- add_dose(comps, dose_level)

    comps <- perform_multiple_comparison_correction(comps, method)


    cat(paste0("USER ALERT: please inspect the following comparisons :"), "\n")
    cat(knitr::kable(comps, format = 'pipe', digits = 2), sep = "\n")
    cat("\n")
    cat("\n")


  }





  #STEP 4: Plot results

  if (is.na(exp_lab)){
    exp_lab <- exposure
  }
  if(is.na(out_lab)){
    out_lab <- outcome
  }

  for (x in seq_along(preds)) {
    comparisons <- data.frame(preds)
    comparisons$term <- gsub(paste0(exposure, "_"), "", comparisons$term)
    comparisons$low_ci <- comparisons$estimate - (1.96 * comparisons$std.error)
    comparisons$high_ci <- comparisons$estimate + (1.96 * comparisons$std.error)
    comparisons$history <- as.factor(comparisons$history)
    comparisons$dose <- as.factor(comparisons$dose)

    comparisons <- comparisons %>%
      arrange(dose) # Order by dose

    if (length(colors) > 1) { # If user input a list of colors
      p <- ggplot(data = comparisons, aes(x = estimate, y = history, color = dose)) +
        geom_point(size = 5) +
        scale_color_manual(values = colors) +
        scale_y_discrete(limits = c(as.character(comparisons$history)), expand = c(0, 0.2)) +
        geom_errorbarh(aes(xmin = low_ci, xmax = high_ci), height = 0.6) +
        xlab(paste0("Predicted ", outcome_labels, " Value")) +
        ylab(paste0(exposure_labels, " Exposure History")) +
        xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci), max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
        theme(text = element_text(size = 14)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

      ggsave(paste0(home_dir, "/plots/", exposure, "-", outcome, ".jpeg"), plot = p)

    } else { # If user lists a palette
      p <- ggplot(data = comparisons, aes(x = estimate, y = history, color = dose)) +
        geom_point(size = 5) +
        scale_color_palette(name = "Dosage") +
        scale_y_discrete(limits = c(as.character(comparisons$history)), expand = c(0, 0.2)) +
        geom_errorbarh(aes(xmin = low_ci, xmax = high_ci), height = 0.6) +
        xlab(paste0("Predicted ", outcome_labels, " Value")) +
        ylab(paste0(exposure_labels, " Exposure History")) +
        xlim(min(comparisons$low_ci) - 1 * sd(comparisons$low_ci), max(comparisons$high_ci) + 1 * sd(comparisons$high_ci)) +
        theme(text = element_text(size = 14)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))

      ggsave(paste0(home_dir, "/plots/", exposure, "-", outcome, ".jpeg"), plot = p)
    }

    message("\n")
    message("See the '/plots/' folder for graphical representations of results.")
  }







  # FUNCTIONS

  #custom contrasts
  get_reference_values <- function(d, reference) {
    ref_vals <- sapply(1:length(unlist(strsplit(reference, "-"))), function(x) {
      d[x,unlist(strsplit(reference, "-"))[x]]
    })
    ref_vals
  }


  get_comparison_values <- function(d, comp_histories) {
    comp_vals <- sapply(comp_histories, function(comp) {
      sapply(1:length(unlist(strsplit(comp, "-"))), function(x) {
        d[x, unlist(strsplit(comp, "-"))[x]]
      })
    })
    return(t(comp_vals))
  }


  create_custom_contrasts <- function(d, reference, comp_histories, exposure, preds) {
    if (is.na(reference) | any(is.na(comp_histories))) {
      return(NULL)  # Invalid input, return early
    }

    ref_vals <- get_reference_values(d, reference)
    comp_vals <- get_comparison_values(d, comp_histories)
    cus_comps <- create_custom_comparisons(preds, ref_vals, comp_vals, exposure)

    comps <- lapply(preds, function(y) {
      y |> marginaleffects::hypotheses(cus_comps)
    })

    return(comps)
  }


  create_custom_comparisons <- function(preds, ref_vals, comp_vals, exposure) {
    cus_comps <- matrix(ncol = nrow(comp_vals), nrow = nrow(as.data.frame(preds[[1]][[1]])))

    ref_pos <- which(apply(as.data.frame(preds[[1]])[grepl(exposure, colnames(as.data.frame(preds[[1]])))], 1,
                           paste, collapse = ",") == paste0(ref_vals, collapse = ","))
    cus_comps[ref_pos, ] <- -1

    for (x in 1:nrow(comp_vals)) {
      c <- paste(comp_vals[x, ], collapse = ",")
      cus_comps[which(apply(as.data.frame(preds[[1]])[grepl(exposure, colnames(as.data.frame(preds[[1]])))], 1,
                            paste, collapse = ",") == paste0(c, collapse = ",")), x] <- 1
    }

    if (nrow(comp_vals) > 1) {
      colnames(cus_comps) <- paste0("(", paste0(paste0(ref_vals, collapse = ","), ") - (", apply(comp_vals, 1, paste, collapse = ",")), ")")
    } else {
      colnames(cus_comps) <- paste0("(", paste0(paste0(ref_vals, collapse = ","), ") - (", paste0(paste0(comp_vals, collapse = ",")), ")"))
    }

    cus_comps[is.na(cus_comps)] <- 0
    return(cus_comps)
  }



  add_histories <- function(p, d){
    if( "list" %in% class(p) & nrow(p) == 1){
      history <- matrix(data = NA, nrow = nrow(p[[1]]), ncol = 1) # Get histories from the first element
      p <- p[[1]]
    }

    if (sum(d$e %in% colnames(p)) == nrow(d)){
      for (i in 1:nrow(p)) {
        vals <- as.data.frame(p)[i, 1:nrow(d)]
        history[i] <- as.data.frame(paste(ifelse(round(vals, 3) == round(d$l, 3), "l", "h"), collapse = "-"))
      }
      history <- unlist(history)
    }

    if("term" %in% colnames(p)){
      history <- matrix(data = NA, nrow = nrow(p), ncol = 1) # Get histories from the first element

      if(grepl("\\=", p$term[1])){
        for (i in 1:nrow(p)){
          vals <- as.numeric(sapply(strsplit(unlist(strsplit(as.character(p$term[i]), "\\,")), "="), "[",2))
          history[i] <- as.data.frame(paste(ifelse(round(vals,3)==round(d$l,3), "l", "h"), collapse="-"))
        }
        history <- unlist(history)

      }else{
        for (i in 1:nrow(p)) {
          temp <- as.character(p$term[i])
          pair <- lapply(1:2, function(y) {
            a <- sapply(strsplit(temp, " - "), "[", y)
            his <- lapply(1:nrow(d), function(z) {
              ifelse(round(as.numeric(gsub("[^0-9.-]", "", sapply(strsplit(a, "\\,"), "[", z))), 3) == round(d[z, "l"], 3), "l", "h")
            })
          })
          history[i, 1] <- paste(sapply(pair, paste, collapse = "-"), collapse = " vs ")
        }
      }

    }

    p <- cbind (p, history = history)

    p
  }


  add_dose <- function(p, dose_level){
    if( length(p$history[1]) == 1 ){
      if(grepl("vs", p$history[1])){
        dose_a <- stringr::str_count(sapply(strsplit(p$history, "vs"), "[", 1), dose_level)
        dose_b <- stringr::str_count(sapply(strsplit(p$history, "vs"), "[", 2), dose_level)
        dose_count <- data.frame(dose = gsub(" ", " vs ", paste(dose_a, dose_b)))
      } else{

        dose_count <- stringr::str_count(p$history, dose_level)

      }
    }

    if (length(p$history[1]) > 1){
      dose_count <- stringr::str_count(p$history, dose_level)
    }

    # dose_count <- rep(list(dose_count), length(p))
    p <- cbind (p, dose_count = dose_count)
    p
  }


  perform_multiple_comparison_correction <- function(comps, method) {
    if (any(is.na(reference) & is.na(comp_histories)) | length(comp_histories) > 1) {
      cat("\n")
      cat(paste0("Conducting multiple comparison correction using the ", method, " method."), "\n")
      cat("\n")
      corr_p <- lapply(comps, function(x) {
        stats::p.adjust(x$p.value, method = method)
      })
      comps <- Map(cbind, comps, p.value_corr = corr_p)
    } else {
      cat(paste0("The user specified comparison only between ", reference, " and ", comp_histories,
                 " so no correction for multiple comparisons will be implemented."), "\n")
    }

    return(comps)
  }



  comps

}




