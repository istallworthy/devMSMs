#' Calculate balance stats based on Jackson paper
#'
#' Calculate weighted or unweighted standardized balance statistics for a given exposure time point,
#' using all relevant confounders. Draws on Jackson, 2016 approaches to
#' assessing balance for time-varying exposures by weighting statistics based on
#' sample distribution in exposure histories.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 guide_axis
#' @importFrom ggplot2 ggsave
#' @importFrom stargazer stargazer
#' @param home_dir (optional) path to home directory (required if save.out = TRUE)
#' @param data data in wide format as: a data frame, path to folder of imputed
#'   .csv files, or mids object
#' @param formulas list of balancing formulas at each time point output from
#'   createFormulas()
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure wass measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param balance_thresh (optional) one or two numbers between 0 and 1 indicating a single
#'   balancingn threshold or thresholds for more and less important confounders,
#'   respectively
#' @param k (optional) imputation number
#' @param weights (optional) list of IPTW weights output from createWeights
#' @param imp_conf (optional) list of variable names reflecting important confounders (required if two balance thresholds are provided)
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and intermediary output locally
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
#'                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                    formulas = f,
#'                    save.out = FALSE)
#'
#' c <- calcBalStats(data = test,
#'                   formulas = f,
#'                   exposure = "A",
#'                   exposure_time_pts = c(1, 2, 3),
#'                   outcome = "D.3",
#'                   balance_thresh = 0.1)
#' c <- calcBalStats(data = test,
#'                   formulas = f,
#'                   exposure = "A",
#'                   exposure_time_pts = c(1, 2, 3),
#'                   outcome = "D.3",
#'                   balance_thresh = c(0.05, 0.1),
#'                   imp_conf = "B2")
#' c <- calcBalStats(data = test,
#'                   formulas = f,
#'                   exposure = "A",
#'                   exposure_time_pts = c(1, 2, 3),
#'                   outcome = "D.3",
#'                   balance_thresh = 0.1,
#'                   weights = w[[1]])

calcBalStats <- function(home_dir = NA, data, formulas, exposure, exposure_time_pts, outcome, balance_thresh, k = 0, weights = NULL, imp_conf = NULL, verbose = TRUE, save.out = TRUE){
  # library(cobalt)

  if(!inherits(formulas, "list")){
    stop("Please provide a list of formulas for each exposure time point", call. = FALSE)
  }
  if (!is.null(weights) & !inherits(weights, "weightitMSM")){
    stop("Please supply a list of weights output from the createWeights function (via WeightIt::WeightItMSM).", call. = FALSE)
  }

  form_name <- sapply(strsplit(names(formulas[1]), "_form"), "[", 1)
  exposure_type <- ifelse(inherits(data[, paste0(exposure, '.', exposure_time_pts[1])], "numeric"), "continuous", "binary")
  weighted = ifelse(!is.null(weights), 1, 0)

  factor_covariates <- colnames(data)[which(sapply(data, class) == "factor")]
  factor_covariates <- factor_covariates[!factor_covariates %in% "ID"]

  if (weighted == 1){
    weights_method = weights$method
    w <- weights$weights #IPTW weights
    data <- data %>% dplyr::mutate(weights = as.numeric(w))
  }
  else{
    weights_method <- "no weights"
  }


  folder <- ifelse(weighted == 0, "prebalance/", "weighted/")

  data_type <- ifelse(k == 0, "single", "imputed")

  if (data_type == "imputed" & verbose){
    cat(paste0("**Imputation ", k, "**"), "\n")
  }


  #creating initial data frames
  #data frame with all sampling weights for all exposures at all exposure time points for all histories
  all_prop_weights <- data.frame("ID" = NA,
                                 exposure = NA,
                                 exp_time = NA,
                                 history = NA)
  all_bal_stats <- data.frame()

  #cycles through user-specified exposure time points
  for (z in seq_len(length(exposure_time_pts))){
    exposure_time_pt <- exposure_time_pts[z]
    lagged_time_pts <- exposure_time_pts[exposure_time_pts<exposure_time_pt]

    # GETS COVARIATES FROM FORM FOR ASSESSING BALANCE
    # full_form <- formulas[[names(formulas)[grepl(paste0("-", exposure_time_pt),
    #                                              names(formulas))]]]
    full_form <- formulas[[names(formulas)[as.numeric(sapply(strsplit(names(formulas), "-"), "[",2)) == exposure_time_pt]]]

    covars <- paste(deparse(full_form[[3]], width.cutoff = 500), collapse = "") # gets covariates
    covar_time <- sapply(strsplit(unlist(strsplit(as.character(covars), "\\+")), "\\."), "[", 2)
    covars <- as.character(unlist(strsplit(covars, "\\+")))
    covars <- gsub(" ", "", covars)


    # GETTING BALANCE STATS FOR T=1 W/ NO EXPOSURE HISTORY (ok to standardize immediately)
    if (length(lagged_time_pts) == 0) {
      temp <- data

      # Unweighted pre-balance checking
      if (!weighted) {
        if (exposure_type == "continuous") {
          bal_stats <- cobalt::col_w_cov(temp[, c(covars)], temp[, paste0(exposure, ".",
                                                                          exposure_time_pt)], std = TRUE) # finding correlation
        }
        else if (exposure_type == "binary") {
          bal_stats <- cobalt::col_w_smd(temp[, c(covars)], temp[, paste0(exposure, ".",
                                                                          exposure_time_pt)], std = TRUE) # finding smd
        }
      }

      # Weighted balance checking
      else if (weighted) {
        if (exposure_type == "continuous") {
          bal_stats <- cobalt::col_w_cov(temp[, c(covars)], temp[, paste0(exposure, ".",
                                                                          exposure_time_pt)], std = TRUE, # finding cor
                                         weights = temp[, "weights"])
        }
        else if (exposure_type == "binary") {
          bal_stats <- cobalt::col_w_smd(temp[, c(covars)], temp[, paste0(exposure, ".",
                                                                          exposure_time_pt)], std = TRUE, # finding smd
                                         weights = temp[, "weights"])
        }
      }
      bal_stats <- as.data.frame(bal_stats)
      colnames(bal_stats) <- "std_bal_stats"
    } #ends lag=0

    # ASSIGNING HISTORIES FOR EXP TIME POINTS T>1
    if (length(lagged_time_pts) > 0) {
      # creating proportion weights based on proportion of individuals in a given exposure history
      prop_weights <- data.frame(id = data[, "ID"],
                                 exposure = exposure,
                                 exp_time = exposure_time_pt,
                                 history = NA)
      colnames(prop_weights)[colnames(prop_weights) == "id"] <- "ID"

      # finding histories up until exp time point T
      histories <- apply(gtools::permutations(2, length(lagged_time_pts), c(1, 0), repeats.allowed = TRUE),
                         1, paste, sep = "", collapse = "-")
      histories <- as.data.frame(histories)

      # cycle through histories to identify individuals with each history
      for (h in seq_len(nrow(histories))) {
        his <- as.data.frame(histories[h, ])
        his <- as.character(unlist(his[1, ]))

        # get wide data variable names for exposures at each time in each history
        exps_time <- apply(expand.grid(exposure, as.character(lagged_time_pts)), 1, paste, collapse = ".", sep = ".")

        # initiate flag marking whether each individual falls into the given history at 0 (t(1)-1)
        data$flag <- 0

        # cycles through times in each history and sequentially flags individuals who meet each criterion
        # (e.g., hi at 6, lo at 15, etc.) for that history
        for (t in seq_len(length(lagged_time_pts))) {
          time <- lagged_time_pts[t]
          exp <- as.numeric(sapply(strsplit(as.character(his), "-"), "[", t)) # hi/lo indicator
          flag <- t - 1

          if (exp == 0) { # low levels/absent
            if (exposure_type == "continuous") {
              data$flag <- ifelse(data[, exps_time[t]] <= median(data[, paste0(exposure, ".", exposure_time_pt)])
                                  & data$flag == flag, t, NA) # finding those w/ vals <= median exp @ time pt & flagged at prev t's
            }
            else { # for binary exp
              data$flag <- ifelse(data[, exps_time[t]] == 0 & data$flag == flag, t, NA) # if exposure is absent & flagged at prev t's
            }

          }

          if (exp == 1) { # hi levels/present
            if (exposure_type == "continuous") {
              data$flag <- ifelse(data[, exps_time[t]] > median(data[, paste0(exposure, ".", exposure_time_pt)])
                                  & data$flag == flag, t, NA) # finding those w/ vals > median exp @ time pt & flagged at prev t's
            }
            else { # binary exp
              data$flag <- ifelse(data[, exps_time[t]] == 1 & data$flag == flag, t, NA) # if exposure is present & flagged at prev t's
            }

          }

        } # ends history's constituent time pts time points (e.g., 6, 15, 24)

        # finding ids who met criteria for that history
        ids <- data[data$flag == t, "ID"]

        # labels those ids w/ that history
        prop_weights[prop_weights[, "ID"] %in% ids, "history"] <- paste(his, collapse = ",")
        data$flag <- NULL # resets flag
      } # ends history loop (e.g., "l-l-l")

      #summarize contributors to each history
      prop_sum <- prop_weights %>%
        dplyr::group_by(as.factor(history)) %>%
        dplyr::summarize(n = dplyr::n())

      # GET BALANCE STATISTICS FOR T>1 (when there is a history to weight on)
      if (length(lagged_time_pts) > 0) {
        # Merge data with prop_weights by ID
        temp <- merge(data, prop_weights, by = "ID", all.x = TRUE)

        # Removing any histories that only have 1 or 0 person contributing (cannot calc bal stats)
        if (sum(prop_sum$n == 1) > 0 | sum(prop_sum$n == 0) > 0) {
          ommitted_histories <- as.character(as.data.frame(prop_sum)[prop_sum$n == 1, 1])

          if (data_type == "imputed"){
            cat(paste0("USER ALERT: the following history/histories, ", ommitted_histories,
                       ", has/have been omitted from balance checking for exposure ", exposure,
                       ", imputation ", k, ", at time point ", exposure_time_pt))
          }
          else{
            cat(paste0("USER ALERT: the following history/histories, ", ommitted_histories,
                       ", has/have been omitted from balance checking for exposure ", exposure,
                       " at time point ", exposure_time_pt))
          }

          temp <- temp[!temp$history %in% ommitted_histories, ]
        } #ends hist exc

        # Unweighted pre-balance checking
        if (weighted == 0) { # no IPTW weighting but weighting on history
          if (exposure_type == "continuous") {
            bal_stats <- sapply(sort(unique(temp$history)), function(i) { # finding balance by history
              temp2 <- temp %>%
                dplyr::filter(history == i)

              cobalt::col_w_cov(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std = FALSE, # finding covariance
                                subset = temp2$history[temp2$history == i] == i) # subsetting by that history
            })

            # getting weighted mean across histories (cols), weighting by proportion of those w/ that same history
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ], tabulate(as.factor(temp$history)) / nrow(temp))
            })

            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))

            # standardizing balance statistics after weighting by history
            bal_stats <- bal_stats %>%
              dplyr::mutate(std_bal_stats = weighted_bal_stats /
                              (sapply(seq(ncol(data[, covars])), function(x) {
                                sd(as.numeric(data[, covars][, x]), na.rm = TRUE) }) *# unweighted covar sd
                                 sd(data[, paste0(exposure, ".", exposure_time_pt)], na.rm = TRUE)))  # exposure SD at that time pt

            bal_stats <- bal_stats %>%
              dplyr::select(contains(c("std")))
          } #ends continuous

          else if (exposure_type == "binary") {
            bal_stats <- sapply(sort(unique(temp$history)), function(i) {
              temp2 <- temp %>%
                dplyr::filter(history == i)

              cobalt::col_w_smd(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std = FALSE, # finding mean difference
                                subset = temp2$history[temp2$history == i] == i) # subsetting by that history
            })

            # getting weighted mean across histories, weighting by proportion of those w/ that same history
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ], tabulate(as.factor(temp$history)) / nrow(temp))
            })

            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))
            # standardizing balance statistics after finding weighted balance stats

            bal_stats <- bal_stats %>%
              dplyr::mutate(std_bal_stats = weighted_bal_stats/sapply(seq(ncol(data[, covars])), function(x){
                sqrt(mean( #dividing by pool SD estimate (unadjusted)
                  var(as.numeric(data[paste0(exposure, ".", exposure_time_pts[1]) == 1, covars[x]])), #treated var
                  var(as.numeric(data[paste0(exposure, ".", exposure_time_pts[1]) == 0, covars[x]]))))})) #untreated var
            bal_stats <- bal_stats %>%
              dplyr::select(contains(c("std")))
          } #ends binary
        } #ends weighted=0


        # Weighted balance checking
        else if (weighted == 1) { # if weighted, use IPTW weights from weightitmsm and weight by history

          if (exposure_type == "continuous") {
            # finds balance for each covariate clustered/subset by history
            bal_stats <- sapply(sort(unique(temp$history)), function(i) {
              temp2 <- temp %>%
                dplyr::filter(history == i)

              cobalt::col_w_cov(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std = FALSE, # finding covariance
                                subset = temp2$history[temp2$history == i] == i, # subsetting by that history
                                weights = temp2[, "weights"]) # adding IPTW weights
            })

            # getting weighted mean across histories, weighting by proportion of those w/ that same history
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ], tabulate(as.factor(temp$history)) / nrow(temp))
            })

            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))

            # standardizing balance statistics after weighting by history
            bal_stats <- bal_stats %>%
              dplyr::mutate(std_bal_stats = weighted_bal_stats /
                              (sapply(seq(ncol(data[, covars])), function(x) {
                                sd(as.numeric(data[, covars][, x]), na.rm = TRUE) # unweighted covar sd
                              }) * sd(data[, paste0(exposure, ".", exposure_time_pt)], na.rm = TRUE)))

            # For a weighted_bal_stat of 0, make std stat also 0 so as not to throw an error
            bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] <- 0

            bal_stats <- bal_stats %>%
              dplyr::select(contains(c("std")))
          } #ends continuous

          else if (exposure_type == "binary") {
            # finds balance for each covariate clustered/subset by history
            bal_stats <- sapply(sort(unique(temp$history)), function(i) {
              temp2 <- temp %>%
                dplyr::filter(history == i)

              cobalt::col_w_smd(temp2[, c(covars)], temp2[, paste0(exposure, ".", exposure_time_pt)], std = FALSE, # finding mean difference
                                subset = temp2$history[temp2$history == i] == i, # subsetting by that history
                                weights = temp2[, "weights"]) # adding IPTW weights
            })

            # getting weighted mean across histories, weighting by proportion of those w/ that same history
            weighted_bal_stats <- sapply(seq(nrow(bal_stats)), function(i) {
              weighted.mean(t(bal_stats[i, ])[1, ], tabulate(as.factor(temp$history)) / nrow(temp))
            })

            bal_stats <- as.data.frame(cbind(bal_stats, weighted_bal_stats))
            # standardizing balance statistics after finding weighted balance stats

            bal_stats <- bal_stats %>%
              dplyr::mutate(std_bal_stats=weighted_bal_stats/ sapply(seq(ncol(data[, covars])), function(x){
                sqrt(mean( #dividing by pool SD estimate (unadjusted)
                  var(as.numeric(data[paste0(exposure, ".",
                                             exposure_time_pts[1]) == 1, covars[x]])), #treated var
                  var(as.numeric(data[paste0(exposure, ".",
                                             exposure_time_pts[1]) == 0, covars[x]]))))}))  #untreated var

            # For a weighted_bal_stat of 0, make std stat also 0 so as not to throw an error
            bal_stats$std_bal_stats[is.nan(bal_stats$std_bal_stats)] <- 0
            bal_stats <- bal_stats %>%
              dplyr::select(contains(c("std")))
          } #ends binary
        } #ends weighted

        # collect proportions for histories at this time point
        all_prop_weights <- rbind(all_prop_weights, prop_weights)
      }
    } # ends lag>0 loops


    # ADDS INFO TO BAL STATS
    bal_stats <- as.data.frame(bal_stats)
    bal_stats$covariate <- rownames(bal_stats)

    # Renames factor covariates
    bal_stats$covariate[sapply(strsplit(bal_stats$covariate, "_"), "[", 1) %in% factor_covariates] <-
      sapply(strsplit(bal_stats$covariate, "_"), "[", 1)[sapply(strsplit(bal_stats$covariate, "_"), "[", 1) %in% factor_covariates]

    #adds custom bal thresh info
    if (!is.null(imp_conf)){
      bal_stats <- bal_stats %>%
        dplyr::mutate(bal_thresh = ifelse(bal_stats$covariate %in% imp_conf, balance_thresh[1], balance_thresh[2]))
      bal_stats <- bal_stats %>%
        dplyr:: mutate(balanced = ifelse(abs(bal_stats$std_bal_stats) < bal_stats$bal_thresh, 1, 0) )
    }
    else{
      bal_stats <- bal_stats %>%
        dplyr::mutate(bal_thresh = balance_thresh)
      bal_stats <- bal_stats %>%
        dplyr:: mutate(balanced = ifelse(abs(bal_stats$std_bal_stats) < bal_stats$bal_thresh, 1, 0) )
    }

    bal_stats <- bal_stats %>%
      dplyr::mutate(exposure = exposure, exp_time = exposure_time_pt) %>%
      dplyr::mutate(covar_time = sapply(strsplit(covariate, "\\."), "[", 2))

    all_bal_stats <- rbind(all_bal_stats, bal_stats)
    all_bal_stats$covar_time[is.na(all_bal_stats$covar_time)] <- 0

    # Make love plot per exposure time point
    make_love_plot(home_dir, folder, exposure, exposure_time_pt, exposure_type, k,
                   form_name, bal_stats, data_type, balance_thresh, weights_method, imp_conf, verbose, save.out)

  }     # Ends exp_time_pt loop


  if (verbose & save.out) {
    if (data_type == "imputed"){
      cat(paste0("For each time point and imputation, ", gsub("/", "", folder), " summary plots for ", form_name, " formulas weighting method ",
                 weights_method, " have now been saved in the '", folder, "plots/' folder."), "\n")
    }
    else {cat(paste0(" For each time point, ", gsub("/", "", folder), " summary plots for ", form_name, " formulas and weighting method ",
                     weights_method, " have now been saved in the '", folder, "plots/' folder."), "\n")
    }
  }


  # Summarizing balance
  bal_summary_exp <- all_bal_stats %>%
    dplyr::group_by(exp_time) %>%
    dplyr::summarize(balanced_n = sum(balanced == 1), # Tallying balanced covariates
                     imbalanced_n = sum(balanced == 0), # Tallying imbalanced covariates
                     n = dplyr::n())


  if (save.out){
    write.csv(bal_summary_exp, paste0(home_dir, "/balance/", folder, form_name, "_", exposure, "_", k, "_",
                                      weights_method, "_balance_stat_summary.csv"))
    write.csv(all_prop_weights, paste0(home_dir, "/balance/", folder, "/", form_name, "_form_", exposure, "_", k,
                                       "_", weights_method, "_history_sample_weight.csv"))
    if (verbose){
      cat("\n")
      if (data_type == "imputed"){
        cat(paste0("Balance statistics using ", form_name, " formulas for ", exposure, ", imputation ", k, ", using ",
                   weights_method, " have been saved in the 'balance/", folder, "' folder"), "\n")

        cat(paste0("Sampling weights ", "using the ", form_name, " for ", exposure, ", imputation ", k,
                   " have been saved in the 'balance/", folder, "' folder"), "\n")
      }
      else{
        cat(paste0("Balance statistics using ", form_name, " formulas for ", exposure, "using ",
                   weights_method, " have been saved in the 'balance/", folder, "' folder"), "\n")

        cat(paste0("Sampling weights ", "using the ", form_name, " for ", exposure,
                   " have been saved in the 'balance/", folder, "' folder"), "\n")
      }
      cat("\n")
    }
  }

  # tallies total possible COVARIATES FROM FORM FOR ASSESSING BALANCE
  all_form <- as.data.frame(do.call(rbind, formulas))
  tot_covars <- deparse(all_form[, 3], width.cutoff = 300)
  tot_covars <- as.character(unlist(strsplit(tot_covars, "\\+")))[!grepl("form", as.character(unlist(strsplit(tot_covars, "\\+"))))]
  tot_covars <- gsub(" ", "", tot_covars)
  tot_covars <- na.omit(sapply(strsplit(tot_covars, "\\."), "[", 1)[!duplicated(sapply(strsplit(tot_covars, "\\."), "[", 1))])

  imbalanced_covars <- sum(bal_summary_exp$imbalanced_n, na.rm = TRUE)
  total_covars <- sum(bal_summary_exp$n, na.rm = TRUE)
  percentage_imbalanced <- round((imbalanced_covars / total_covars) * 100, 0)

  remaining_imbalanced_domains <- length(sapply(strsplit(all_bal_stats[all_bal_stats$balanced == 0, "covariate"], "\\."),
                                                "[", 1)[!duplicated(sapply(strsplit(all_bal_stats[all_bal_stats$balanced == 0, "covariate"],
                                                                                    "\\."), "[", 1))])
  total_domains <- length(tot_covars)

  # remaining_avg_abs_corr <- round(mean(abs(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"]), na.rm = TRUE), 2)
  remaining_avg_abs_corr <- round(median(abs(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"]), na.rm = TRUE), 2)
  remaining_corr_range <- paste0(round(min(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"], na.rm = TRUE), 2),
                                 "-", round(max(all_bal_stats[all_bal_stats$balanced == 0, "std_bal_stats"], na.rm = TRUE), 2))



  if (data_type == "imputed"){

    cat(paste0("USER ALERT: For exposure ", exposure, " imputation ", k, " using ", weights_method, " and ", form_name, " formulas: "), "\n")
    cat(paste0("The median absolute value relation between exposure and confounder is ", round(median(abs(all_bal_stats$std_bal_stats)), 2), " (range = ",
               round(min(all_bal_stats$std_bal_stats), 2), "-", round(max(all_bal_stats$std_bal_stats), 2), ")."), "\n")

    cat(paste0("As shown below, ", imbalanced_covars, " out of ", total_covars, " (", percentage_imbalanced,
               "%) covariates across time points, corresponding to ",
               remaining_imbalanced_domains, " out of ", total_domains,
               " domains, remain imbalanced with a remaining median absolute value correlation/std mean difference of ",
               remaining_avg_abs_corr, " (range= ", remaining_corr_range, "):"), "\n")
    cat("\n")
    cat(knitr::kable(bal_summary_exp, caption = paste0("Imbalanced Covariates for imputation ", k, " using ",
                                                       weights_method, " and ", form_name, " formulas"), format = 'pipe'), sep = "\n")
    cat("\n")
    cat("\n")
  } else {
    cat(paste0("As shown below, for exposure ", exposure, " using ", weights_method, ", and ", form_name, " formulas, ",
               imbalanced_covars, " out of ", total_covars, " (", percentage_imbalanced, "%) covariates across time points corresponding to ",
               remaining_imbalanced_domains, " out of ", total_domains,
               " domains remain imbalanced with a remaining average absolute value correlation/std mean difference of ",
               remaining_avg_abs_corr, " (range= ", remaining_corr_range, ") :"), "\n")
    cat("\n")
    cat(knitr::kable(bal_summary_exp, caption = paste0("Imbalanced covariates using ",
                                                       weights_method, " and ", form_name, " formulas"), format = 'pipe'), sep = "\n")
    cat("\n")
    cat("\n")

  }

  rownames(all_bal_stats) <- NULL

  all_bal_stats
}
