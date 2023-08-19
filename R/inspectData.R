
inspectData <-function(data, home_dir, exposure, outcome, tv_confounders, ti_confounders, epochs = NULL, hi_lo_cut = NULL){

  # Error checking
  if (!class(data) %in% c("mids", "data.frame", "character")) {
    stop("Please provide either a 'mids' object, a data frame, or a directory with imputed csv files in the 'data' field.")
  }

  ID <- "ID"
  time_invar_covars <- ti_confounders
  time_var_covars <- tv_confounders
  time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))

  if (exposure_type == "continuous"){
    if (is.null(hi_lo_cut)){
      hi_lo_cut <- c(0.75, 0.25)
    }
  }

  #
  # Load necessary packages
  library(readr)
  library(dplyr)
  library(tidyr)

  if (class(data) == "mids"){
    data <-as.data.frame(mice::complete(data,1))
  }


  if (class(data) == "character") {
    if (!dir.exists(data)) {
      stop("Please provide a valid directory path with imputed datasets, a data frame, or a 'mids' object for the 'data' field.")
    }

    # List imputed files
    files <- list.files(data, full.names = TRUE, pattern = "\\.csv")

    # Read and process 1st imp
    data <- read.csv(files[1])
  }


  # long format to wide
  if("WAVE" %in% colnames(data)){
    v <- sapply(strsplit(tv_confounders, "\\."), "[",1)
    v <- v[!duplicated(v)]
    data_wide <- stats::reshape(data = data_long, idvar = "ID", v.names = v, timevar = "WAVE",
                                direction = "wide")

    #removing all NA cols (i.e., when data were not collected)
    data_wide <- data_wide[,colSums(is.na(data_wide))<nrow(data_wide)]
    data <- data_wide
  }

  exposure_type <- ifelse(class(data[, paste0(exposure, '.', exposure_time_pts[1])]) == "numeric", "continuous", "binary")

  # Exposure summary
  exposure_summary <- data %>%
    dplyr:: select(colnames(data)[grepl(exposure, colnames(data))])
  exposure_summary <- psych::describe(exposure_summary, fast = TRUE)



  cat(knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"), format = 'pipe'), sep = "\n")
  knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"), format = 'html') %>%
    kableExtra::kable_styling() %>%
    kableExtra::save_kable(file = file.path(home_dir, paste0("/", exposure, "_exposure_info.html")))

  cat(paste0(exposure, " exposure descriptive statistics have now been saved in the home directory"), "\n")
  cat("\n")


  # Outcome summary
  outcome_summary <- data %>%
    dplyr:: select(contains(sapply(strsplit(outcome, "\\."), "[",1)))
  outcome_summary <- psych::describe(outcome_summary, fast = TRUE)


  cat(knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ",
                                                     sapply(strsplit(outcome, "\\."), "[",1), " Information"), format = 'pipe'), sep = "\n")
  knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ",
                                                 sapply(strsplit(outcome, "\\."), "[",1), " Information"), format = 'html') %>%
    kableExtra::kable_styling() %>%
    kableExtra::save_kable(file = file.path(home_dir, paste0("/", sapply(strsplit(outcome, "\\."), "[",1), "_outcome_info.html")))

  cat(paste0(sapply(strsplit(outcome, "\\."), "[",1), " outcome descriptive statistics have now been saved in the home directory"), "\n")



  # Confounder summary
  potential_covariates <- colnames(data)[!(colnames(data) %in% c(ID))]

  if (sum(tv_confounders %in% potential_covariates) != length(tv_confounders)){
    stop(paste(tv_confounders[!tv_confounders %in% potential_covariates]), " time-varying confounders are not present in the dataset.")
  }

  if (sum(ti_confounders %in% potential_covariates) != length(ti_confounders)){
    stop(paste(ti_confounders[!ti_confounders %in% potential_covariates]), " time invariant confounders are not present in the dataset.")
  }

  all_potential_covariates <- c(time_invar_covars, time_var_covars)
  # all_potential_covariates <- all_potential_covariates[!(all_potential_covariates %in% c(paste(outcome, outcome_time_pt, sep = "."), time_var_exclude))]
  all_potential_covariates <- all_potential_covariates[order(all_potential_covariates)]

  # Format for table output to visualize available covariates by time point
  covar_table <- data.frame(variable = sapply(strsplit(all_potential_covariates, "\\."), "[", 1),
                            time_pt = sapply(strsplit(all_potential_covariates, "\\."), "[", 2)) %>%
    dplyr::arrange(time_pt, variable) %>%
    dplyr::group_by(time_pt) %>%
    dplyr::summarize(variable = toString(variable))

  write.csv(covar_table, glue::glue("{home_dir}/balance/{exposure}-{outcome}_covariates_considered_by_time_pt.csv"), row.names = FALSE)

  unique_vars <- length(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))

  test <- data.frame(matrix(nrow = length(time_pts), ncol = unique_vars))
  colnames(test) <- unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."),
                                                       "[", 1)))[order(unique(c(time_invar_covars,
                                                                                sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))]
  rownames(test) <- time_pts

  for (l in 1:nrow(test)) {
    test[l, c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]), all_potential_covariates)], "\\."), "[", 1), time_invar_covars)] <- 1
  }

  test <- test[, colnames(test)[!(colnames(test) %in% c(ID))]]
  NumTimePts <- data.frame(NumTimePts = colSums(test, na.rm = TRUE))
  test <- rbind(test, t(NumTimePts))
  NumVars <- data.frame(NumVars = rowSums(test, na.rm = TRUE))
  test[1:nrow(test), ncol(test) + 1] <- NumVars
  write.csv(test, glue::glue("{home_dir}/balance/{exposure}-{outcome}_matrix_of_covariates_considered_by_time_pt.csv"), row.names = TRUE)

  message(glue::glue("See the 'balance/' folder for a table and matrix displaying all covariates confounders considered at each exposure time point for {exposure} and {outcome}."), "\n")
  cat("\n")

  #-2 to exclude ID and WAVE
  message(glue::glue("USER ALERT: Below are the {as.character(length(all_potential_covariates) - 2)} variables spanning {unique_vars - 2} unique domains that will be treated as confounding variables for the relation between {exposure} and {outcome}."), "\n",
          "Please inspect this list carefully. It should include all time-varying covariates, time invariant covariates, as well as lagged levels of exposure and outcome variables if they were collected at time points earlier than the outcome time point.", "\n")
  cat("\n")
  print(all_potential_covariates[!(all_potential_covariates %in% c(ID))])

  #covariate correlations
  covariates_to_include <- all_potential_covariates
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    tibble(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor = cormat[ut],
      p = pmat[ut]
    )
  }

  # Creates final dataset with only relevant variables
  covariates_to_include <- covariates_to_include[order(covariates_to_include)]
  variables_to_include <- unique(c(ID,  outcome, covariates_to_include, time_var_covars))
  data2 <- data %>%
    select(all_of(variables_to_include))
  data_to_impute <- data2

  # Makes correlation table
  corr_matrix <- cor(as.data.frame(lapply(data_to_impute[, colnames(data_to_impute) != ID], as.numeric)), use = "pairwise.complete.obs")
  ggcorrplot::ggcorrplot(corr_matrix,  type = "lower")+
    ggplot2::theme(axis.text.x = element_text(size = 5, margin = ggplot2::margin(-2,0,0,0)),  # Order: top, right, bottom, left
          axis.text.y = element_text(size = 5, margin = ggplot2::margin(0,-2,0,0))) +
    ggplot2::geom_vline(xintercept = 1:ncol(mtcars) - 0.5, colour="white", size = 2) +
    ggplot2::geom_hline(yintercept = 1:ncol(mtcars) - 0.5, colour="white", size = 2)

  # Save correlation plot
  pdf(file = paste0(home_dir, "/", exposure, "-", outcome, "_all_vars_corr_plot.pdf"))
  print(ggplot2::last_plot())
  dev.off()

  cat("\n")
  cat("A correlation plot of all variables in the dataset has been saved in the home directory", "\n")
  cat("\n")



  # Exposure history summary
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(time_pts),
                         values = time_pts)
  }




}
