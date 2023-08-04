
inspectConfounders <-function(data, home_dir, exposure, outcome, tv_confounders, ti_confounders){
  ID <- "ID"
  time_invar_covars <- ti_confounders
  time_var_covars <- tv_confounders

  #
  # Load necessary packages
  library(readr)
  library(dplyr)
  library(tidyr)

  # Identifying potential covariates (time-varying and time invariant)
  potential_covariates <- colnames(data)[!(colnames(data) %in% c(ID, "WAVE"))]


  all_potential_covariates <- c(time_invar_covars, time_var_covars)
  # all_potential_covariates <- all_potential_covariates[!(all_potential_covariates %in% c(paste(outcome, outcome_time_pt, sep = "."), time_var_exclude))]
  all_potential_covariates <- all_potential_covariates[order(all_potential_covariates)]

  # Format for table output to visualize available covariates by time point
  covar_table <- data.frame(variable = sapply(strsplit(all_potential_covariates, "\\."), "[", 1),
                            time_pt = sapply(strsplit(all_potential_covariates, "\\."), "[", 2)) %>%
    arrange(time_pt, variable) %>%
    group_by(time_pt) %>%
    summarize(variable = toString(variable))

  write_csv(covar_table, glue::glue("{home_dir}balance/{exposure}-{outcome}_covariates_considered_by_time_pt.csv"), row.names = FALSE)

  unique_vars <- length(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))

  test <- data.frame(matrix(nrow = length(time_pts), ncol = unique_vars))
  colnames(test) <- unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."),
                                                       "[", 1)))[order(unique(c(time_invar_covars,
                                                                                sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))]
  rownames(test) <- time_pts

  for (l in 1:nrow(test)) {
    test[l, c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]), all_potential_covariates)], "\\."), "[", 1), time_invar_covars)] <- 1
  }

  test <- test[, colnames(test)[!(colnames(test) %in% c(ID, "WAVE"))]]
  NumTimePts <- data.frame(NumTimePts = colSums(test, na.rm = TRUE))
  test <- rbind(test, t(NumTimePts))
  NumVars <- data.frame(NumVars = rowSums(test, na.rm = TRUE))
  test[1:nrow(test), ncol(test) + 1] <- NumVars
  write_csv(test, glue::glue("{home_dir}balance/{exposure}-{outcome}_matrix_of_covariates_considered_by_time_pt.csv"), row.names = TRUE)

  cat(glue::glue("See the 'balance/' folder for a table and matrix displaying all covariates confounders considered at each exposure time point for {exposure} and {outcome}."), "\n")
  cat("\n")

  #-2 to exclude ID and WAVE
  cat(glue::glue("USER ALERT: Below are the {as.character(length(all_potential_covariates) - 2)} variables spanning {unique_vars - 2} unique domains that will be treated as confounding variables for the relation between {exposure} and {outcome}."), "\n",
      "Please inspect this list carefully. It should include all time-varying covariates, time invariant covariates, as well as lagged levels of exposure and outcome variables if they were collected at time points earlier than the outcome time point.", "\n")
  cat("\n")
  print(all_potential_covariates[!(all_potential_covariates %in% c(ID, "WAVE"))])


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
  variables_to_include <- unique(c(ID, "WAVE", exposure, outcome, covariates_to_include, time_var_covars))
  data2 <- data %>%
    select(all_of(variables_to_include))
  data_to_impute <- data2

  # Makes correlation table
  corr_matrix <- cor(as.data.frame(lapply(data_to_impute[, colnames(data_to_impute) != ID], as.numeric)), use = "pairwise.complete.obs")
  ggcorrplot::ggcorrplot(corr_matrix, method = "color", order = 'alphabet', diag = FALSE, type = "lower", tl.cex = 0.5, tl.col = "black")

  # Save correlation plot
  pdf(file = paste0(home_dir, exposure, "-", outcome, "_all_vars_corr_plot.pdf"))
  print(ggcorrplot::last_plot())
  dev.off()

  cat("A correlation plot of all variables in the dataset has been saved in the home directory", "\n")
  cat("\n")

  # Remove exposure/outcome from covariate correlation assessment
  data2 <- data2 %>%
    select(-c(exposure, outcome))

  # Inspect correlations among covariates to check for redundancy and opportunities to simplify the model
  hi_corr_covars <- suppressWarnings(Hmisc::rcorr(as.matrix(data2[, 3:ncol(data2)])))
  hi_corr_covars <- flattenCorrMatrix(hi_corr_covars$r, hi_corr_covars$P)

  View_hi_corr_covars <- hi_corr_covars %>%
    filter(abs(cor) > 0.7) %>%
    filter(!(row %in% "WAVE") & !(column %in% "WAVE"))

  cat("USER ALERT: To simplify the balancing models, consider removing any highly correlated:", "\n")
  cat(knitr::kable(View_hi_corr_covars, caption = "Correlated Covariates", format = 'pipe'), sep = "\n")
  cat("\n")


}
