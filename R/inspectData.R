
#' Inspect long/wide/imputed data
#'
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param home_dir (optional) path to home directory (required if save.out = TRUE)
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix
#' @param ti_confounders list of time invariant confounders
#' @param epochs (optional) data frame of exposure epoch labels and values
#' @param hi_lo_cut (optional) list of two numbers indicating quantile values
#'   that reflect high and low values, respectively, for continuous exposure
#'   (default is median split)
#' @param reference (optional)string of "-"-separated "l" and "h" values
#'   indicative of a reference exposure history to which to compare comparison,
#'   required if comparison is specified
#' @param comparison (optional)list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is specified
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return none
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
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             hi_lo_cut = c(0.8, 0.2),
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             hi_lo_cut = c(0.8, 0.2),
#'             reference = "l-l-l",
#'             comparison = "h-h-h",
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
#'                                 values = I(list(c(1, 2), c(3)))),
#'             save.out = FALSE)

inspectData <- function(data, home_dir, exposure, exposure_time_pts, outcome, tv_confounders, ti_confounders, epochs = NULL,
                        hi_lo_cut = NULL, reference = NA, comparison = NULL, verbose = TRUE, save.out = TRUE){

  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if(!dir.exists(home_dir)) {
      stop("Please provide a valid home directory path if you wish to save output locally.", call. = FALSE)
    }
  }
  if (missing(data)){
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  if (missing(exposure)){
    stop("Please supply a single exposure.", call. = FALSE)
  }
  if (missing(outcome)){
    stop("Please supply a single outcome.", call. = FALSE)
  }
  if (missing(exposure_time_pts)){
    stop("Please supply the exposure time points at which you wish to create weights.", call. = FALSE)
  }
  if (missing(tv_confounders)){
    stop("Please supply a list of time-varying confounders.", call. = FALSE)
  }
  if (missing(ti_confounders)){
    stop("Please supply a list of time invariant confounders.", call. = FALSE)
  }

  if (!mice::is.mids(data) & !is.data.frame(data) & !inherits(data, "list")) {
    stop("Please provide either a 'mids' object, a data frame, or a list of imputed csv files in the 'data' field.", call. = FALSE)
  }



  ID <- "ID"
  time_invar_covars <- ti_confounders
  time_var_covars <- tv_confounders
  time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, tv_confounders)] , "\\."), "[",2))

  if (mice::is.mids(data)){
    data <-as.data.frame(mice::complete(data,1))
  }

  else if (inherits(data, "list")) { #just inspects frist imputed dataset
    data <- data[[1]]

  }

  # long format to wide
  if("WAVE" %in% colnames(data)){
    v <- sapply(strsplit(tv_confounders, "\\."), "[", 1)
    v <- v[!duplicated(v)]
    data_wide <- stats::reshape(data = data_long, idvar = "ID", v.names = v, timevar = "WAVE",
                                direction = "wide")

    #removing all NA cols (i.e., when data were not collected)
    data_wide <- data_wide[,colSums(is.na(data_wide)) < nrow(data_wide)]
    data <- data_wide
  }


  if(!inherits(data, "data.frame")){
    warning(paste0("Your data is a ", class(data), ". Convert to data frame before running devMSMs."),
            call. = FALSE)
  }


  exposure_type <- ifelse(inherits(data[, paste0(exposure, '.', exposure_time_pts[1])],
                                   "numeric"), "continuous", "binary")


  # Confounder summary
  potential_covariates <- colnames(data)[!(colnames(data) %in% c(ID))]

  if (sum(tv_confounders %in% potential_covariates) != length(tv_confounders)){
    stop(paste(tv_confounders[!tv_confounders %in% potential_covariates]),
         " time-varying confounders are not present in the dataset.", call. = FALSE)
  }

  if (sum(ti_confounders %in% potential_covariates) != length(ti_confounders)){
    stop(paste(ti_confounders[!ti_confounders %in% potential_covariates]),
         " time invariant confounders are not present in the dataset.", call. = FALSE)
  }

  all_potential_covariates <- c(time_invar_covars, time_var_covars)
  all_potential_covariates <- all_potential_covariates[order(all_potential_covariates)]

  # Format for table output to visualize available covariates by time point
  covar_table <- data.frame(variable = sapply(strsplit(all_potential_covariates, "\\."), "[", 1),
                            time_pt = sapply(strsplit(all_potential_covariates, "\\."), "[", 2)) %>%
    dplyr::arrange(time_pt, variable) %>%
    dplyr::group_by(time_pt) %>%
    dplyr::summarize(variable = toString(variable))

  if(save.out){
    write.csv(covar_table, glue::glue("{home_dir}/{exposure}-{outcome}_covariates_considered_by_time_pt.csv"),
              row.names = FALSE)
  }

  unique_vars <- length(unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."), "[", 1))))

  test <- data.frame(matrix(nrow = length(time_pts), ncol = unique_vars))
  colnames(test) <- unique(c(time_invar_covars, sapply(strsplit(all_potential_covariates, "\\."),
                                                       "[", 1)))[order(unique(c(time_invar_covars,
                                                                                sapply(strsplit(all_potential_covariates,
                                                                                                "\\."), "[", 1))))]
  rownames(test) <- time_pts

  for (l in seq_len(nrow(test))) {
    z = c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]),
                                                         all_potential_covariates)], "\\."), "[", 1), time_invar_covars)
    z = z[!duplicated(z)]
    test[l, z ] <- 1
  }

  test <- test[, colnames(test)[!(colnames(test) %in% c(ID))]]
  NumTimePts <- data.frame(NumTimePts = colSums(test, na.rm = TRUE))
  test <- rbind(test, t(NumTimePts))
  NumVars <- data.frame(NumVars = rowSums(test, na.rm = TRUE))
  test[seq_len(nrow(test)), ncol(test) + 1] <- NumVars

  if(save.out){
    write.csv(test, glue::glue("{home_dir}/{exposure}-{outcome}_matrix_of_covariates_considered_by_time_pt.csv"),
              row.names = TRUE)

    if(verbose){
      print(glue::glue("See the home directory for a table and matrix displaying all covariates confounders considered at each exposure time point for {exposure} and {outcome}."), "\n")

      #-2 to exclude ID and WAVE
      print(glue::glue("USER ALERT: Below are the {as.character(length(all_potential_covariates) - 2)} variables spanning {unique_vars - 2} unique domains that will be treated as confounding variables for the relation between {exposure} and {outcome}."),
            "Please inspect this list carefully. It should include all time-varying covariates, time invariant covariates, as well as lagged levels of exposure and outcome variables if they were collected at time points earlier than the outcome time point.", "\n")
      print(all_potential_covariates[!(all_potential_covariates %in% c(ID))])
    }
  }


  # Data type
  if(verbose){
    cat("\n")
    cat("The following variables are designated as numeric:", "\n")
    print(paste(colnames(data)[sapply(data, class) == "numeric"], sep = ",", collapse = ", "))
    cat("\n")

    cat("The following variables are designated as factors:", "\n")
    print(paste(colnames(data)[sapply(data, class) == "factor"], sep = ",", collapse = ", "))
    cat("\n")

    #temporary warning re: factor levels
    cat("*temp: please inspect the levels of your factors below. at present, excluding ID, the code can only accept 2-level factors. set rest to numeric", "\n")
    print(sapply(data[,colnames(data)[sapply(data, class) == "factor"]], nlevels))

    oth <- data.frame(variable = names(sapply(data, class)) [!sapply(data, class) %in% c("numeric", "factor")],
                      type = sapply(data, class) [!sapply(data, class) %in% c("numeric", "factor")])
    if(nrow(oth) > 0 ){
      cat(knitr::kable(oth, caption = "Other variable types",
                       format = 'pipe'), sep = "\n")
      cat("\n")
    }

    if(sum(sapply(data, is.character)) > 0){
      warning(paste0(paste(names(data)[sapply(data, is.character)], sep = ", ", collapse = ", "),
                     " are of class character.", " The package cannot accept character variables."), call. = FALSE)
    }
  }
  #covariate correlations
  covariates_to_include <- all_potential_covariates

  # Creates final dataset with only relevant variables
  covariates_to_include <- covariates_to_include[order(covariates_to_include)]
  variables_to_include <- unique(c(ID,  outcome, covariates_to_include, time_var_covars))
  data2 <- data %>%
    select(all_of(variables_to_include))

  # Makes correlation table
  corr_matrix <- cor(as.data.frame(lapply(data2[, colnames(data2) != ID],
                                          as.numeric)), use = "pairwise.complete.obs")

  if(save.out){
    ggcorrplot::ggcorrplot(corr_matrix,  type = "lower")+
      ggplot2::theme(axis.text.x = element_text(size = 5, margin = ggplot2::margin(-2, 0, 0, 0)),  # Order: top, right, bottom, left
                     axis.text.y = element_text(size = 5, margin = ggplot2::margin(0, -2, 0, 0))) +
      ggplot2::geom_vline(xintercept = seq_len(ncol(mtcars)) - 0.5, colour="white", size = 2) +
      ggplot2::geom_hline(yintercept = seq_len(ncol(mtcars)) - 0.5, colour="white", size = 2)

    # Save correlation plot
    pdf(file = paste0(home_dir, "/", exposure, "-", outcome, "_all_vars_corr_plot.pdf"))
    print(ggplot2::last_plot())
    dev.off()

    if(verbose){
      cat("\n")
      cat("A correlation plot of all variables in the dataset has been saved in the home directory", "\n")
      cat("\n")
    }
  }


  # Exposure summary
  exposure_summary <- data %>%
    dplyr:: select(colnames(data)[grepl(exposure, colnames(data))])
  exposure_summary <- sapply(exposure_summary, as.numeric)
  exposure_summary <- psych::describe(exposure_summary, fast = TRUE)


  if (save.out){
    knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"),
                 format = 'html') %>%
      kableExtra::kable_styling() %>%
      kableExtra::save_kable(file = file.path(home_dir, paste0("/", exposure, "_exposure_info.html")))
    if(verbose){
      cat(knitr::kable(exposure_summary, caption = paste0("Summary of ", exposure, " Exposure Information"),
                       format = 'pipe'), sep = "\n")
      cat(paste0(exposure, " exposure descriptive statistics have now been saved in the home directory"), "\n")
      cat("\n")
    }
  }

  eval_hist(data = data2, exposure, tv_confounders, epochs,
            exposure_time_pts, hi_lo_cut, ref = reference, comps = comparison, verbose)

  # Exposure history summary
  if( is.null(epochs)){ #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(time_pts),
                         values = time_pts)
  }
  else{
    if( !is.data.frame(epochs) | ncol(epochs) != 2 | sum(colnames(epochs) == c("epochs", "values")) != ncol(epochs)){
      stop("If you supply epochs, please provide a dataframe with two columns of epochs and values.",
           call. = FALSE)
    }
    if(sum(is.na(epochs$values)) > 0){
      stop("Please provide one or a list of several values for each epoch.", call. = FALSE)
    }
  }

  # Outcome summary
  outcome_summary <- data[, grepl(sapply(strsplit(outcome, "\\."),
                                         "[", 1), colnames(data))]
  outcome_summary <- psych::describe(outcome_summary, fast = TRUE)

  if(save.out){
    knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ",
                                                   sapply(strsplit(outcome, "\\."), "[", 1), " Information"), format = 'html') %>%
      kableExtra::kable_styling() %>%
      kableExtra::save_kable(file = file.path(home_dir, paste0("/", sapply(strsplit(outcome, "\\."), "[", 1), "_outcome_info.html")))

    if (verbose){
      cat(knitr::kable(outcome_summary, caption = paste0("Summary of Outcome ",
                                                         sapply(strsplit(outcome, "\\."), "[", 1), " Information"),
                       format = 'pipe'), sep = "\n")

      cat(paste0(sapply(strsplit(outcome, "\\."), "[", 1), " outcome descriptive statistics have now been saved in the home directory"), "\n")
    }
  }
}
