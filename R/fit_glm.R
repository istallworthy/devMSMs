# Internal call
fit_glm <- function(
  data, weights, outcome, exposure, epoch,
  model, int_order, family, covariates, sep
) {

  # Create epoch-averages of exposure variables
  epoch_vars <- exposure
  if (any(exposure != epoch)) {
    epoch_df <- .generate_epoch_exposures(data, exposure, epoch, sep)
    epoch_vars <- colnames(epoch_df)
    missing_columns <- !(colnames(epoch_df) %in% colnames(data))
    data <- cbind(data, epoch_df[, missing_columns])
  }
  if (model == "m2" || model == "m3") {
    interactions <- .get_interaction_terms(epoch_vars, int_order)
  }

  # H_0: intercept only
  if (model == "int") {
    covs <- "1"
    # H_0: covariates only
  } else if (model == "covs") {
    covs <- covariates
    # H_A: main effects only
  } else if (model == "m0") {
    covs <- epoch_vars
    # H_A: main effects + covariates
  } else if (model == "m1") {
    covs <- c(epoch_vars, covariates)
    # H_A: main effects + interactions
  } else if (model == "m2") {
    covs <- c(epoch_vars, interactions)
    # H_A: main effects + covariates + interactions
  } else if (model == "m3") {
    covs <- c(epoch_vars, interactions, covariates)
  }

  model <- WeightIt::glm_weightit(
    formula = reformulate(covs, response = outcome),
    data = data,
    family = family,
    weightit = weights
  )
  return(model)
}
