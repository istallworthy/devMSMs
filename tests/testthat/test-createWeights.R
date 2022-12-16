run_code_weights <- function(){
  quiet <- function(x) { #quiets cat output
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")

  weights_models <- quiet(createWeights(object, wide_long_datasets, forms))


  save(list=c("weights_models"),
       file="/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/weights_models.Rdata")
}


test_that("there are weights for each imputation, exposure time point, and exposure-outcome pair ", {
  run_code_weights()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/weights_models.Rdata")

  expect_true(
    length(weights_models)==nrow(expand.grid(object$exposures, object$exposure_time_pts, object$m, object$outcomes))
  )

})


test_that("weights are not NAs", {

  run_code_weights()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
  # load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/weights_models.Rdata")

  expect_false(
    sum(as.logical(unlist(lapply(weights_models, function(x) is.na(x$weights)))))>0
  )

})
