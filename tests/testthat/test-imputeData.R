
run_code_imp <- function(){
  quiet <- function(x) { #quiets cat output
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")

  imputed_datasets <- quiet(suppressWarnings(imputeData(object, data_to_impute, read_imps_from_file="no")))
  save(list=c("imputed_datasets"),
       file="/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
}


test_that("imputed data is length of original data", {

  run_code_imp()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")

  expect_true(nrow(imputed_datasets[[1]]) ==nrow(data))

  # expect_match(names(imputed_datasets[[1]]), covariates_to_include[[1]]$row, ignore.case=T)
})

test_that("confounders, exposures, outcomes in imputed dataset", {

  run_code_imp()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")

  expect_true(sum(sapply(strsplit(covariates_to_include[[1]]$row, "\\."), "[",1) %in% names(imputed_datasets[[1]]))==nrow(covariates_to_include[[1]]))
  expect_true(sum(object$exposures %in% names(imputed_datasets[[1]]))==length(object$exposures))
  expect_true(sum(object$outcomes %in% names(imputed_datasets[[1]]))==length(object$outcomes))


})
