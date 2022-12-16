
run_code_f_w <- function(){
  quiet <- function(x) { #quiets cat output
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")

  wide_long_datasets <- quiet(formatForWeights(object, data, imputed_datasets))


  save(list=c("wide_long_datasets"),
       file="/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
}


test_that("formatted data contain exposures and outcomes", {

  run_code_f_w()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")

  expect_true(
    sum(as.logical(unlist(
    lapply(wide_long_datasets, function(x) c(apply(expand.grid(object$exposures, as.character(as.numeric(object$exposure_time_pts))), 1, paste0, sep="", collapse="."),
                apply(expand.grid(object$outcomes, object$outcome_time_pt), 1, paste0, sep="", collapse=".")) %in%
                colnames(x)))))==length(c(apply(expand.grid(object$exposures, as.character(as.numeric(object$exposure_time_pts))), 1, paste0, sep="", collapse="."),
                                                             apply(expand.grid(object$outcomes, object$outcome_time_pt), 1, paste0, sep="", collapse=".")))
  )


})
