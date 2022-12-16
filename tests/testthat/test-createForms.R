
run_code_forms <- function(){
  quiet <- function(x) { #quiets cat output
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")

  forms <- quiet(createForms(object, wide_long_datasets, covariates_to_include))


  save(list=c("forms"),
       file="/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")
}




test_that("forms DO NOT contain excluded variables", {

  run_code_forms()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")

  if(length(object$exclude_covariates)==0){
    skip("No exclude_covariates")}

  expect_false(
    sum(as.logical(unlist(lapply(forms, function(x) object$exclude_covariates %in% gsub(" ", "", c(paste0(unlist(strsplit(paste0(x),"\\+")))))))))>0
  )

})

test_that("forms DO contain mandatory keep variables", {

  run_code_forms()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")

  if(length(object$mandatory_keep_covariates)==0){
    skip("No mandatory_keep_covariates")}

  expect_true(
    sum(as.logical(unlist(lapply(forms, function(x) object$mandatory_keep_covariates %in% gsub(" ", "", c(paste0(unlist(strsplit(paste0(x),"\\+")))))))))==length(forms)
  )

})


test_that("forms DO NOT contain exposure at that exposure time point", {

  run_code_forms()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/imputed_data_sets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/wide_long_datasets.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/forms.Rdata")


  expect_false(

    sum(as.logical(lapply(seq_along(forms), function(i) unlist(strsplit((paste(gsub("form_", "", sapply(strsplit(names(forms)[[i]], "-"), "[", 1)), sapply(strsplit(names(forms)[[i]], "-"), "[", 3), sep=".", collapse=" ")), " ")) %in%
             gsub(" ", "", c(paste0(unlist(strsplit(paste0(forms[[i]]),"\\+"))))))))>0
  )

})


