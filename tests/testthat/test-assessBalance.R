run_code_balance <- function(){
  quiet <- function(x) { #quiets cat output
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/weights_models.Rdata")

  unbalanced_covariates_for_models <- quiet(suppressWarnings(assessBalance(object, data, weights_models)))


  save(list=c("unbalanced_covariates_for_models"),
       file="/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/unbalanced_covariates_for_models.Rdata")
}


test_that("the only variables missing from covariate table vs unbalanced covars for model are those re: exposure", {
  run_code_balance()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/weights_models.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/unbalanced_covariates_for_models.Rdata")

  for (x in 1:length(unbalanced_covariates_for_models)){

    if (length(unlist(strsplit(unbalanced_covariates_for_models[[x]], "\\+")))>0){
      table=read.csv(paste0(object$home_dir, "balance/comparison values/all_sig_post_balance_cors_", names(unbalanced_covariates_for_models[1]), ".csv"))


      expect_true(sum(grepl(sapply(strsplit(names(unbalanced_covariates_for_models[1]), "-"), "[", 1),
            table$covariate[!gsub("\\.", "_", table$covariate) %in% gsub(" ", "", unlist(strsplit(unbalanced_covariates_for_models[[x]], "\\+")))]))==
        sum(gsub("\\.", "_", table$covariate) %in% gsub(" ", "", unlist(strsplit(unbalanced_covariates_for_models[[x]], "\\+")))==F))
    }
  }
})


test_that("unbalanced covariates have residual correlation/std mean diff > balance threshold ", {
  run_code_balance()
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/weights_models.Rdata")
  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/unbalanced_covariates_for_models.Rdata")

  for (x in 1:length(unbalanced_covariates_for_models)){

    if (length(unlist(strsplit(unbalanced_covariates_for_models[[x]], "\\+")))>0){
  table=read.csv(paste0(object$home_dir, "balance/comparison values/all_sig_post_balance_cors_", names(unbalanced_covariates_for_models[1]), ".csv"))

  expect_true(sum(abs(table$mean_corr[which(gsub("\\.", "_", table$covariate) %in% gsub(" ", "", unlist(strsplit(unbalanced_covariates_for_models[[x]], "\\+"))))])>
        object$balance_thresh)==length(which(gsub("\\.", "_", table$covariate) %in% gsub(" ", "", unlist(strsplit(unbalanced_covariates_for_models[[x]], "\\+"))))))
}}
})





