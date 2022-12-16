

run_code <- function(){
  quiet <- function(x) { #quiets cat output
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  object=msmObject(data_path ="/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/MSM_ALL_LONG_bin_8.25.20.csv",
                   home_dir = "/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/", #string
                   ID = "s_id", #string
                   time_pts=c(6, 15, 24, 35, 58, 90, 154), #list of integers
                   time_var="WAVE", #string
                   missing=-9999, #string or integer(s)
                   time_varying_variables= c("AGE", "WIND", "PCX_POS", "PCX_NEG", "BSI", "INR_bin", "INR", "TIMEofDAY", "CORTBASE", "sAABASE", "ESETA1", "HOMEETA1", "CTSETA1"),
                   continuous_variables=c("GrosPay1", "RMomAgeU", "RMAge1st", "SwghtLB","RWghtLb", "IBQDnovm", "MDI", "AGE", "WIND", "PCX_POS", "PCX_NEG", "BSI", "TIMEofDAY", "CORTBASE", "sAABASE", "ESETA1", "HOMEETA1", "CTSETA1", "INR"),
                   factor_covariates=c("NC","TcSex2", "TcBlac2", "INR_bin"), #list of characters
                   exposures=c("INR", "HOMEETA1"), #list of strings
                   exposure_time_pts=c(6, 15), #list of integers
                   exposure_epochs=data.frame(epochs=c("Infancy", "Toddlerhood", "Childhood"), #user-created names of epochs as a list of strings
                                              values=I(list(c(6,15), c(24,35), c(58, 90, 154)))), #time point(s) encompassed in each epoch that correspond to data as a list of numeric lists
                   outcomes=c("CORTBASE"), #list of strings
                   m=1,
                   outcome_time_pt=154)

  data <- quiet(suppressWarnings(formatDataStruct(object)))
  potential_covariates <- quiet(identifyCovariates(object, data))
  covariates_to_include <-quiet(identifyPotentialConfounds(object))
  data_to_impute <- quiet(dataToImpute(object, covariates_to_include))

  # return(covariates_to_include)
  save(list=c("object", "data", "potential_covariates", "covariates_to_include", "data_to_impute"),
       file="/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")
}


test_that("potential_confounders has a list for each exposure-outcome pair", {

  run_code()

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")

  expect_length(covariates_to_include, nrow(expand.grid(object$exposures, object$outcomes)))

  # expect_match(covariates_to_include[[1]], potential_covariates, ignore.case=T)
})

test_that("potential_confounders have correlations above user-specified value", {

  run_code()

  load("/Users/isabella/Documents/Github/MSMs/tests/testthat/fixtures/data_to_impute.Rdata")

  expect_true(
    sum(unlist(lapply(covariates_to_include, function(x) sum(abs(x$cor) >object$balance_thresh)==length(x$cor))))==nrow(expand.grid(object$exposures, object$outcomes)))

})
