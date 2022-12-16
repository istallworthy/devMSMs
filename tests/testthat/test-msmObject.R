

#one or more test_that cases, each test has one or more expectations
test_that("msm object is created as a list with all required inputs", {
  # expect_equal(2 * 2, 4)

  #creates a list
  expect_type(msmObject(data_path ="/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/MSM_ALL_LONG_bin_8.25.20.csv",
                        home_dir = "/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/", #string
                        ID = "s_id", #string
                        time_pts=c(6, 15, 24, 35, 58, 90, 154), #list of integers
                        time_var="WAVE", #string
                        missing=-9999, #string or integer(s)
                        time_varying_variables= c("AGE", "WIND", "PCX_POS", "PCX_NEG", "BSI", "INR_bin", "INR", "TIMEofDAY", "CORTBASE", "sAABASE", "ESETA1", "HOMEETA1", "CTSETA1"),
                        factor_covariates=c("NC","TcSex2", "TcBlac2", "INR_bin"), #list of characters
                        exposures=c("INR", "HOMEETA1"), #list of strings
                        exposure_time_pts=c(6, 15), #list of integers
                        exposure_epochs=data.frame(epochs=c("Infancy", "Toddlerhood", "Childhood"), #user-created names of epochs as a list of strings
                                                   values=I(list(c(6,15), c(24,35), c(58, 90, 154)))), #time point(s) encompassed in each epoch that correspond to data as a list of numeric lists
                        outcomes=c("CORTBASE"), #list of strings
                        outcome_time_pt=154),
              "list")})



#throws error in absence of required input; e.g., wrong class, NA
test_that("throws error in absence of a required input",{
  expect_error(msmObject(data_path ="/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/MSM_ALL_LONG_bin_8.25.20.csv",
                         home_dir = "/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/", #string
                         ID = "s_id", #string
                         time_pts=c(6, 15, 24, 35, 58, 90, 154), #list of integers
                         time_var="WAVE", #string
                         missing=-9999, #string or integer(s)
                         time_varying_variables= c("AGE", "WIND", "PCX_POS", "PCX_NEG", "BSI", "INR_bin", "INR", "TIMEofDAY", "CORTBASE", "sAABASE", "ESETA1", "HOMEETA1", "CTSETA1"),
                         factor_covariates=c("NC","TcSex2", "TcBlac2", "INR_bin"), #list of characters
                         exposures=c("INR", "HOMEETA1"), #list of strings
                         exposure_time_pts=c(6, 15), #list of integers
                         exposure_epochs=data.frame(epochs=c("Infancy", "Toddlerhood", "Childhood"), #user-created names of epochs as a list of strings
                                                    values=I(list(c(6,15), c(24,35), c(58, 90, 154)))), #time point(s) encompassed in each epoch that correspond to data as a list of numeric lists
                         outcomes=c("CORTBASE"))
               # outcome_time_pt=154),
  )})



#throws warning
expect_warning()

expect_message()

expect_match(string, "Testing", ignore.case=T)

expect_length(object, n)

expect_setequal(x,y) #tests that every element of x occurs in y, and that every element of y occurs in x. But it won’t fail if x and y happen to have their elements in a different order.

expect_s3_class(model, "lm")

expect_true()

expect_false()



})
