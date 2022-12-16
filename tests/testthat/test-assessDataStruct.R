
#things to test
#class


test_that("data frame is produced", {
  object=msmObject(data_path ="/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/MSM_ALL_LONG_bin_8.25.20.csv",
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
            outcome_time_pt=154)

  expect_s3_class(formatDataStruct(object), "data.frame")

  })


test_that("warning is produced", {
  object=msmObject(data_path ="/Users/isabella/Library/CloudStorage/Box-Box/BSL General/Isa/MSMs/Meriah test/MSM_ALL_LONG_bin_8.25.20.csv",
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
                   outcome_time_pt=154)

  expect_warning(formatDataStruct(object))

})


