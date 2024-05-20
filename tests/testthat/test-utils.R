library(testthat)
test_that(".extract_time_pts_from_vars", {
  x1 = .extract_time_pts_from_vars(
    vars = c("A.1", "A.2", "C"), 
    sep = "\\."
  ) 
  expect_equal(unname(x1), c(1, 2, NA))

  x2 = .extract_time_pts_from_vars(
    vars = c("A_w1", "A_w2", "C"), 
    sep = "_w"
  ) 
  expect_equal(unname(x2), c(1, 2, NA))
})
