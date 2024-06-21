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


test_that(".create_hypotheses_mat", {
  mat = .create_hypotheses_mat(
    histories = c("a", "b", "c"), 
    reference = "c", 
    comparison = c("a", "b")
  )

  expect_equal(mat[, 1], c(1, 0, -1))
  expect_equal(mat[, 2], c(0, 1, -1))
  expect_equal(colnames(mat), c("(a) - (c)", "(b) - (c)"))
})
