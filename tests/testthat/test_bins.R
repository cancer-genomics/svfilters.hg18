context("bins")
test_that("bins", {
  data(bins1kb)
  expect_identical(length(bins1kb),
                   2650192L)
})
