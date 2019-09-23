# test w4mclassfilter::w4m_filter_imputation

library(base)
library(testthat)
library(w4mclassfilter)

test_that("zero-imputation test", {
  my_input <- matrix(c(NA, 1, -1, 2), ncol = 2, nrow = 2)
  my_expected <- matrix(c(0, 1, 0, 2), ncol = 2, nrow = 2)
  my_output <- w4m_filter_zero_imputation(my_input)
  expect_equivalent(my_output, my_expected, info = "zero-imputation validation")
})

test_that("no imputation test", {
  my_input <- matrix(c(NA, 1, -1, 2), ncol = 2, nrow = 2)
  my_expected <- matrix(c(NA, 1, 0, 2), ncol = 2, nrow = 2)
  my_output <- w4m_filter_no_imputation(my_input)
  expect_equivalent(my_output, my_expected, info = "no-imputation validation")
})

test_that("center-imputation test", {
  interpolate_row_median <- function(m) {
    t_result <- sapply(as.data.frame(t(m)), function(x){x[is.na(x)] <- median(x,na.rm=TRUE); x})
    rownames(t_result) <- colnames(m)
    return (t(t_result))
  }
  my_input <- matrix(c(
    NA, 1, 2, 3,
    5, 7, 11, 13,
    17, 19, 23, 29,
    31, 37, 41, 43), ncol = 4, nrow = 4)
  my_expected <- matrix(c(
    17, 1, 2, 3,
    5, 7, 11, 13,
    17, 19, 23, 29,
    31, 37, 41, 43), ncol = 4, nrow = 4)
  my_output <- w4m_filter_median_imputation(my_input)
  #my_output <- interpolate_row_median(my_input)
  expect_equivalent(my_output, my_expected, info = "center-imputation validation")
})
