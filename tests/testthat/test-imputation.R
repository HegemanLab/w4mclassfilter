# test w4mclassfilter::w4m_filter_imputation

library(base)
library(testthat)
library(w4mclassfilter)

test_that("imputation test", {
  my_input <- matrix(c(NA,1,-1,2), ncol = 2, nrow = 2)
  my_expected <- matrix(c(0,1,0,2), ncol = 2, nrow = 2)
  my_output <- w4m_filter_imputation(my_input)
  expect_equivalent(my_output, my_expected, info = "imputation validation")
})
