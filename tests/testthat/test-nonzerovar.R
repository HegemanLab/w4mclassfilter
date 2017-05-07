# test w4mclassfilter::w4m_filter_imputation

library(base)
library(testthat)
library(w4mclassfilter)

test_that("nonzero variance rows and columns test", {
  m <- matrix(
    c(
      1, 2, 3, 4,
      3, 3, 3, 4,
      5, 7, 11, 4,
      13, 17, 19, 4
    )
  , nrow = 4
  , ncol = 4
  , byrow = TRUE
  )
  rownames(m) <- c("A", "B", "C", "D")
  colnames(m) <- c("W", "X", "Y", "Z")

  expected <- matrix(
    c(
      1, 2, 3,
      5, 7, 11,
      13, 17, 19
    )
    , nrow = 3
    , ncol = 3
    , byrow = TRUE
  )
  rownames(expected) <- c("A", "C", "D")
  colnames(expected) <- c("W", "X", "Y")

  expect_equivalent(w4m__nonzero_var(m), expected, info = "nonzero variance rows and columns validation")
})
