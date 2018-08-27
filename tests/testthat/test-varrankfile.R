# test w4mclassfilter::w4m_filter_imputation

library(base)
library(testthat)
library(w4mclassfilter)

test_that("variance by rank or file test", {
  m <- matrix(
    c(
      1, 2, 3,
      5, 7, 11,
      13, 17, 19
    )
  , nrow = 3
  , ncol = 3
  , byrow = TRUE
  )
  rowvars <- w4m__var_by_rank_or_file(m = m, dim = 1)
  expecteds <- c(var(c(1, 2, 3)), var(c(5, 7, 11)), var(c(13, 17, 19)))
  all.equal(rowvars, expecteds)
  expect_equivalent(rowvars, expecteds, info = "variance by row validation")

  colvars <- w4m__var_by_rank_or_file(m = m, dim = 2)
  expecteds <- c(var(c(1, 5, 13)), var(c(2, 7, 17)), var(c(3, 11, 19)))
  all.equal(colvars, expecteds)
  expect_equivalent(colvars, expecteds, info = "variance by column validation")
})
