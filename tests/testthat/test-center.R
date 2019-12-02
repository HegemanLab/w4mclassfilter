# test w4mclassfilter::w4m_filter_by_sample_class

# # how the test data were generated
# sm <- data.frame(
#   sample=c("a1","a2","a3","b1","b2","b3","b4","b5","c1","c2","c3","c4"),
#   trt=c("a","a","a","b","b","b","b","b","c","c","c","c")
#   )
# 
# centers <- c(rep.int(4,3), rep.int(6,5), rep.int(7,4))
# dm.df <- data.frame(
#   v1 = centers + rnorm(12),
#   v2 = centers + rnorm(12),
#   v3 = centers + rnorm(12),
#   v4 = centers + rnorm(12),
#   v5 = centers + rnorm(12)
#   )
# 
# dm <- t(as.matrix(dm.df))
# 
# colnames(dm)<-sm$sample
# 
# vm <- data.frame(
#   variableMetadata = rownames(dm),
#   rt = 100 * rnorm(n=5,mean=5),
#   mz = 100 * rnorm(n=5,mean=5)
#   )
# 
# write.table(vm,"center_vm.tsv",sep="\t",quote=FALSE, row.names = FALSE)
# write.table(sm,"center_sm.tsv",sep="\t",quote=FALSE, row.names = FALSE)
# write.table(dm,"center_dm.tsv",sep="\t",quote=FALSE, col.names = NA)

expect_equivalent_length <- function(target, current, info) {
  expect_equivalent(length(target), length(current), info = info)
}
expect_all_equal <- function(target, current, info) {
  expect_equivalent_length(target = target, current = current, info = info)
  expect_true(all.equal(target = target, current = current, info = info))
}

read_data_frame <- function(file_path, kind_string, failure_action = print) {
  # ---
  # read in the data frame
  my.env <- new.env()
  my.env$success <- FALSE
  my.env$msg <- sprintf("no message reading %s", kind_string)
  tryCatch(
    expr = {
      my.env$data    <- utils::read.delim( fill = FALSE, file = file_path )
      my.env$success <- TRUE
    }
  , error = function(e) {
     my.env$ msg <- sprintf("%s read failed", kind_string)
    }
  )
  if (!my.env$success) {
    failure_action(my.env$msg)
    return ( FALSE )
  }
  return (my.env)
}

run_center_test <- function(
  classes_to_filter,
  class_column,
  samplename_column = "sampleMetadata",
  false_to_exclude_classes_in_filter, 
  centering = c("none","median","representative")
  ) {
  # set up variables
  variableMetadata_in  <- "input_center_vm.tsv"
  variableMetadata_out <- "output_center_vm.tsv"
  variableMetadata_exp <- "expected_center_vm.tsv"
  sampleMetadata_in <- "input_center_sm.tsv"
  sampleMetadata_out <- switch(centering,
                              "none" = "output_center_sm.tsv",
                              "median" = "output_center_median_sm.tsv",
                              "representative" = "output_center_representative_sm.tsv")
  sampleMetadata_exp <- switch(centering,
                              "none" = "expected_center_sm.tsv",
                              "median" = "expected_center_median_sm.tsv",
                              "representative" = "expected_center_representative_sm.tsv")
  dataMatrix_in <- "input_center_dm.tsv"
  dataMatrix_out <- switch(centering,
                           "none" = "output_center_dm.tsv",
                           "median" = "output_center_median_dm.tsv",
                           "representative" = "output_center_representative_dm.tsv")
  dataMatrix_exp <- switch(centering,
                           "none" = "expected_center_dm.tsv",
                           "median" = "expected_center_median_dm.tsv",
                           "representative" = "expected_center_representative_dm.tsv")
  # test input files
  data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  expect_true(data_matrix_input_env$success, info = "read data matrix input")
  rm(data_matrix_input_env)
  sample_metadata_input_env <- read_data_frame(sampleMetadata_in, "sample metadata input")
  expect_true(sample_metadata_input_env$success, info = "read sample metadata input")
  rm(sample_metadata_input_env)
  variable_metadata_input_env <- read_data_frame(variableMetadata_in, "variable metadata input")
  expect_true(variable_metadata_input_env$success, info = "read variable metadata input")
  rm(variable_metadata_input_env)
  # filter, impute, and write output
  filter_result <- w4m_filter_by_sample_class(
    dataMatrix_in = dataMatrix_in
  , dataMatrix_out = dataMatrix_out
  , variableMetadata_in = variableMetadata_in
  , variableMetadata_out = variableMetadata_out
  , sampleMetadata_out = sampleMetadata_out
  , sampleMetadata_in = sampleMetadata_in
  , samplename_column = samplename_column
  , classes = classes_to_filter
  , include = false_to_exclude_classes_in_filter
  , class_column = class_column
  , centering = centering
  )
  expect_true(filter_result, info = "filter_result should be true")
  # read actual output files
  data_matrix_output_env <- read_data_frame(dataMatrix_out, "data matrix output")
  expect_true(data_matrix_output_env$success, info = "read data matrix output")
  sample_metadata_output_env <- read_data_frame(sampleMetadata_out, "sample metadata output")
  expect_true(sample_metadata_output_env$success, info = "read sample metadata output")
  variable_metadata_output_env <- read_data_frame(variableMetadata_out, "variable metadata output")
  expect_true(variable_metadata_output_env$success, info = "read variable metadata output")
  # read expected output files
  data_matrix_expected_env <- read_data_frame(dataMatrix_exp, "data matrix expected")
  expect_true(data_matrix_expected_env$success, info = "read data matrix expected")
  sample_metadata_expected_env <- read_data_frame(sampleMetadata_exp, "sample metadata expected")
  expect_true(sample_metadata_expected_env$success, info = "read sample metadata expected")
  variable_metadata_expected_env <- read_data_frame(variableMetadata_exp, "variable metadata expected")
  expect_true(variable_metadata_expected_env$success, info = "read variable metadata expected")
  # compare actuals with expecteds
  expect_equivalent_length(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")
  expect_equivalent_length(sample_metadata_output_env$data, sample_metadata_expected_env$data, info = "validate sample metadata")
  expect_equivalent_length(variable_metadata_output_env$data, variable_metadata_expected_env$data, info = "validate variable metadata")
  expect_equivalent(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")
  expect_equivalent(sample_metadata_output_env$data, sample_metadata_expected_env$data, info = "validate sample metadata")
  expect_equivalent(variable_metadata_output_env$data, variable_metadata_expected_env$data, info = "validate variable metadata")
}

#' @import testthat w4mclassfilter
#' @export
test_that("center none test", {
  run_center_test(
    classes_to_filter = c()
  , class_column = "trt"
  , samplename_column = "sampleMetadata"
  , false_to_exclude_classes_in_filter = TRUE
  , centering = "none"
  )
})

#' @import testthat w4mclassfilter
#' @export
test_that("center median test", {
  run_center_test(
    classes_to_filter = c()
    , class_column = "trt"
    , samplename_column = "sampleMetadata"
    , false_to_exclude_classes_in_filter = TRUE
    , centering = "median"
  )
})

#' @import testthat w4mclassfilter
#' @export
test_that("center representative test", {
  run_center_test(
    classes_to_filter = c()
    , class_column = "trt"
    , samplename_column = "sampleMetadata"
    , false_to_exclude_classes_in_filter = TRUE
    , centering = "representative"
  )
})
