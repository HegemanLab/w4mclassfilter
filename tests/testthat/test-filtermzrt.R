# test w4mclassfilter::w4m_filter_by_sample_class

#require(base)
#require(testthat)
#require(w4mclassfilter)

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


#' @import testthat w4mclassfilter
#' @export
test_that("filter mz rt",{
  #expect_true(FALSE, info = "abort filter mz rt")
  # set up variables
  variableMetadata_in  <- "input_variableMetadata.tsv"
  variableMetadata_out <- "output_mzrt_variableMetadata.tsv"
  variableMetadata_exp <- "expected_mzrt_variableMetadata.tsv"
  sampleMetadata_in <- "input_sampleMetadata.tsv"
  sampleMetadata_out <- "output_mzrt_sampleMetadata.tsv"
  sampleMetadata_exp <- "expected_mzrt_sampleMetadata.tsv"
  dataMatrix_in <- "input_dataMatrix.tsv"
  dataMatrix_out <- "output_mzrt_dataMatrix.tsv"
  dataMatrix_exp <- "expected_mzrt_dataMatrix.tsv"
  classes_to_filter <- c("M")
  class_column <- "gender"
  false_to_exclude_classes_in_filter <- TRUE
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
    , classes = classes_to_filter
    , include = false_to_exclude_classes_in_filter 
    , class_column = class_column
    , variable_range_filter = c("mz:125:850","rt:250:850")
  )
  expect_true(filter_result, info = "filter_result should be true")
  #expect_true(FALSE, info = "first checkpoint - filter mz rt")
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
})

