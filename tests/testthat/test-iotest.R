# test w4mclassfilter::w4m_filter_by_sample_class

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
      my.env$data    <- utils::read.delim(
        fill = FALSE,
        file = file_path,
        stringsAsFactors = TRUE
        )
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

# test matrix outputs and inputs
#' @import testthat w4mclassfilter
#' @export
test_that("io test - matrix", {
  # set up variables
  out_env <- new.env()
  variableMetadata_in  <- "input_variableMetadata.tsv"
  variableMetadata_out <- out_env
  variableMetadata_exp <- "expected_variableMetadata.tsv"
  sampleMetadata_in <- "input_sampleMetadata.tsv"
  sampleMetadata_out <- out_env
  sampleMetadata_exp <- "expected_sampleMetadata.tsv"
  dataMatrix_in <- "input_dataMatrix.tsv"
  dataMatrix_out <- "output_iotest_dataMatrix.tsv"
  dataMatrix_exp <- "expected_filterio_dataMatrix.tsv"
  classes_to_filter <- c("M")
  class_column <- "gender"
  false_to_exclude_classes_in_filter <- TRUE

  # test input ...

  data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  expect_true(data_matrix_input_env$success, info = "read data matrix input")
  #   dataMatrix_in as data.frame:
  dataMatrix_in <- data_matrix_input_env$data
  data_matrix_nrow <- nrow(dataMatrix_in)
  data_matrix_ncol <- ncol(dataMatrix_in) - 1
  rm(data_matrix_input_env)

  sample_metadata_input_env <- read_data_frame(sampleMetadata_in, "sample metadata input")
  expect_true(sample_metadata_input_env$success, info = "read sample metadata input")
  #   sampleMetadata_in as list:
  sampleMetadata_in <- list(sampleMetadata = sample_metadata_input_env$data)
  rm(sample_metadata_input_env)

  variable_metadata_input_env <- read_data_frame(variableMetadata_in, "variable metadata input")
  expect_true(variable_metadata_input_env$success, info = "read variable metadata input")
  #   variableMetadata_in as environment:
  variableMetadata_in <- list2env(list(variableMetadata = variable_metadata_input_env$data))
  rm(variable_metadata_input_env)

  # first of two passes: create an unfiltered matrix as output

  data_matrix_env <- new.env()

  # filter, impute, and write output
  filter_result <- w4m_filter_by_sample_class(
    dataMatrix_in = dataMatrix_in
  , dataMatrix_out = data_matrix_env
  , variableMetadata_in = variableMetadata_in
  , variableMetadata_out = variableMetadata_out
  , sampleMetadata_out = sampleMetadata_out
  , sampleMetadata_in = sampleMetadata_in
  , classes = c()                                 # no classes to exclude
  , include = !false_to_exclude_classes_in_filter # FALSE - exclude only 'no classes to exclude'
  , class_column = class_column
  )
  expect_true(filter_result, info = "filter_result should be true")
  expect_equivalent(ls_env <- ls(data_matrix_env), "dataMatrix", info = sprintf("expected ls(data_matrix_env) == 'dataMatrix' but instead was '%s'", ls_env))
  expect_true(is.matrix(data_matrix_env$dataMatrix), info = "is.matrix(data_matrix_env$dataMatrix) should be true")
  expect_equivalent(nrow(data_matrix_env$dataMatrix), 15, info = "nrow(data_matrix_env$dataMatrix) is not 15 as expected")
  expect_equivalent(ncol(data_matrix_env$dataMatrix), 19, info = "ncol(data_matrix_env$dataMatrix) is not 19 as expected")

  # first of two passes: read an unfiltered matrix producing a filtered data file for dataMatrix

  # filter, impute, and write output
  filter_result <- w4m_filter_by_sample_class(
    dataMatrix_in = data_matrix_env
  , dataMatrix_out = dataMatrix_out
  , variableMetadata_in = variableMetadata_in
  , variableMetadata_out = variableMetadata_out
  , sampleMetadata_out = sampleMetadata_out
  , sampleMetadata_in = sampleMetadata_in
  , classes = classes_to_filter
  , include = false_to_exclude_classes_in_filter
  , class_column = class_column
  )
  expect_true(filter_result, info = "filter_result should be true")

  # read actual outputs and compare to expecteds
  data_matrix_output_env <- read_data_frame(dataMatrix_out, "data matrix output")
  expect_true(data_matrix_output_env$success, info = "read data matrix output")
  data_matrix_expected_env <- read_data_frame(dataMatrix_exp, "data matrix expected")
  expect_true(data_matrix_expected_env$success, info = "read data matrix expected")
  expect_equivalent_length(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")
  expect_equivalent(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")

  sample_metadata_expected_env <- read_data_frame(sampleMetadata_exp, "sample metadata expected")
  expect_true(sample_metadata_expected_env$success, info = "read sample metadata expected")
  expect_true(all.equal(out_env$sampleMetadata, sample_metadata_expected_env$data), info = {
    x <- all.equal(sample_metadata_expected_env$data, out_env$sampleMetadata)
    sprintf("all.equal(sample_metadata_expected_env$data, out_env$sampleMetadata) = %s", x)
  })

  variable_metadata_expected_env <- read_data_frame(variableMetadata_exp, "variable metadata expected")
  expect_true(variable_metadata_expected_env$success, info = "read variable metadata expected")
  expect_true(all.equal(out_env$variableMetadata, variable_metadata_expected_env$data), info = {
    x <- all.equal(variable_metadata_expected_env$data, out_env$variableMetadata)
    sprintf("all.equal(variable_metadata_expected_env$data, out_env$variableMetadata) = %s", x)
  })
})

# test alternative inputs and outputs
#' @import testthat w4mclassfilter
#' @export
test_that("io test", {
  # set up variables
  out_env <- new.env()
  variableMetadata_in  <- "input_variableMetadata.tsv"
  variableMetadata_out <- out_env
  variableMetadata_exp <- "expected_variableMetadata.tsv"
  sampleMetadata_in <- "input_sampleMetadata.tsv"
  sampleMetadata_out <- out_env
  sampleMetadata_exp <- "expected_sampleMetadata.tsv"
  dataMatrix_in <- "input_dataMatrix.tsv"
  dataMatrix_out <- "output_iotest2_dataMatrix.tsv"
  dataMatrix_exp <- "expected_filterio_dataMatrix.tsv"
  classes_to_filter <- c("M")
  class_column <- "gender"
  false_to_exclude_classes_in_filter <- TRUE

  # test input ...

  data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  expect_true(data_matrix_input_env$success, info = "read data matrix input")
  #   dataMatrix_in as data.frame:
  dataMatrix_in <- data_matrix_input_env$data
  rm(data_matrix_input_env)

  sample_metadata_input_env <- read_data_frame(sampleMetadata_in, "sample metadata input")
  expect_true(sample_metadata_input_env$success, info = "read sample metadata input")
  #   sampleMetadata_in as list:
  sampleMetadata_in <- list(sampleMetadata = sample_metadata_input_env$data)
  rm(sample_metadata_input_env)

  variable_metadata_input_env <- read_data_frame(variableMetadata_in, "variable metadata input")
  expect_true(variable_metadata_input_env$success, info = "read variable metadata input")
  #   variableMetadata_in as environment:
  variableMetadata_in <- list2env(list(variableMetadata = variable_metadata_input_env$data))
  rm(variable_metadata_input_env)

  # filter, impute, and write output
  filter_result <- w4m_filter_by_sample_class(
    dataMatrix_in = dataMatrix_in
  , dataMatrix_out = dataMatrix_out
  , variableMetadata_in = variableMetadata_in
  , variableMetadata_out = variableMetadata_out
  , sampleMetadata_out = sampleMetadata_out
  , sampleMetadata_in = sampleMetadata_in
  , classes = classes_to_filter                  # filter out this list of class names
  , include = false_to_exclude_classes_in_filter # TRUE - exclude 'classes_to_filter' in previous line
  , class_column = class_column
  )
  expect_true(filter_result, info = "filter_result should be true")

  # read actual outputs and compare to expecteds
  data_matrix_output_env <- read_data_frame(dataMatrix_out, "data matrix output")
  expect_true(data_matrix_output_env$success, info = "read data matrix output")
  data_matrix_expected_env <- read_data_frame(dataMatrix_exp, "data matrix expected")
  expect_true(data_matrix_expected_env$success, info = "read data matrix expected")
  expect_equivalent_length(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")
  expect_equivalent(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")

  sample_metadata_expected_env <- read_data_frame(sampleMetadata_exp, "sample metadata expected")
  expect_true(sample_metadata_expected_env$success, info = "read sample metadata expected")
  expect_true(all.equal(out_env$sampleMetadata, sample_metadata_expected_env$data), info = {
    x <- all.equal(sample_metadata_expected_env$data, out_env$sampleMetadata)
    sprintf("all.equal(sample_metadata_expected_env$data, out_env$sampleMetadata) = %s", x)
  })

  variable_metadata_expected_env <- read_data_frame(variableMetadata_exp, "variable metadata expected")
  expect_true(variable_metadata_expected_env$success, info = "read variable metadata expected")
  expect_true(all.equal(out_env$variableMetadata, variable_metadata_expected_env$data), info = {
    x <- all.equal(variable_metadata_expected_env$data, out_env$variableMetadata)
    sprintf("all.equal(variable_metadata_expected_env$data, out_env$variableMetadata) = %s", x)
  })
})
