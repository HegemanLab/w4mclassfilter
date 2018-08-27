# test w4mclassfilter::w4m_filter_by_sample_class

expect_equivalent_length <- function(target, current, info) {
  expect_equivalent(length(target), length(current), info = info)
}
expect_all_equal <- function(target, current, info) {
  expect_equivalent_length(target = target, current = current, info = info)
  expect_true(all.equal(target = target, current = current, info = info))
}

my_read_df <- function(file_path, kind_string, my_failure_action = print) {
  # ---
  # read in the data frame
  my.env <- new.env()
  my.env$success <- FALSE
  my.env$msg <- sprintf("no message reading %s", kind_string)
  if (!file.exists(file_path)) {
    my.env$msg <- sprintf("file '%s' not found when trying to read %s", file_path, kind_string)
    my_failure_action(my.env$msg)
    return ( NULL )
  }
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
    my_failure_action(paste("my_read_df - ", my.env$msg))
    return ( NULL )
  }
  return (my.env)
}

run_nofilter_test <- function(classes_to_filter, class_column, samplename_column = "sampleMetadata", false_to_exclude_classes_in_filter) {
  # set up variables
  variableMetadata_in  <- "missing_input_variableMetadata.tsv"
  variableMetadata_out <- "output_variableMetadata.tsv"
  variableMetadata_exp <- "expected_nofilter_variableMetadata.tsv"
  sampleMetadata_in <- "input_sampleMetadata.tsv"
  sampleMetadata_out <- "output_sampleMetadata.tsv"
  sampleMetadata_exp <- "expected_nofilter_sampleMetadata.tsv"
  dataMatrix_in <- "input_dataMatrix.tsv"
  dataMatrix_out <- "output_nofilter_dataMatrix.tsv"
  dataMatrix_exp <- "expected_nofilter_dataMatrix.tsv"
  # test input files
  data_matrix_input_env <- my_read_df(dataMatrix_in, "data matrix input")
  expect_true(data_matrix_input_env$success, info = "read data matrix input")
  rm(data_matrix_input_env)
  sample_metadata_input_env <- my_read_df(sampleMetadata_in, "sample metadata input")
  expect_true(sample_metadata_input_env$success, info = "read sample metadata input")
  rm(sample_metadata_input_env)
  variable_metadata_input_env <- my_read_df(
    variableMetadata_in
  , "variable metadata input"
  , my_failure_action = invisible
  )
  expect_true(is.null(variable_metadata_input_env), info = "read variable metadata input")
  rm(variable_metadata_input_env)
  error_list <- c("Error list:")
  add_error <- function(...) {
    error_list <<- c(error_list, paste(...))
  }
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
  , failure_action = add_error
  )
  expect_false(filter_result, info = "filter_result should be false in the case of a missing file")
  if (filter_result) {
    cat(paste("\n", error_list))
  }
}

#' @import testthat w4mclassfilter
#' @export
test_that("missing file test", {
  run_nofilter_test(classes_to_filter = c("M"), class_column = "", samplename_column = "sampleMetadata", false_to_exclude_classes_in_filter = TRUE)
})
