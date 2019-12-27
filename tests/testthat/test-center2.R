# test centering for w4mclassfilter::w4m_filter_by_sample_class


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
  centering = c("none","centroid","median","medoid"),
  order_smpl,
  order_vrbl
  ) {
  #### print("order_smpl")
  #### print(order_smpl)
  # set up variables
  variableMetadata_in <- "in_cent_vm.tsv"
  variableMetadata_out <- switch(centering,
                               "none"     = "out_cent_none_vm.tsv",
                               "centroid" = "out_cent_centroid_vm.tsv",
                               "median"   = "out_cent_median_vm.tsv",
                               "medoid"   = "out_cent_medoid_vm.tsv")
  variableMetadata_exp <- switch(centering,
                               "none"     = "exp_cent_none_vm.tsv",
                               "centroid" = "exp_cent_centroid_vm.tsv",
                               "median"   = "exp_cent_median_vm.tsv",
                               "medoid"   = "exp_cent_medoid_vm.tsv")
  sampleMetadata_in <- "in_cent_sm.tsv"
  sampleMetadata_out <- switch(centering,
                              "none"     = "out_cent_none_sm.tsv",
                              "centroid" = "out_cent_centroid_sm.tsv",
                              "median"   = "out_cent_median_sm.tsv",
                              "medoid"   = "out_cent_medoid_sm.tsv")
  sampleMetadata_exp <- switch(centering,
                              "none"     = "exp_cent_none_sm.tsv",
                              "centroid" = "exp_cent_centroid_sm.tsv",
                              "median"   = "exp_cent_median_sm.tsv",
                              "medoid"   = "exp_cent_medoid_sm.tsv")
  dataMatrix_in <- "in_cent_dm.tsv"
  dataMatrix_out <- switch(centering,
                           "none"        = "out_cent_none_dm.tsv",
                           "centroid"    = "out_cent_centroid_dm.tsv",
                           "median"      = "out_cent_median_dm.tsv",
                           "medoid"      = "out_cent_medoid_dm.tsv")
  dataMatrix_exp <- switch(centering,
                           "none"        = "exp_cent_none_dm.tsv",
                           "centroid"    = "exp_cent_centroid_dm.tsv",
                           "median"      = "exp_cent_median_dm.tsv",
                           "medoid"      = "exp_cent_medoid_dm.tsv")
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
  , order_smpl = order_smpl
  , order_vrbl = order_vrbl
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
  expect_equivalent_length(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix length")
  expect_equivalent_length(sample_metadata_output_env$data, sample_metadata_expected_env$data, info = "validate sample metadata length")
  expect_equivalent_length(variable_metadata_output_env$data, variable_metadata_expected_env$data, info = "validate variable metadata length")
  expect_equivalent(data_matrix_output_env$data, data_matrix_expected_env$data, info = "validate data matrix")
  expect_equivalent(rownames(data_matrix_output_env$data), rownames(data_matrix_expected_env$data), info = "validate data matrix rownames")
  expect_equivalent(colnames(data_matrix_output_env$data), colnames(data_matrix_expected_env$data), info = "validate data matrix colnames")
  expect_equivalent(sample_metadata_output_env$data, sample_metadata_expected_env$data, info = "validate sample metadata")
  expect_equivalent(variable_metadata_output_env$data, variable_metadata_expected_env$data, info = "validate variable metadata")
}

#' @import testthat w4mclassfilter
#' @export
test_that("center2 none test", {
  #print("*** center2 none test ***")
  run_center_test(
    classes_to_filter = c()
  , class_column = "gender"
  , samplename_column = "sampleMetadata"
  , false_to_exclude_classes_in_filter = TRUE
  , centering = "none"
  , order_smpl = c("gender")
  , order_vrbl = c("mz")
  )
})

#' @import testthat w4mclassfilter
#' @export
test_that("center2 centroid test", {
  #print("*** center2 centroid test ***")
  run_center_test(
    classes_to_filter = c()
    , class_column = "gender"
    , samplename_column = "sampleMetadata"
    , false_to_exclude_classes_in_filter = TRUE
    , centering = "centroid"
    , order_smpl = c("sampleMetadata")
    , order_vrbl = c("variableMetadata")
  )
})

#' @import testthat w4mclassfilter
#' @export
test_that("center2 median test", {
  #print("*** center2 median test ***")
  run_center_test(
    classes_to_filter = c()
    , class_column = "gender"
    , samplename_column = "sampleMetadata"
    , false_to_exclude_classes_in_filter = TRUE
    , centering = "median"
    , order_smpl = c("sampleMetadata")
    , order_vrbl = c("mz")
  )
})

#' @import testthat w4mclassfilter
#' @export
test_that("center2 medoid test", {
  #print("*** center2 medoid test ***")
  run_center_test(
    classes_to_filter = c()
    , class_column = "gender"
    , samplename_column = "sampleMetadata"
    , false_to_exclude_classes_in_filter = TRUE
    , centering = "medoid"
    , order_smpl = c("gender")
    , order_vrbl = c("rt")
  )
})
