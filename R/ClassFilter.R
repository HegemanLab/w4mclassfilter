#' @title
#' Impute missing intensities to zero for W4M data matrix (deprecated)
#'
#' @description
#' Substitute zero for missing or negative intensity values in W4M data matrix (synonym for w4m_filter_zero_imputation, deprecated)
#'
#' @param m  matrix: W4M data matrix potentially containing NA or negative values
#'
#' @return matrix: input data matrix with zeros substituted for negative or NA values
#'
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @examples
#' # input contains negative and missing values
#' my_input <- matrix(c(NA,1,-1,2), ncol = 2, nrow = 2)
#'
#' # expected output converts negative and missing values to zero
#' my_expected <- matrix(c(0,1,0,2), ncol = 2, nrow = 2)
#'
#' # run the imputation method to generate actual output
#' my_output <- w4m_filter_imputation(my_input)
#'
#' # validate actual output against expected output
#' all.equal(my_output, my_expected, check.attributes = FALSE)
#'
#' @export
w4m_filter_imputation <-
  zero_imputation <-
  function(m) {
    # replace NA values with zero
    m[is.na(m)] <- 0
    # replace negative values with zero, if applicable (It should never be applicable!)
    m[m < 0] <- 0
    # return matrix as the result
    return (m)
  }

#' @title
#' Impute missing values to zero for W4M data matrix
#'
#' @description
#' Substitute zero for missing or negative intensity values in W4M data matrix
#'
#' @param m  matrix: W4M data matrix potentially containing NA or negative values
#'
#' @return matrix: input data matrix with zeros substituted for negative or NA values
#'
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @examples
#' # input contains negative and missing values
#' my_input <- matrix(c(NA,1,-1,2), ncol = 2, nrow = 2)
#'
#' # expected output converts negative and missing values to zero
#' my_expected <- matrix(c(0,1,0,2), ncol = 2, nrow = 2)
#'
#' # run the imputation method to generate actual output
#' my_output <- w4m_filter_zero_imputation(my_input)
#'
#' # validate actual output against expected output
#' all.equal(my_output, my_expected, check.attributes = FALSE)
#'
#' @export
w4m_filter_zero_imputation <-
  zero_imputation

#' @title
#' Do not impute missing intensities to zero for W4M data matrix, but convert negative intensities to zero
#'
#' @description
#' Substitute zero for negative intensity values in W4M data matrix, but not for missing intensity values
#'
#' @param m  matrix: W4M data matrix potentially containing negative values
#'
#' @return matrix: input data matrix with zeros substituted for negative values
#'
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @examples
#' # input contains negative and missing values
#' my_input <- matrix(c(NA,1,-1,2), ncol = 2, nrow = 2)
#'
#' # expected output converts negative and missing values to zero
#' my_expected <- matrix(c(NA,1,0,2), ncol = 2, nrow = 2)
#'
#' # run the imputation method to generate actual output
#' my_output <- w4m_filter_no_imputation(my_input)
#'
#' # validate actual output against expected output
#' all.equal(my_output, my_expected, check.attributes = FALSE)
#'
#' @export
w4m_filter_no_imputation <-
  function(m) {
    # replace negative values with zero, if applicable
    m[m < 0] <- 0
    return (m)
  }

#' @title
#' Impute missing intensities to median for W4M data matrix
#'
#' @description
#' Substitute median feature intensity (across all samples) for missing values and zero for negative values in W4M data matrix
#'
#' @param m  matrix: W4M data matrix potentially containing NA or negative values
#'
#' @return matrix: input data matrix with zeros substituted for negative values and median substituted for missing values
#'
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @examples
#' # input contains negative and missing values
#' my_input <- matrix(c(NA,-1,3,2), ncol = 2, nrow = 2)
#'
#' # expected output converts negative and missing values to zero
#' my_expected <- matrix(c(3,0,3,2), ncol = 2, nrow = 2)
#'
#' # run the imputation method to generate actual output
#' my_output <- w4m_filter_median_imputation(my_input)
#'
#' # validate actual output against expected output
#' all.equal(my_output, my_expected, check.attributes = FALSE)
#'
#' @export
w4m_filter_median_imputation <-
  function(m) {
    # Substitute NA with median for the row.
    # For W4M datamatrix:
    #   - each row has intensities for one feature
    #   - each column has intensities for one sample
    interpolate_row_median <- function(m) {
      # ref: https://stats.stackexchange.com/a/28578
      #   - Create a data.frame whose columns are features and rows are samples.
      #   - For each feature, substitute NA with the median value for the feature.
      t_result <- sapply(
          as.data.frame(t(m))
        , function(x) {
            x[is.na(x)] <- median(x, na.rm = TRUE)
            x
          }
        , simplify = TRUE
        )
      #   - Recover the rownames discarded by sapply.
      rownames(t_result) <- colnames(m)
      #   - Transform result so that rows are features and columns are samples.
      m <- t(t_result)
      # eliminate negative values
      m[m < 0] <- 0
      return (m)
    }
    return (interpolate_row_median(m))
  }

#' @title
#' Support function to compute variances of matrix rows or columns
#'
#' @description
#' (w4mclassfilter support function) Compute variances of rows or columns of a W4M data matrix
#'
#' @param m  matrix: W4M data matrix for which variance must be computed for rows or columns
#' @param dim  integer: For variances of rows, dim == 1, for variances of columns, dim == 2
#'
#' @return vector of numeric: variances for rows or columns
#'
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://stackoverflow.com/a/25100036}
#'
#' @examples
#'
#' m <- base::matrix(
#'   c(
#'     1, 2, 3,
#'     5, 7, 11,
#'     13, 17, 19
#'   )
#' , nrow = 3
#' , ncol = 3
#' , byrow = TRUE
#' )
#' rowvars <- w4m__var_by_rank_or_file(m = m, dim = 1)
#' expecteds <- c(stats::var(c(1,2,3)),stats::var(c(5,7,11)),stats::var(c(13,17,19)))
#' base::all.equal(rowvars, expecteds)
#' colvars <- w4m__var_by_rank_or_file(m = m, dim = 2)
#' expecteds <- c(stats::var(c(1,5,13)),stats::var(c(2,7,17)),stats::var(c(3,11,19)))
#' base::all.equal(colvars, expecteds)
#'
#' @export
w4m__var_by_rank_or_file <- function(m, dim = 1) {
  if (dim == 1) {
    dim_x_2 <- dim(m)[2]
    if ( dim_x_2 == 0 )
      stop("w4m__var_by_rank_or_file: there are zero columns")
  }
  else if (dim == 2) {
    dim_x_1 <- dim(m)[1]
    if ( dim_x_1 == 0 ) {
      stop("w4m__var_by_rank_or_file: there are zero rows")
    }
    m <- t(m)
  }
  else {
    stop("w4m__var_by_rank_or_file: dim is invalid; for rows, use dim = 1; for colums, use dim = 2")
  }
  return(rowSums( (m - rowMeans(m)) ^ 2) / (dim(m)[2] - 1) )
}

# produce matrix from matrix xpre where all rows and columns having zero variance have been removed
#' @title
#' Support function to eliminate rows or columns that have zero variance
#'
#' @description
#' (w4mclassfilter support function) Produce matrix from matrix xpre where all rows and columns having zero variance have been removed
#'
#' @param m  matrix: W4M data matrix potentially having rows or columns having zero variance
#'
#' @return matrix: input data matrix after removal of rows or columns having zero variance
#'
#' @examples
#'
#' m <- matrix(
#'   c(
#'     1, 2, 3, 4,
#'     3, 3, 3, 4,
#'     5, 7, 11, 4,
#'     13, 17, 19, 4
#'   )
#' , nrow = 4
#' , ncol = 4
#' , byrow = TRUE
#' )
#' rownames(m) <- c("A", "B", "C", "D")
#' colnames(m) <- c("W", "X", "Y", "Z")
#' expected <- matrix(
#'   c(
#'     1, 2, 3,
#'     5, 7, 11,
#'     13, 17, 19
#'   )
#'   , nrow = 3
#'   , ncol = 3
#'   , byrow = TRUE
#' )
#' rownames(expected) <- c("A", "C", "D")
#' colnames(expected) <- c("W", "X", "Y")
#' all.equal(w4m__nonzero_var(m), expected)
#'
#' @export
w4m__nonzero_var <- function(m) {
  nonzero_var_internal <- function(x) {
    if (nrow(x) == 0) {
      utils::str(x)
      stop("matrix has no rows")
    }
    if (ncol(x) == 0) {
      utils::str(x)
      stop("matrix has no columns")
    }
    # exclude any rows or columns with zero variance
    if ( is.numeric(x) ) {
      # exclude any rows with zero variance
      if ( nrow(x) > 0) {
        row_vars <- w4m__var_by_rank_or_file(x, dim = 1)
        nonzero_row_vars <- row_vars > 0
        nonzero_rows <- row_vars[nonzero_row_vars]
        if ( length(rownames(x)) != length(rownames(nonzero_rows)) ) {
          row.names <- attr(nonzero_rows, "names")
          x <- x[ row.names, , drop = FALSE ]
        }
      }
      # exclude any columns with zero variance
      if ( ncol(x) > 0) {
        column_vars <- w4m__var_by_rank_or_file(x, dim = 2)
        nonzero_column_vars <- column_vars > 0
        nonzero_columns <- column_vars[nonzero_column_vars]
        if ( length(colnames(x)) != length(colnames(nonzero_columns)) ) {
          column.names <- attr(nonzero_columns, "names")
          x <- x[ , column.names, drop = FALSE ]
        }
      }
    }
    return (x)
  }

  # purge rows and columns that have zero variance until there are nor more changes
  #   rationale: nonzero_var_internal first purges rows with zero variance, then columns with zero variance,
  #              so there exists the possibility that a row's variance might change to zero when a column is removed;
  #              therefore, iterate till there are no more changes
  if ( is.numeric(m) ) {
    my.nrow <- 0
    my.ncol <- 0
    while ( nrow(m) != my.nrow || ncol(m) != my.ncol ) {
      my.nrow <- nrow(m)
      my.ncol <- ncol(m)
      m <- nonzero_var_internal(m)
    }
  }
  return (m)
}

#' @title
#' Filter W4M data matrix by sample-class
#'
#' @description
#' Filter a set of retention-corrected W4M files (dataMatrix, sampleMetadata, variableMetadata) by sample-class or feature-attributes
#'
#' @details
#' The W4M files dataMatrix, sampleMetadata, and variableMetadata must be a consistent set, i.e.,
#' there must be metadata in the latter two files for all (and only for) the samples and variables named in the columns and rows of dataMatrix.
#'
#' For multivariate statistics functions, samples and variables with zero variance must be eliminated, and missing values are problematic.
#'
#' Furthermore, frequently, it is desirable to analyze a subset of samples (or features) in the dataMatrix.
#'
#' This function manipulates produces a set of files with imputed missing values, omitting features and samples that are not consistently present within the set or have zero variance.
#' Secondly, it provides a selection-capability for samples based on whether their sample names match a regular expression pattern; this capability can be used either to select for samples with matching sample names or to exclude them.
#' Thirdly, it provides a selection-capability for features based on whether their metadata lie within the ranges specified by 'variable_range_filter'.
#'
#' Finally, this function provides as an advanced option to compute one of three types of centers for each treatment:
#' * "centroid" - Return only treatment-centers computed for each treatment as the mean intensity for each feature.
#' * "median" - Return only treatment-centers computed for each treatment as the median intensity for each feature.
#' * "medoid" - Return only treatment-centers computed for each treatement as the sample most similar to the other samples (the medoid).
#'   * By definition, the medoid is the sample having the smallest sum of its distances from other samples in the treatment.
#'   * Distances computed in principal-components space.
#'     * Principal components are uncorrelated, so they are used here to minimize the distortion of computed distances by correlated features.
#' * "none" - Return all samples; do not computing centers
#'
#' Inputs (dataMatrix_in, sampleMetadata_in, variableMetadata_in) may be:
#' * character: path to input tab-separated-values-file (TSV)
#' * data.frame
#' * matrix: allowed for dataMatrix_in only
#' * list: must have a member named "dataMatrix", "sampleMetadata", or "variableMetadata" for dataMatrix_in, sampleMetadata_in, or variableMetadata_in, respectively.
#' * environment: must have a member named "dataMatrix", "sampleMetadata", or "variableMetadata" for dataMatrix_in, sampleMetadata_in, or variableMetadata_in, respectively.
#'
#' Outputs (dataMatrix_out, sampleMetadata_out, variableMetadata_out) may be:
#' * character: path to write a tab-separated-values-file (TSV)
#' * list: will add a member named "dataMatrix", "sampleMetadata", or "variableMetadata" for dataMatrix_out, sampleMetadata_out, or variableMetadata_out, respectively.
#' * environment: will add a member named "dataMatrix", "sampleMetadata", or "variableMetadata" for dataMatrix_out, sampleMetadata_out, or variableMetadata_out, respectively.
#'
#' Please see the package vignette for further details.
#'
#' @param dataMatrix_in          input  data matrix (rows are feature names, columns are sample names
#' @param sampleMetadata_in      input  sample metadata (rows are sample names, one column's name matches class_column)
#' @param variableMetadata_in    input  variable metadata (rows are variable names)
#' @param dataMatrix_out         output data matrix (rows are feature names, columns are sample names
#' @param sampleMetadata_out     output sample metadata (rows are sample names, one column's name matches class_column)
#' @param variableMetadata_out   output variable metadata (rows are variable names)
#' @param classes                character array: names of sample classes to include or exclude; default is an empty array
#' @param include                logical: TRUE, include named sample classes; FALSE (the default), exclude named sample classes
#' @param class_column           character: name of "class" column, defaults to "class"
#' @param samplename_column      character: name of column with sample name, defaults to "sampleMetadata"
#' @param name_varmetadata_col1  logical: TRUE, name column 1 of variable metadata as "variableMetadata"; FALSE, no change; default is TRUE
#' @param name_smplmetadata_col1 logical: TRUE, name column 1 of sample metadata as "sampleMetadata"; FALSE, no change; default is TRUE
#' @param variable_range_filter  character array: array of filters specified as 'variableMetadataColumnName:min:max'; default is empty array
#' @param data_imputation        function(m): default imputation method for 'intb' data, where intensities have background subtracted - impute zero for NA
#' @param order_vrbl             character: name of column of variableMetadata on which to sort, defaults to "variableMetadata" (i.e., the first column)
#' @param order_smpl             character: name of column of sampleMetadata on which to sort, defaults to "sampleMetadata" (i.e., the first column)
#' @param centering              character: center samples by class column (which names treatment).  Possible choices: "none", "centroid", "medoid", or "median"
#' @param failure_action         function(x, ...): action to take upon failure - defaults to 'print(x,...)'
#'
#' @return logical: TRUE only if filtration succeeded
#'
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @importFrom utils read.delim write.table str
#' @importFrom stats median prcomp dist
#'
#' @examples
#' \dontrun{
#'   # set the paths to your input files
#'   dataMatrix_in <- "tests/testthat/input_dataMatrix.tsv"
#'   sampleMetadata_in <- "tests/testthat/input_sampleMetadata.tsv"
#'   variableMetadata_in <- "tests/testthat/input_variableMetadata.tsv"
#'
#'   # set the paths to your (nonexistent) output files
#'   #    in a directory that DOES need to exist
#'   dataMatrix_out <- "tests/testthat/output_dataMatrix.tsv"
#'   sampleMetadata_out <- "tests/testthat/output_sampleMetadata.tsv"
#'   variableMetadata_out <- "tests/testthat/output_variableMetadata.tsv"
#'
#'   # Example: running the filter to exclude only unwanted samples
#'   #   include = FALSE means exclude samples with class blankpos
#'   w4m_filter_by_sample_class(
#'     dataMatrix_in = dataMatrix_in
#'   , dataMatrix_out = dataMatrix_out
#'   , variableMetadata_in = variableMetadata_in
#'   , variableMetadata_out = variableMetadata_out
#'   , sampleMetadata_out = sampleMetadata_out
#'   , sampleMetadata_in = sampleMetadata_in
#'   , classes = c("M")
#'   , include = TRUE
#'   , class_column = "gender"
#'   , samplename_column = "sampleMetadata"
#'   , name_varmetadata_col1 = TRUE
#'   , name_smplmetadata_col1 = TRUE
#'   , variable_range_filter = c()
#'   , data_imputation = w4m_filter_zero_imputation
#'   , order_vrbl = "variableMetadata"
#'   , order_smpl = "sampleMetadata"
#'   , centering  = "none"
#'   , failure_action = function(...) { cat(paste(..., SEP = "\n")) }
#'   )
#' }
#'
#' @export
w4m_filter_by_sample_class <- function(
  dataMatrix_in                           # character:          path to input file containing data matrix (tsv, rows are feature names, columns are sample names)
, sampleMetadata_in                       # character:          path to input file containing sample metadata (tsv, rows are sample names, one column is "class")
, variableMetadata_in                     # character:          path to input file containing variable metadata (tsv, rows are variable names)
, dataMatrix_out                          # character:          path to output file containing data matrix (tsv, rows are feature names, columns are sample names)
, sampleMetadata_out                      # character:          path to output file containing sample metadata (tsv, rows are sample names, one column is "class")
, variableMetadata_out                    # character:          path to output file containing variable metadata (tsv, rows are variable names)
, classes = c()                           # character array:    names of sample classes to include or exclude; default is an empty array
, include = FALSE                         # logical:            TRUE, include named sample classes; FALSE (the default), exclude named sample classes
, class_column = "class"                  # character:          name of "class" column, defaults to "class"
, samplename_column = "sampleMetadata"    # character:          name of column with sample name, defaults to "sampleMetadata"
, name_varmetadata_col1 = TRUE            # logical:            TRUE, name column 1 of variable metadata as "variableMetadata"; FALSE, no change; default is TRUE
, name_smplmetadata_col1 = TRUE           # logical:            TRUE, name column 1 of sample metadata as "sampleMetadata"; FALSE, no change; default is TRUE
, variable_range_filter = c()             # character array:    array of filters specified as 'variableMetadataColumnName:min:max'; default is empty array
, data_imputation = w4m_filter_zero_imputation   # function(m):   default imputation method is for 'intb' data, where intensities have background subtracted - impute zero for NA or negative
, order_vrbl = "variableMetadata"         # character:          order variables by column whose name is supplied here
, order_smpl = "sampleMetadata"           # character:          order samples by column whose name is supplied here
, centering  = c("none", "centroid", "median", "medoid")[1]   # character: center samples by class column (which names treatment)
, failure_action = function(...) { cat(paste(..., SEP = "\n")) }   # function(x, ...):   action to take upon failure - defaults to 'print(x,...)'
) {

  (my_failure_action <- (failure_action))
  # ---
  # define internal functions
  # ...

  # read_data_frame - read a w4m data frame, with error handling
  #   e.g., data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  read_data_frame <- function(file_path, kind_string, failure_action = my_failure_action) {
    # ---
    # read in the data frame
    my.env <- new.env()
    my.env$success <- FALSE
    my.env$msg <- sprintf("no message reading %s", kind_string)
    if (!file.exists(file_path)) {
      my.env$msg <- sprintf("file '%s' not found when trying to read %s", file_path, kind_string)
      failure_action(my.env$msg)
      return ( NULL )
    }
    tryCatch(
      expr = {
        my.env$data    <- utils::read.delim( fill = FALSE, file = file_path, stringsAsFactors = FALSE )
        my.env$success <- TRUE
      }
    , error = function(e) {
       my.env$ msg <- sprintf("%s read failed", kind_string)
      }
    )
    if (!my.env$success) {
      failure_action(my.env$msg)
      return ( NULL )
    }
    return (my.env)
  }

  # return FALSE if any paths are exact duplicates
  #   N.B. This does not check for different relative paths that resolve to the same file.
  my_list <- list(dataMatrix_in, dataMatrix_out, sampleMetadata_in, sampleMetadata_out, variableMetadata_in, variableMetadata_out)
  my.paths <- unlist(my_list[sapply(my_list, is.character)])
  if ( length(my.paths) != length(unique(my.paths)) ) {
    failure_action("some paths are duplicated")
    for ( my.path in my.paths ) {
      failure_action(my.path)
    }
    stop("some paths are duplicated")
    return (FALSE)
  }

  # tryCatchFunc wraps an expression that produces a value if it does not stop:
  #   tryCatchFunc produces a list
  #   On success of expr(), tryCatchFunc produces
  #     list(success TRUE, value = expr(), msg = "")
  #   On failure of expr(), tryCatchFunc produces
  #     list(success = FALSE, value = NA, msg = "the error message")
  tryCatchFunc <- function(expr) {
    # format error for logging
    format_error <- function(e) {
      return(
        paste(
          c("Error { message:"
          , e$message
          , ", call:"
          , e$call
          , "}"
          )
        , collapse = " "
        )
      )
    }
    retval <- NULL
    tryCatch(
      expr = {
        retval <- ( list( success = TRUE, value = eval(expr = expr), msg = "" ) )
      }
      , error = function(e) {
        retval <<- list( success = FALSE, value = NA, msg = format_error(e) )
      }
    )
    return (retval)
  }

  # read one of three XCMS data elements: dataMatrix, sampleMetadata, variableMetadata
  # returns respectively: matrix, data.frame, data.frame, or FALSE if there is a failure
  read_xcms_data_element <- function(xcms_data_in, xcms_data_type, failure_action = stop) {
    my_failure_action <- function(...) {
      failure_action("w4mclassfilter::w4m_filter_by_sample_class::read_xcms_data_element: ", ...)
    }
    # xcms_data_type must be in c("sampleMetadata", "variableMetadata", "dataMatrix")
    if ( ! is.character(xcms_data_type) ) {
      my_failure_action(sprintf("bad parameter xcms_data_type '%s'", deparse(xcms_data_type)))
      return ( NULL )
    }
    if ( 1 != length(xcms_data_type)
         || ! ( xcms_data_type %in% c("sampleMetadata", "variableMetadata", "dataMatrix") )
    ) {
      my_failure_action( sprintf("bad parameter xcms_data_type '%s'", xcms_data_type) )
      return ( NULL )
    }
    if ( is.character(xcms_data_in) ){
      # case: xcms_data_in is a path to a file
      xcms_data_input_env <- read_data_frame(
        xcms_data_in
      , sprintf("%s input", xcms_data_type)
      )
      if (is.null(xcms_data_input_env)) {
        action_msg <- sprintf(
            "read_data_frame failed for '%s', data type '%s'"
          , toString(xcms_data_in)
          , xcms_data_type)
        my_failure_action(action_msg)
        return ( NULL )
      }
      if (!xcms_data_input_env$success) {
        my_failure_action(xcms_data_input_env$msg)
        my_failure_action("w4mclassfilter::w4m_filter_by_sample_class::read_xcms_data_element is returning NULL")
        return ( NULL )
      }
      return (xcms_data_input_env$data)
    } else if ( is.data.frame(xcms_data_in) || is.matrix(xcms_data_in) ) {
      # case: xcms_data_in is a data.frame or matrix
      return(xcms_data_in)
    } else if ( is.list(xcms_data_in) || is.environment(xcms_data_in) ) {
      # NOTE WELL: is.list succeeds for data.frame, so the is.data.frame test must appear before the is.list test
      # case: xcms_data_in is a list
      if ( ! exists(xcms_data_type, where = xcms_data_in) ) {
        my_failure_action(
          sprintf(
            "%s xcms_data_in is missing member '%s'"
          , ifelse(
              is.environment(xcms_data_in)
            , "environment"
            , "list"
            )
          , xcms_data_type
          )
        )
        return (NULL)
      }
      prospect <- getElement(name = xcms_data_type, object = xcms_data_in)
      if ( ! is.data.frame(prospect) && ! is.matrix(prospect) ) {
        utils::str("list - str(prospect)")
        utils::str(prospect)
        if ( is.list(xcms_data_in) ) {
          my_failure_action(sprintf("the first member of xcms_data_in['%s'] is neither a data.frame nor a matrix but is a %s", xcms_data_type, typeof(prospect)))
        } else {
          my_failure_action(sprintf("the first member of xcms_data_in$%s is neither a data.frame nor a matrix but is a %s", xcms_data_type, typeof(prospect)))
        }
        return (prospect)
      }
      return (prospect)
    } else {
      # case: xcms_data_in is invalid
      my_failure_action( sprintf("xcms_data_in has unexpected type %s", typeof(xcms_data_in)) )
      return (NULL)
    }
  }

  # ---
  # entrypoint
  # ...

  # ---
  # read in the sample metadata
  read_data_result <- tryCatchFunc(
    expr = {
      read_xcms_data_element(
        xcms_data_in = sampleMetadata_in
      , xcms_data_type = "sampleMetadata"
      )
    }
  )
  if (is.null(read_data_result)) {
    failure_action("read_xcms_data_element returnd null; aborting w4m_filter_by_sample_class")
  }
  if ( read_data_result$success ) {
    smpl_metadata <- read_data_result$value
  } else {
    failure_action(read_data_result$msg)
    return (FALSE)
  }

  # extract rownames
  if (name_smplmetadata_col1) {
    colnames(smpl_metadata)[1] <- "sampleMetadata"
  }
  rownames(smpl_metadata) <- smpl_metadata[ , samplename_column]

  if (nchar(class_column) > 0 && length(classes) > 0) {
    # select the first column of the rows indicated by classes, include, & class_column, but don't drop dimension
    #   > Reduce(`|`,list(c(TRUE,FALSE,FALSE),c(FALSE,TRUE,FALSE),c(FALSE,FALSE,FALSE)))
    #   [1]  TRUE  TRUE FALSE
    #   > Reduce(`|`,lapply(X=c("[aC]", "[Ab]"), FUN = function(pattern) {grepl(pattern = pattern, x = c("b", "ba", "c", "ab", "bcb")) }))
    #   [1]  TRUE  TRUE FALSE  TRUE  TRUE
    selected_rows <- smpl_metadata[
      xor(
        !include
      , Reduce(
          `|`
        , lapply(X = classes, FUN = function(pattern) {
            grepl(pattern = pattern, x = smpl_metadata[ , class_column])
          })
        )
      )
    , 1
    , drop = FALSE
    ]

    # obtain the row names
    sample_names <- rownames( selected_rows )
  } else {
    sample_names <- rownames( smpl_metadata )
  }
  # ...

  # ---
  # read in the variable metadata
  read_data_result <- tryCatchFunc(
    expr = {
      read_xcms_data_element(xcms_data_in = variableMetadata_in, xcms_data_type = "variableMetadata")
    }
  )
  if (is.null(read_data_result)) {
    failure_action("read_xcms_data_element returnd null; aborting w4m_filter_by_sample_class")
  }
  if ( read_data_result$success ) {
    vrbl_metadata <- read_data_result$value
  } else {
    failure_action(read_data_result$msg)
    return (FALSE)
  }


  # extract rownames (using make.names to handle degenerate feature names)
  err.env <- new.env()
  err.env$success <- FALSE
  err.env$msg <- "no message setting vrbl_metadata rownames"
  tryCatch(
    expr = {
      rownames(vrbl_metadata) <- make.names( vrbl_metadata[ , 1 ], unique = TRUE )
      vrbl_metadata[ , 1 ] <- rownames(vrbl_metadata)
      if (name_varmetadata_col1) {
        colnames(vrbl_metadata)[1] <- "variableMetadata"
      }
      err.env$success     <- TRUE
    }
  , error = function(e) {
     err.env$ msg <- sprintf("failed to set rownames for vrbl_metadata read because '%s'", e$message)
    }
  )
  if (!err.env$success) {
    failure_action(err.env$msg)
    return ( FALSE )
  }
  # ...

  # ---
  # read in the data matrix
  #   For W4M, each row has intensities for one feature and each column has intensities for one sample
  read_data_result <- tryCatchFunc(
    expr = {
      read_xcms_data_element(xcms_data_in = dataMatrix_in, xcms_data_type = "dataMatrix")
    }
  )
  if (is.null(read_data_result)) {
    failure_action("read_xcms_data_element returnd null; aborting w4m_filter_by_sample_class")
  }
  if ( read_data_result$success ) {
    data_matrix <- read_data_result$value
  } else {
    failure_action(read_data_result$msg)
    return (FALSE)
  }

  if ( ! is.matrix(data_matrix) ) {
    # extract rownames (using make.names to handle degenerate feature names)
    err.env <- new.env()
    err.env$success <- FALSE
    err.env$msg <- "no message setting data_matrix rownames"
    tryCatch(
      expr = {
        rownames(data_matrix) <- make.names( data_matrix[ , 1 ], unique = TRUE )
        err.env$success     <- TRUE
      }
    , error = function(e) {
       err.env$msg <- sprintf("failed to set rownames for data_matrix read because '%s'", e$message)
      }
    )
    if (!err.env$success) {
      failure_action(err.env$msg)
      return ( FALSE )
    }

    # remove rownames column
    data_matrix <- data_matrix[ , 2:ncol(data_matrix) ]

    # select the subset of samples indicated by classes, include, & class_column
    data_matrix <- data_matrix[ , intersect(sample_names, colnames(data_matrix)), drop = FALSE ]

    # convert data_matrix to matrix from data.frame
    data_matrix <- as.matrix(data_matrix)
  }
  # ...

  data_matrix_old <- data_matrix
  # Impute missing values with supplied or default method
  #   (necessary for w4m__nonzero_var)
  data_matrix <- zero_imputation(data_matrix)
  # The user-supplied data_imputation function may include transformation, so it's necessary
  #   to apply it here.  The MAXFEAT filter won't work correctly without this step.
  #   (This addresses https://github.com/HegemanLab/w4mclassfilter/issues/5)
  data_matrix <- data_imputation(data_matrix)

  # ---
  # purge unwanted data from data_matrix
  some_data_were_eliminated <- TRUE
  while (some_data_were_eliminated) {
    # count rows and columns before elimination
    nrow_before <- nrow(data_matrix)
    ncol_before <- ncol(data_matrix)

    # run filters for variable metadata and maximum intensity for each feature
    if (length(variable_range_filter) > 0) {
      # filter variables having out-of-range metadata or intensity maximum
      for (variable_range_filter_string in variable_range_filter) {
        variable_range_filter_string <- sub(":$", ":NA",  variable_range_filter_string)
        variable_range_filter_string <- sub("::", ":NA:", variable_range_filter_string)
        split_list <- strsplit(x = variable_range_filter_string, split = ":", fixed = TRUE)
        if ( length(split_list) == 1 ) {
          split_strings <- split_list[[1]]
          if ( length(split_strings) == 3 ) {
            filter_col <- split_strings[1]
            filter_min <- tryCatch({ as.numeric(split_strings[2]) }, warning = function(w){ -Inf })
            filter_max <- tryCatch({ as.numeric(split_strings[3]) }, warning = function(w){ Inf })
            vrbl_colnames <- colnames(vrbl_metadata)
            if ( filter_col %in% vrbl_colnames ) {
              row_value <- vrbl_metadata[filter_col]
              if (filter_min <= filter_max) {
                # filter specifies an inclusion range
                keep_row <- row_value >= filter_min & row_value <= filter_max
              } else {
                # filter specifies an exclusion range
                keep_row <- row_value > filter_min | row_value < filter_max
              }
              vrbl_metadata <- vrbl_metadata[ keep_row , ]
            } else if (filter_col == "FEATMAX") {
              # apply the function 'max' to rows (1, columns would be 2) of data_matrix
              row_maxima <- apply(data_matrix, 1, max, na.rm = TRUE)
              if (filter_min <= filter_max) {
                # filter specifies an inclusion range
                keep_row <- row_maxima >= filter_min & row_maxima <= filter_max
                data_matrix <- data_matrix[ keep_row , ]
              } else {
                warning("w4m_filter_by_sample_class: FEATMAX filter specified but not applied")
              }
            }
          } else {
            warning("w4m_filter_by_sample_class: split_list is not of the expected length")
          }
        }
      }
    }

    # purge data_matrix of rows and columns that have zero variance
    data_matrix <- w4m__nonzero_var(data_matrix)

    # count rows and columns after elimination
    nrow_after <- nrow(data_matrix)
    ncol_after <- ncol(data_matrix)
    some_data_were_eliminated <- nrow_before != nrow_after | ncol_before != ncol_after
  }

  # purge smpl_metadata and vrbl_metadata of irrelevant rows
  # column names
  sample_names <- intersect(sample_names, colnames(data_matrix))
  sample_order <- order(smpl_metadata[sample_names, order_smpl])
  sample_names <- sample_names[sample_order]
  smpl_metadata <- smpl_metadata[sample_names, ]
  # row names
  variable_names <- intersect( rownames(vrbl_metadata), rownames(data_matrix) )
  variable_order <- order(vrbl_metadata[variable_names, order_vrbl])
  variable_names <- variable_names[variable_order]
  vrbl_metadata <- vrbl_metadata[variable_names, ]

  # Impute missing values with supplied or default method and the ORIGINAL dataMatrix
  #   This is to avoid biasing median-imputation toward the center of the selected features and samples.
  data_matrix <- data_imputation(data_matrix_old)
  # filter out undesired features and samples
  data_matrix <- data_matrix[variable_names, sample_names, drop = FALSE ]
  # ...

  # ---
  if (centering == "centroid") {
    treatments <- smpl_metadata[class_column][[1]]
    nrow_dm <- nrow(data_matrix)
    unitrts <- unique(treatments)
    ntrts <- length(unitrts)
    smpl_metadata <- data.frame(
      trt = unitrts,
      n = sapply(X = unitrts, FUN = function(x) sum(x == treatments)),
      stringsAsFactors = FALSE
    )
    sample_names <- unitrts[order(unitrts)]
    # for each treatment, calculate the mean intensity for each feature
    new_df <- as.data.frame(
      sapply(
        X = unitrts,
        FUN = function(x) {
          unitrt <- x
          sapply(
            X = 1:nrow_dm,
            FUN = function(x) {
              mean(data_matrix[x, unitrt == treatments])
            }
          )
        }
      ),
      stringsAsFactors = FALSE
    )
    rownames(new_df) <- rownames(data_matrix)
    data_matrix <- as.matrix(new_df)
  }
  else if (centering == "median") {
    treatments <- smpl_metadata[class_column][[1]]
    nrow_dm <- nrow(data_matrix)
    unitrts <- unique(treatments)
    ntrts <- length(unitrts)
    smpl_metadata <- data.frame(
      trt = unitrts,
      n = sapply(X = unitrts, FUN = function(x) sum(x == treatments)),
      stringsAsFactors = FALSE
    )
    sample_names <- unitrts[order(unitrts)]
    # for each treatment, calculate the median intensity for each feature
    new_df <- as.data.frame(
      sapply(
        X = unitrts,
        FUN = function(x) {
          unitrt <- x
          sapply(
            X = 1:nrow_dm,
            FUN = function(x) {
              median(data_matrix[x, unitrt == treatments])
            }
          )
        }
      ),
      stringsAsFactors = FALSE
    )
    rownames(new_df) <- rownames(data_matrix)
    data_matrix <- as.matrix(new_df)
  }
  else if (centering == "medoid") {
    # compute medoid (ref: https://www.biostars.org/p/11987/#11989)
    #medoid_col <- function(trt) names(which.min(rowSums(as.matrix(dist(t(trt))))))
    medoid_row  <- function(trt) names(which.min(rowSums(as.matrix(dist(  trt )))))

    treatments <- smpl_metadata[,class_column]
    unitrts <- unique(treatments)
    # When computing principal components with prcomp, set scale. to TRUE because, according to
    #   https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html,
    #   "in general scaling is advisable".
    my_pca <- prcomp(t(data_matrix), scale. = TRUE, tol = sqrt(.Machine$double.eps))
    # Extract eigenvalues to determine how many are < 1
    # ref for extraction: https://stat.ethz.ch/pipermail/r-help/2005-August/076610.html
    ev <- my_pca$sdev^2
    # The cut-off for the scree is somewhat arbitrary,
    #   https://en.wikipedia.org/wiki/Scree_plot, which cites
    #   Norman and Steiner, Biostatistics: The Bare Essentials, p. 201
    #   (https://books.google.com/books?id=8rkqWafdpuoC&pg=PA201)
    # To be conservative, limit the number of PCs to twice the number of eigenvalues that are greater than 1.
    #   It might be better instead to keep adding components until the residual approaches some threshold.
    my_rank <- min(length(ev),2*sum(ev > 1))
    my_scores <- my_pca$x
    my_scores <- my_scores[,1:min(ncol(my_scores),my_rank)]
    # For each treatment, calculate the medoid, i.e.,
    #   the sample with the minimum distance to the other samples in the trt
    my_sapply_result <- sapply(
      X = unitrts,
      FUN = function(x) {
        unitrt <- x
        my_trt_scores <- my_scores[treatments == unitrt,]
        my_trt_medoid <- medoid_row(my_trt_scores)
        return (my_trt_medoid)
      }
    )
    data_matrix <- data_matrix[ , my_sapply_result ]
    # rewrite smpl_metadata:
    #   - rename column 1 as "medoid"
    #   - copy class_column into column 1 as "sampleMetadata"
    smpl_metadata <- smpl_metadata[my_sapply_result,]
    colnames(smpl_metadata)[1] <- "medoid"
    smpl_metadata_colnames <- colnames(smpl_metadata)
    smpl_metadata$sampleMetadata <- smpl_metadata[,class_column]
    smpl_metadata <- smpl_metadata[c("sampleMetadata",smpl_metadata_colnames)]
    rownames(smpl_metadata) <- smpl_metadata$sampleMetadata
    # rename data_matrix columns as class
    colnames(data_matrix) <- smpl_metadata[,class_column]
    # reset sample_names
    sample_names <- smpl_metadata$sampleMetadata
    sample_order <- order(smpl_metadata[sample_names, order_smpl])
    sample_names <- sample_names[sample_order]
  }
  # ...

  # ---
  # write out the results
  err.env <- new.env()
  err.env$success <- FALSE
  err.env$msg <- "no message writing output files"
  err.env$trace <- "trace string:"
  tryCatch(
    expr = {
      err.env$trace <- paste(err.env$trace, "A")
      sub_matrix <- data_matrix[ rownames(data_matrix) %in% variable_names    # row selector
                    , colnames(data_matrix) %in% sample_names      # column selector
                    , drop = FALSE                                 # keep two dimensions
                    ]
      err.env$trace <- paste(err.env$trace, "B")
      # sort matrix to match order of variable_names and sample_names
      sorted_matrix <- sub_matrix[variable_names, sample_names]
      err.env$trace <- paste(err.env$trace, "C")
      # write the data matrix
      if ( is.character(dataMatrix_out) ){

        # write the results
        utils::write.table(  x = sorted_matrix
                           , file = dataMatrix_out
                           , sep = "\t"
                           , quote = FALSE
                           , row.names = TRUE
                           , col.names = NA
                           )
      } else if ( is.environment(dataMatrix_out) || (is.list(dataMatrix_out) && ! is.matrix(dataMatrix_out)) ) {
        dataMatrix_out$dataMatrix <- sorted_matrix
      } else {
        stop(sprintf("w4m_filter_by_sample_class: dataMatrix_out has unexpected type %s"), typeof(dataMatrix_out))
        return (FALSE)
      }

      # write the sample metadata
      if ( is.character(sampleMetadata_out) ){
        utils::write.table( x = smpl_metadata
                             [ sample_names # row selector
                             ,              # column selector (select all)
                             , drop = FALSE # keep two dimensions
                             ]
                           , file = sampleMetadata_out
                           , sep = "\t"
                           , quote = FALSE
                           , row.names = FALSE
                           )
      } else if ( is.environment(sampleMetadata_out) || (is.list(sampleMetadata_out) && ! is.matrix(sampleMetadata_out)) ) {
        sampleMetadata_out$sampleMetadata <-
          smpl_metadata [ sample_names # row selector
                        ,              # column selector (select all)
                        , drop = FALSE # keep two dimensions
                        ]
        rownames(sampleMetadata_out$sampleMetadata) <- 1:nrow(sampleMetadata_out$sampleMetadata)
        sampleMetadata_out$sampleMetadata$sampleMetadata <- as.factor(sampleMetadata_out$sampleMetadata$sampleMetadata)
        sampleMetadata_out$sampleMetadata <- droplevels(sampleMetadata_out$sampleMetadata)
      } else {
        stop(sprintf("w4m_filter_by_sample_class: sampleMetadata_out has unexpected type %s"), typeof(sampleMetadata_out))
        return (FALSE)
      }

      # write the variable metadata
      if ( is.character(variableMetadata_out) ){
        utils::write.table( x = vrbl_metadata
                             [ variable_names # row selector
                             ,                # column selector (select all)
                             , drop = FALSE   # keep two dimensions
                             ]
                           , file = variableMetadata_out
                           , sep = "\t"
                           , quote = FALSE
                           , row.names = FALSE
                           )
      } else if ( is.environment(variableMetadata_out) || (is.list(variableMetadata_out) && ! is.matrix(variableMetadata_out)) ) {
        variableMetadata_out$variableMetadata <-
          vrbl_metadata[ variable_names # row selector
                       ,                # column selector (select all)
                       , drop = FALSE   # keep two dimensions
                       ]
        rownames(variableMetadata_out$variableMetadata) <- 1:nrow(variableMetadata_out$variableMetadata)
        variableMetadata_out$variableMetadata$variableMetadata <- as.factor(variableMetadata_out$variableMetadata$variableMetadata)
        variableMetadata_out$variableMetadata <- droplevels(variableMetadata_out$variableMetadata)
      } else {
        stop(sprintf("w4m_filter_by_sample_class: variableMetadata_out has unexpected type %s"), typeof(variableMetadata_out))
        return (FALSE)
      }

      err.env$success     <- TRUE
    }
  , error = function(e) {
     err.env$ msg <- sprintf("failed to write output files because '%s'; %s", e$message, err.env$trace)
    }
  )
  # ...

  # ---
  # report results
  if (!err.env$success) {
    failure_action(err.env$msg)
    return ( FALSE )
  } else {
    return (TRUE)
  }
  # ...
}

# vim: sw=2 ts=2 et ai :
