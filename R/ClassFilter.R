#require(utils)

#' Impute W4M Data Matrix Values
#'
#' Impute missing or negative values in W4M files as zero
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
#' my_input <- matrix(c(NA,1,-1,2), ncol = 2, nrow = 2)
#' my_expected <- matrix(c(0,1,0,2), ncol = 2, nrow = 2)
#' my_output <- w4m_filter_imputation(my_input)
#' all.equal(my_output, my_expected, check.attributes = FALSE)
#'
#' @export
w4m_filter_imputation <- 
  function(m) {
    # replace NA values with zero
    m[is.na(m)] <- 0
    # replace negative values with zero, if applicable (It should never be applicable!)
    if (min(m) < 0) {
      m <- matrix(sapply(X = m, FUN = function(z) {max(z,0)}), nrow = nrow(m) )
    }
    # return matrix as the result
    return (m)
  }

#' Filter W4M Samples by Class of Sample
#'
#' Filter set of retention-corrected W4M files (dataMatrix, sampleMetadata) by feature
#'
#' @param dataMatrix_in        character: path to input file containing data matrix (tsv, rows are feature names, columns are sample names
#' @param sampleMetadata_in    character: path to input file containing sample metadata (tsv, rows are sample names, one column's name matches class_column)
#' @param variableMetadata_in  character: path to input file containing variable metadata (tsv, rows are variable names)
#' @param dataMatrix_out       character: path to output file containing data matrix (tsv, rows are feature names, columns are sample names
#' @param sampleMetadata_out   character: path to output file containing sample metadata (tsv, rows are sample names, one column's name matches class_column)
#' @param variableMetadata_out character: path to input file containing variable metadata (tsv, rows are variable names)
#' @param classes              character array: names of sample classes to include or exclude; default is an empty array
#' @param include              logical: TRUE, include named sample classes; FALSE (the default), exclude named sample classes
#' @param class_column         character: name of "class" column, defaults to "class"
#' @param samplename_column    character: name of column with sample name, defaults to "sampleMetadata"
#' @param data_imputation      function(m): default imputation method for 'intb' data, where intensities have background subtracted - impute zero for NA
#' @param failure_action       function(x, ...): action to take upon failure - defaults to 'print(x,...)'
#'
#' @return logical: TRUE only if filtration succeeded
#' 
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @keywords multivariate
#' @seealso \url{https://github.com/HegemanLab/w4mclassfilter}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @importFrom utils read.delim write.table str
#'
#' @examples
#' print(' # set the paths to your input files                                  ')
#' print(' dataMatrix_in <- "tests/testthat/input_dataMatrix.tsv"               ')
#' print(' sampleMetadata_in <- "tests/testthat/input_sampleMetadata.tsv"       ')
#' print(' variableMetadata_in <- "tests/testthat/input_variableMetadata.tsv"   ')
#' print('                                                                      ')
#' print(' # set the paths to your (nonexistent) output files                   ')
#' print(' #    in a directory that DOES need to exist                          ')
#' print(' dataMatrix_out <- "tests/testthat/output_dataMatrix.tsv"             ')
#' print(' sampleMetadata_out <- "tests/testthat/output_sampleMetadata.tsv"     ')
#' print(' variableMetadata_out <- "tests/testthat/output_variableMetadata.tsv" ')
#' print('                                                                      ')
#' print(' # Example: running the filter to exclude only unwanted samples       ')
#' print(' #   include = FALSE means exclude samples with class blankpos        ')
#' print(' w4m_filter_by_sample_class(                                          ')
#' print('   dataMatrix_in = dataMatrix_in                                      ')
#' print(' , dataMatrix_out = dataMatrix_out                                    ')
#' print(' , variableMetadata_in = variableMetadata_in                          ')
#' print(' , variableMetadata_out = variableMetadata_out                        ')
#' print(' , sampleMetadata_out = sampleMetadata_out                            ')
#' print(' , sampleMetadata_in = sampleMetadata_in                              ')
#' print(' , classes = c("blankpos")                                            ')
#' print(' , include = FALSE                                                    ')
#'
#' @export
w4m_filter_by_sample_class <- function(
  dataMatrix_in                           # character: path to input file containing data matrix (tsv, rows are feature names, columns are sample names)
, sampleMetadata_in                       # character: path to input file containing sample metadata (tsv, rows are sample names, one column is "class")
, variableMetadata_in                     # character: path to input file containing variable metadata (tsv, rows are variable names)
, dataMatrix_out                          # character: path to output file containing data matrix (tsv, rows are feature names, columns are sample names)
, sampleMetadata_out                      # character: path to output file containing sample metadata (tsv, rows are sample names, one column is "class")
, variableMetadata_out                    # character: path to output file containing variable metadata (tsv, rows are variable names)
, classes = c()                           # array of character: names of sample classes to include or exclude; default is an empty array
, include = FALSE                         # logical: TRUE, include named sample classes; FALSE (the default), exclude named sample classes
, class_column = "class"                  # character: name of "class" column, defaults to "class"
, samplename_column = "sampleMetadata"    # character: name of column with sample name, defaults to "sampleMetadata"
, data_imputation = w4m_filter_imputation # function(m): default imputation method is for 'intb' data, where intensities have background subtracted - impute zero for NA
, failure_action = print                  # function(x, ...): action to take upon failure - defaults to 'print(x,...)'
) {
  # ---
  # define internal functions

  # read_data_frame - read a w4m data frame, with error handling
  #   e.g., data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  read_data_frame <- function(file_path, kind_string, failure_action = failure_action) {
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
  
  # MatVar - Compute variance of rows or columns of a matrix
  # ref: http://stackoverflow.com/a/25100036
  # For row variance, dim == 1, for col variance, dim == 2
  MatVar <- function(x, dim = 1) {
    if (dim == 1) {
      dim.x.2 <- dim(x)[2]
      if ( dim.x.2 == 0 )
        stop("MatVar: there are zero columns")
      if ( dim.x.2 == 1 ) {
        stop("MatVar: a single column is insufficient to calculate a variance")
        # return ( rep.int(x = 0, times = nrow(x)) )
      } else {
        return ( rowSums( (x    - rowMeans(x))^2 ) / ( dim(x)[2] - 1 ) )
      }
    } else if (dim == 2) {
      dim.x.1 <- dim(x)[1]
      if ( dim.x.1 == 0 ) {
        stop("MatVar: there are zero rows")
      }
      if ( dim.x.1 == 1 ) {
        stop("MatVar: a single row is insufficient to calculate a variance")
        # return ( rep.int(x = 0, times = ncol(x)) )
      } else {
        return ( rowSums( (t(x) - colMeans(x))^2 ) / ( dim(x)[1] - 1 ) )
      }
    } else stop("Please enter valid dimension, for rows, dim = 1; for colums, dim = 2")
  }

  # get names of columns that do not have only NA
  nonempty_column_names <-
    function(x) {
      # compute column sums; result is zero for columns having no non-NA values
      column_sum      <- sapply(1:ncol(x), function(i) sum(x[,i], na.rm = TRUE))
      
      # return names of columns 
      return ( as.character( colnames(x)[column_sum > 0.0] ) )
    }

  # produce matrix from matrix xpre where all rows and columns having zero variance have been removed 
  nonzero_var <- function(xpre) {
    nonzero_var_internal <- function(x) {
      if (nrow(x) == 0) {
          utils::str(x)
          stop("matrix has no rows")
      }
      if (ncol(x) == 0) {
          utils::str(x)
          stop("matrix has no columns")
      }
      if ( is.numeric(x) ) {
        # exclude any rows with zero variance
        row.vars <- MatVar(x, dim = 1)
        nonzero.row.vars <- row.vars > 0
        nonzero.rows <- row.vars[nonzero.row.vars]
        if ( length(rownames(x)) != length(rownames(nonzero.rows)) ) {
          row.names <- attr(nonzero.rows,"names")
          x <- x[ row.names, , drop = FALSE ]
        }

        # exclude any columns with zero variance
        column.vars <- MatVar(x, dim = 2)
        nonzero.column.vars <- column.vars > 0
        nonzero.columns <- column.vars[nonzero.column.vars]
        if ( length(colnames(x)) != length(colnames(nonzero.columns)) ) {
          column.names <- attr(nonzero.columns,"names")
          x <- x[ , column.names, drop = FALSE ]
        }
      }
      return (x)
    }

    # purge rows and columns that have zero variance until there are nor more changes
    #   rationale: nonzero_var_internal first purges rows with zero variance, then columns with zero variance,
    #              so there exists the possibility that a row's variance might change to zero when a column is removed;
    #              therefore, iterate till there are no more changes
    if ( is.numeric(xpre) ) {
      my.nrow <- 0
      my.ncol <- 0
      while ( nrow(xpre) != my.nrow || ncol(xpre) != my.ncol ) {
        my.nrow <- nrow(xpre)
        my.ncol <- ncol(xpre)
        xpre <- nonzero_var_internal(xpre)
      }
    }
    return (xpre)
  }
  # ...

  # return FALSE if any paths are exact duplicates
  my.paths <- c(dataMatrix_in, dataMatrix_out, sampleMetadata_in, sampleMetadata_out, variableMetadata_in, variableMetadata_out)
  if ( length(my.paths) != length(unique(my.paths)) ) {
    failure_action("some paths are duplicated")
    for ( my.path in my.paths ) {
      failure_action(my.path)
    }
    stop("some paths are duplicated")
    return (FALSE)
  }

  # ---
  # read in the sample metadata
  smpl_metadata_input_env <- read_data_frame(sampleMetadata_in, "sample metadata input")
  if (!smpl_metadata_input_env$success) {
    failure_action(smpl_metadata_input_env$msg)
    return ( FALSE )
  }
  smpl_metadata <- smpl_metadata_input_env$data

  # extract rownames
  rownames(smpl_metadata) <- smpl_metadata[,samplename_column]
  
  # select the first column of the rows indicated by classes, include, & class_column, but don't drop dimension
  selected_rows           <- smpl_metadata[ xor( !include, smpl_metadata[,class_column] %in% classes ), 1, drop = FALSE ]
  # obtain the row names
  sample_names            <- rownames( selected_rows )
  # ...

  # ---
  # read in the variable metadata
  vrbl_metadata_input_env <- read_data_frame(variableMetadata_in, "variable metadata input")
  if (!vrbl_metadata_input_env$success) {
    failure_action(vrbl_metadata_input_env$msg)
    return ( FALSE )
  }
  vrbl_metadata <- vrbl_metadata_input_env$data

  # extract rownames (using make.names to handle degenerate feature names)
  err.env <- new.env()
  err.env$success <- FALSE
  err.env$msg <- "no message setting vrbl_metadata rownames"
  tryCatch(
    expr = {
      rownames(vrbl_metadata) <- make.names( vrbl_metadata[,1], unique = TRUE )
      vrbl_metadata[,1] <- rownames(vrbl_metadata)
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
  data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  if (!data_matrix_input_env$success) {
    failure_action(data_matrix_input_env$msg)
    return ( FALSE )
  }
  data_matrix <- data_matrix_input_env$data

  # extract rownames (using make.names to handle degenerate feature names)
  err.env <- new.env()
  err.env$success <- FALSE
  err.env$msg <- "no message setting data_matrix rownames"
  tryCatch(
    expr = {
      rownames(data_matrix) <- make.names( data_matrix[,1], unique = TRUE )
      err.env$success     <- TRUE
    }
  , error = function(e) {
     err.env$ msg <- sprintf("failed to set rownames for data_matrix read because '%s'", e$message) 
    }
  )
  if (!err.env$success) {
    failure_action(err.env$msg)
    return ( FALSE )
  }
  # ...

  # remove rownames column
  data_matrix <- data_matrix[,2:ncol(data_matrix)]

  # select the subset of samples indicated by classes, include, & class_column
  data_matrix <- data_matrix[,intersect(sample_names,colnames(data_matrix)), drop = FALSE]

  # impute missing values with supplied or default method
  data_matrix <- data_imputation(data_matrix)

  # purge data_matrix of rows and columns that have zero variance
  data_matrix <- nonzero_var(data_matrix)

  # purge smpl_metadata and vrbl_metadata of irrelevant rows
  sample_names <- intersect(sample_names,colnames(data_matrix))
  variable_names <- intersect( rownames(vrbl_metadata), rownames(data_matrix) )

  # ---
  # write out the results
  err.env <- new.env()
  err.env$success <- FALSE
  err.env$msg <- "no message writing output files"
  tryCatch(
    expr = {
      utils::write.table( x = data_matrix  
                           [ rownames(data_matrix) %in% variable_names    # row selector
                           , colnames(data_matrix) %in% sample_names      # column selector
                           , drop = FALSE                                 # keep two dimensions
                           ]
                         , file = dataMatrix_out
                         , sep = "\t"
                         , quote = FALSE
                         , row.names = TRUE
                         , col.names = NA
                         )
      utils::write.table( x = smpl_metadata
                           [ sample_names                                 # row selector
                           ,                                              # column selector (select all)
                           , drop = FALSE                                 # keep two dimensions
                           ]
                         , file = sampleMetadata_out
                         , sep = "\t"
                         , quote = FALSE
                         , row.names = FALSE
                         )
      utils::write.table( x = vrbl_metadata
                           [ rownames(vrbl_metadata) %in% variable_names  # row selector
                           ,                                              # column selector (select all)
                           , drop = FALSE                                 # keep two dimensions
                           ]
                         , file = variableMetadata_out
                         , sep = "\t"
                         , quote = FALSE
                         , row.names = FALSE
                         )
      err.env$success     <- TRUE
    }
  , error = function(e) {
     err.env$ msg <- sprintf("failed to set write output files because '%s'", e$message) 
    }
  )
  if (!err.env$success) {
    failure_action(err.env$msg)
    return ( FALSE )
  }
  return (TRUE)
  # ...
}



  
