#' @title
#' Impute values for W4M data matrix
#'
#' @description
#' Substitute zero for missing or negative values in W4M data matrix
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
  function(m) {
    # replace NA values with zero
    m[is.na(m)] <- 0
    # replace negative values with zero, if applicable (It should never be applicable!)
    m[m<0] <- 0
    # return matrix as the result
    return (m)
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
    dim.x.2 <- dim(m)[2]
    if ( dim.x.2 == 0 )
      stop("w4m__var_by_rank_or_file: there are zero columns")
    if ( dim.x.2 == 1 ) {
      stop("w4m__var_by_rank_or_file: a single column is insufficient to calculate a variance")
    }
  }
  else if (dim == 2) {
    dim.x.1 <- dim(m)[1]
    if ( dim.x.1 == 0 ) {
      stop("w4m__var_by_rank_or_file: there are zero rows")
    }
    if ( dim.x.1 == 1 ) {
      stop("w4m__var_by_rank_or_file: a single row is insufficient to calculate a variance")
    }
    m <- t(m)
  }
  else {
    stop("w4m__var_by_rank_or_file: dim is invalid; for rows, use dim = 1; for colums, use dim = 2")
  }
  return(rowSums((m - rowMeans(m))^2)/(dim(m)[2] - 1))
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
    if ( is.numeric(x) ) {
      # exclude any rows with zero variance
      row.vars <- w4m__var_by_rank_or_file(x, dim = 1)
      nonzero.row.vars <- row.vars > 0
      nonzero.rows <- row.vars[nonzero.row.vars]
      if ( length(rownames(x)) != length(rownames(nonzero.rows)) ) {
        row.names <- attr(nonzero.rows,"names")
        x <- x[ row.names, , drop = FALSE ]
      }

      # exclude any columns with zero variance
      column.vars <- w4m__var_by_rank_or_file(x, dim = 2)
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
#' Filter a set of retention-corrected W4M files (dataMatrix, sampleMetadata, variableMetadata) by sample-class
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
#'   , classes = c("blankpos")
#'   , include = FALSE
#'   )
#' }
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
  
  # get names of columns that do not have only NA
  nonempty_column_names <-
    function(x) {
      # compute column sums; result is zero for columns having no non-NA values
      column_sum      <- sapply(1:ncol(x), function(i) sum(x[,i], na.rm = TRUE))
      
      # return names of columns 
      return ( as.character( colnames(x)[column_sum > 0.0] ) )
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

  # read in the sample metadata
  if ( is.character(sampleMetadata_in) ){
    smpl_metadata_input_env <- read_data_frame(sampleMetadata_in, "sample metadata input")
    if (!smpl_metadata_input_env$success) {
      failure_action(smpl_metadata_input_env$msg)
      return ( FALSE )
    }
    smpl_metadata <- smpl_metadata_input_env$data
  } else if ( is.environment(smpl_metadata) ||  is.list(smpl_metadata) ) {
    if ( Reduce(`|`,names(sampleMetadata_in) == "sampleMetadata") ) {
      smpl_metadata <- as.list(sampleMetadata_in)$sampleMetadata
    }
    else {
      stop("sampleMetadata_in has no member 'sampleMetadata'")
      return (FALSE)
    }
  } else if ( is.data.frame(sampleMetadata_in) ) {
    smpl_metadata <- sampleMetadata_in
  } else {
    stop(sprintf("sampleMetadata_in has unexpected type %s"), typeof(sampleMetadata_in))
    return (FALSE)
  }
  

  # extract rownames
  rownames(smpl_metadata) <- smpl_metadata[,samplename_column]
  
  if (nchar(class_column) > 0 && length(classes) > 0) {
    # select the first column of the rows indicated by classes, include, & class_column, but don't drop dimension
    #selected_rows <- smpl_metadata[ xor( !include, smpl_metadata[,class_column] %in% classes ), 1, drop = FALSE ]
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
            grepl(pattern = pattern, x = smpl_metadata[,class_column])
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

  # read in the variable metadata
  if ( is.character(variableMetadata_in) ){
    vrbl_metadata_input_env <- read_data_frame(variableMetadata_in, "variable metadata input")
    if (!vrbl_metadata_input_env$success) {
      failure_action(vrbl_metadata_input_env$msg)
      return ( FALSE )
    }
    vrbl_metadata <- vrbl_metadata_input_env$data
  } else if ( is.environment(vrbl_metadata) ||  is.list(vrbl_metadata) ) {
    if ( Reduce(`|`,names(variableMetadata_in) == "variableMetadata") ) {
      vrbl_metadata <- as.list(variableMetadata_in)$variableMetadata
    }
    else {
      stop("variableMetadata_in has no member 'variableMetadata'")
      return (FALSE)
    }
  } else if ( is.data.frame(variableMetadata_in) ) {
    vrbl_metadata <- variableMetadata_in
  } else {
    stop(sprintf("variableMetadata_in has unexpected type %s"), typeof(variableMetadata_in))
    return (FALSE)
  }
  

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

  # read in the data matrix
  if ( is.character(dataMatrix_in) ){
    data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
    if (!data_matrix_input_env$success) {
      failure_action(data_matrix_input_env$msg)
      return ( FALSE )
    }
    data_matrix <- data_matrix_input_env$data
  } else if ( is.environment(data_matrix) ||  is.list(data_matrix) ) {
    if ( Reduce(`|`,names(dataMatrix_in) == "dataMatrix") ) {
      data_matrix <- as.list(dataMatrix_in)$dataMatrix
    }
    else {
      stop("dataMatrix_in has no member 'dataMatrix'")
      return (FALSE)
    }
  } else if ( is.data.frame(dataMatrix_in) || is.matrix(dataMatrix_in) ) {
    data_matrix <- as.matrix(dataMatrix_in)
  } else {
    stop(sprintf("dataMatrix_in has unexpected type %s"), typeof(dataMatrix_in))
    return (FALSE)
  }

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
     err.env$msg <- sprintf("failed to set rownames for data_matrix read because '%s'", e$message) 
    }
  )
  if (!err.env$success) {
    failure_action(err.env$msg)
    return ( FALSE )
  }

  # remove rownames column
  data_matrix <- data_matrix[,2:ncol(data_matrix)]

  # select the subset of samples indicated by classes, include, & class_column
  data_matrix <- data_matrix[,intersect(sample_names,colnames(data_matrix)), drop = FALSE]

  # impute missing values with supplied or default method
  data_matrix <- data_imputation(data_matrix)

  # convert data_matrix to matrix from data.frame
  data_matrix <- as.matrix(data_matrix)

  # purge data_matrix of rows and columns that have zero variance
  data_matrix <- w4m__nonzero_var(data_matrix)

  # purge smpl_metadata and vrbl_metadata of irrelevant rows
  sample_names <- intersect(sample_names,colnames(data_matrix))
  variable_names <- intersect( rownames(vrbl_metadata), rownames(data_matrix) )

  # write out the results
  err.env <- new.env()
  err.env$success <- FALSE
  err.env$msg <- "no message writing output files"
  tryCatch(
    expr = {
      # write the data matrix
      if ( is.character(dataMatrix_out) ){
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
      } else if ( is.environment(dataMatrix_out) ||  is.list(dataMatrix_out) ) {
        dataMatrix_out$dataMatrix <- data_matrix
      } else {
        stop(sprintf("dataMatrix_out has unexpected type %s"), typeof(dataMatrix_out))
        return (FALSE)
      }

      # write the sample metadata
      if ( is.character(sampleMetadata_out) ){
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
      } else if ( is.environment(sampleMetadata_out) ||  is.list(sampleMetadata_out) ) {
        sampleMetadata_out$sampleMetadata <- smpl_metadata
      } else {
        stop(sprintf("sampleMetadata_out has unexpected type %s"), typeof(sampleMetadata_out))
        return (FALSE)
      }

      # write the variable metadata
      if ( is.character(variableMetadata_out) ){
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
      } else if ( is.environment(variableMetadata_out) ||  is.list(variableMetadata_out) ) {
        variableMetadata_out$variableMetadata <- vrbl_metadata
      } else {
        stop(sprintf("variableMetadata_out has unexpected type %s"), typeof(variableMetadata_out))
        return (FALSE)
      }

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



  
