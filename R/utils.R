# UTILITY FUNCTIONS ----
# These are used to read and write the directories and feather files


## MISC ----
#' Make Directory
#'
#' Creates or overwrites a directory depending if it already exists or not
#'
#' @param fp File path
make_dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
}


## WRITING ----
### Write a sparse integer matrix to a feather file ----
#' Write Sparse Integer Matrix
#'
#' Write out a sparse matrix with integer elements as 3 separate
#' feather files:
#' 'x': the data written out as a coordinate matrix.
#' 'colnames': the names of the columns.
#' 'rownames': the names of the rows.
#'
#' @param mat A matrix with dimension names
#' @param path A file path to write to
#' @param compression The type of compression to use
#' @param compression_level The level of compression to use if "zstd" is specified
write_sparse_int <- function(mat, path, compression = "lz4", compression_level = NULL) {
  # Check for correct sparse matrix type !!
  if (class(mat) != "dgTMatrix") {
    message("Converting matrix of class ", class(mat), " to dgTMatrix(AKA a COOrdinate Matrix)")

    mat <- as(mat, "TsparseMatrix")
  }

  # Convert to arrow tables
  colnames <- arrow::arrow_table(colnames = colnames(mat))
  rownames <- arrow::arrow_table(rownames = rownames(mat))
  mat <- arrow::arrow_table(data = mat@x, i_index = mat@i, j_index = mat@j,
                            schema = arrow::schema( data = arrow::int64(), # Elements are integers
                                                    i_index = arrow::uint64(),
                                                    j_index = arrow::uint64()
                            )
  )


  # Add metadata
  mat$metadata$class <- c("class" = "COOmatrix")

  # Write out data
  arrow::write_feather(colnames, paste(path, 'colnames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(rownames, paste(path, 'rownames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(mat, paste(path, 'x', sep='/'),
                       compression = compression, compression_level = compression_level)
}


### Write a sparse float matrix to a feather file ----
#' Write Sparse Float Matrix
#'
#' Write out a sparse matrix with float elements as 3 separate
#' feather files:
#' 'x': the data written out as a coordinate matrix.
#' 'colnames': the names of the columns.
#' 'rownames': the names of the rows.
#'
#' @param mat A matrix with dimension names
#' @param path A file path to write to
#' @param compression The type of compression to use
#' @param compression_level The level of compression to use if "zstd" is specified
write_sparse_float <- function(mat, path, compression = "lz4",  compression_level = NULL) {
  # Check for correct sparse matrix type !!
  if (class(mat) != "dgTMatrix") {
    message("Converting matrix of class ", class(mat), " to dgTMatrix(AKA a COOrdinate Matrix)")

    mat <- as(mat, "TsparseMatrix")
  }

  # Convert to arrow table
  colnames <- arrow::arrow_table(colnames = colnames(mat))
  rownames <- arrow::arrow_table(rownames = rownames(mat))
  mat <- arrow::arrow_table(data = mat@x, i_index = mat@i, j_index = mat@j,
                            schema = arrow::schema(data = arrow::float64(), # Elements are floats
                                                   i_index = arrow::uint64(),
                                                   j_index = arrow::uint64()
                            )
  )

  # Add metadata
  mat$metadata$class <- c("class" = "COOmatrix")

  # Write out data
  arrow::write_feather(colnames, paste(path, 'colnames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(rownames, paste(path, 'rownames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(mat, paste(path, 'x', sep='/'),
                       compression = compression, compression_level = compression_level)
}


### Write a dense integer matrix to a feather file ----
#' Write Dense Integer Matrix
#'
#' Write out a dense matrix with integer elements as 3 separate
#' feather files:
#' 'x': the data written out as a single vector.
#' 'colnames': the names of the columns.
#' 'rownames': the names of the rows.
#'
#' @param mat A matrix with dimension names
#' @param path A file path to write to
#' @param compression The type of compression to use
#' @param compression_level The level of compression to use if "zstd" is specified
write_dense_int <- function(mat, path, compression = "lz4", compression_level = NULL) {
  # Check for dense matrix !!
  if (!setequal(class(mat), c("matrix", "array"))) {
    message("Converting matrix of class ", class(mat), " to dense matrix")

    mat <- as(mat, "matrix")
  }

  # Convert to arrow table
  colnames <- arrow::arrow_table(colnames = colnames(mat))
  rownames <- arrow::arrow_table(rownames = rownames(mat))
  mat <- arrow::arrow_table(data = mat,
                            schema = arrow::schema(data = arrow::int64() # Elements are integers
                            )
  )

  # Remove and add metadata
  mat$metadata <- NULL
  mat$metadata$class <- c("class" = "matrix")

  # Write
  arrow::write_feather(colnames, paste(path, 'colnames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(rownames, paste(path, 'rownames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(mat, paste(path, 'x', sep='/'),
                       compression = compression, compression_level = compression_level)
}


### Write a dense float matrix to a feather file ----
#' Write Dense Float Matrix
#'
#' Write out a dense matrix with float elements as 3 separate
#' feather files:
#' 'x': the data written out as a single vector.
#' 'colnames': the names of the columns.
#' 'rownames': the names of the rows.
#'
#' @param mat A matrix with dimension names
#' @param path A file path to write to
#' @param compression The type of compression to use
#' @param compression_level The level of compression to use if "zstd" is specified
write_dense_float <- function(mat, path, compression = "lz4", compression_level = NULL) {
  # Check for dense matrix !!
  if (!setequal(class(mat), c("matrix", "array"))) {
    message("Converting matrix of class ", class(mat), " to dense matrix")

    mat <- as(mat, "matrix")
  }

  # Convert to arrow table
  colnames <- arrow::arrow_table(colnames = colnames(mat))
  rownames <- arrow::arrow_table(rownames = rownames(mat))
  mat <- arrow::arrow_table(data = mat,
                            schema = arrow::schema(data = arrow::float64() # Elements are decimals
                            )
  )

  # Remove and add metadata
  mat$metadata <- NULL
  mat$metadata$class <- c("class" = "matrix")

  # Write
  arrow::write_feather(colnames, paste(path, 'colnames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(rownames, paste(path, 'rownames', sep='/'),
                       compression = compression, compression_level = compression_level)
  arrow::write_feather(mat, paste(path, 'x', sep='/'),
                       compression = compression, compression_level = compression_level)
}

### Write a dataframe to a feather file ----
#' Write Dataframe
#'
#' Write out a dataframe as a feather file:
#' 'x': the data written out as a table.
#'
#' @param df A dataframe
#' @param path A file path to write to
#' @param name_of_rows The name of the column that the rownames of the dataframe will have
#' @param compression The type of compression to use
#' @param compression_level The level of compression to use if "zstd" is specified
write_dataframe <- function(df, path, name_of_rows = "rownames", compression = "lz4", compression_level = NULL) {
  # Check for dataframe !!
  if (class(df) != "data.frame") {
    message("Converting dataframe of class ", class(df), " to data.frame")

    df <- as(df, "data.frame")
  }

  # Add row names to dataframe
  df <- cbind(rownames(df), df)
  colnames(df)[1] <- name_of_rows

  # Convert to arrow table
  df <- arrow::arrow_table(df)

  # Add metadata
  df$metadata$class <- c("class" = "data.frame")

  # Write
  arrow::write_feather(df, paste(path, 'x', sep='/'), compression = compression, compression_level = compression_level)
}


## READING ----
### Load feather file into R ----
#' Load Feather File
#'
#' Read a feather file, determine the objects class based on
#' its metadata, and return the recreated the object
#'
#' @param path File path
load_feather <- function(path) {
  # Check if object file exists !!
  if (!file.exists(paste(path, 'x', sep ='/'))) {
    stop("Feather file 'x' does not exist in directory: ", path ,"\nPlease ensure that you inputted the correct filepath.")
  }

  # Read feather file as an arrow table
  obj <- arrow::read_feather(paste(path, 'x', sep='/'), as_data_frame = FALSE)

  # Check class !!
  class_name <- obj$metadata$class
  if (is.null(class_name)) {
    stop("Loaded object does not have an associated 'class' metadata key(obj$metadata$class is NULL).")
  }

  # For matrices
  if (grepl("matrix", class_name)) {

    # Check to see if dimension names exist !!
    if (!file.exists(paste(path, "colnames", sep='/')) || !file.exists(paste(path, "rownames", sep='/'))) {
      stop("Dimension name files are not found with the matrix data.")
    }

    # Sparse matrices
    if (class_name == "COOmatrix") {

      # Load dimension names
      colnames <- arrow::read_feather(paste(path, "colnames", sep='/'))$colnames
      rownames <- arrow::read_feather(paste(path, "rownames", sep='/'))$rownames

      # Create sparse matrix
      obj <- as.data.frame(obj)
      obj <- Matrix::sparseMatrix(x = obj$data, i = obj$i_index, j = obj$j_index,
                                  dimnames = list(rownames, colnames), index1 = FALSE)

      return(obj)
    }

    # Dense matrices
    else if (class_name == "matrix") {

      # Load dimension names
      colnames <- arrow::read_feather(paste(path, "colnames", sep='/'))$colnames
      rownames <- arrow::read_feather(paste(path, "rownames", sep='/'))$rownames

      # Create dense matrix
      obj <- matrix(data = obj$data, length(rownames), length(colnames), dimnames = list(rownames, colnames))

      return(obj)
    }

    else {
      stop("Matrix does not have class 'COOmatrix' or 'matrix' as its class in the object's metadata.")
    }

  }

  else if (class_name == "data.frame") {
    # Convert to dataframe
    obj <- as(as.data.frame(obj), "data.frame") # Required twice because first conversion is to tibble

    # Add row names
    rownames(obj) <- obj[,1]

    # Return object
    return(obj)
  }

  else {
    stop("Object does not have a recognized class.")
  }


}
