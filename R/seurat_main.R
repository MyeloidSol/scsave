# SEURAT INTERACTION ----


## Create AssayV5 object ----
#' Create Assay V5 Object
#'
#' Given a named list of matrices, create a Seurat Assay V5 Object
#'
#' @param assay_data A named list of matrices
#'
#' @return Returns a Seurat Assay V5 Object
#' @export
create_assayv5 <- function(assay_data) {

  a5obj <- SeuratObject::CreateAssay5Object(counts = assay_data[[1]]) # Create temporary assay with initial layer

  # Unfortunately, Seurat does not allow you to change layer names so this set of code is needed to correct for naming
  if (names(assay_data)[1] != "counts") {
    init_layer <- names(assay_data)[1]

    a5obj[[init_layer]] <- assay_data[[1]] # Add data with correct layer name
    DefaultLayer(a5obj) <- init_layer

    a5obj$counts <- NULL # Remove old layer
  }

  assay_data[[1]] <- NULL; gc() # Free up space

  if (length(assay_data) != 0) { # If there are other layers, add them
    for (layer in names(assay_data)) {
      a5obj[[layer]] <- assay_data[[layer]]
    }
  }

  return(a5obj)
}

## Save Seurat object ----
#' Save Seurat Object
#'
#' Save a Seurat object to a directory
#' @param sobj A Seurat object
#' @param dir_path A file path to a directory
#' @param compression The type of compression to use, default is "zstd"
#' @param compression_level If compression "zstd" is specified, the level of compression may also be specified
#'
#' @return Outputs a folder with the file path specified
#' @export
save_seurat <- function(sobj, dir_path, compression = "lz4", compression_level = NULL) {
  # Update dir_path
  dir_path <- paste(dir_path, "scdata", sep = '/')

  # Create main folder
  make_dir(dir_path)


  ### Metadata ----
  message("Writing out cell metadata...")
  subdir_path <- paste(dir_path, "cell_metadata", sep = '/') # Path to current sub directory

  # Create metadata directory
  make_dir(subdir_path)

  # Write out cell metadata
  write_dataframe(sobj[[]], path = subdir_path, name_of_rows = "cell_names")


  ### Assays ----
  message("Writing out assays...")
  subdir_path <- paste(dir_path, "assays", sep = '/') # Path to current sub directory

  # Create assays directory
  make_dir(subdir_path)

  assay_names <- SeuratObject::Assays(sobj) # Grab assay names
  for (assay in assay_names) { # Iterate over all assays
    # Path to current assay
    assay_path <- paste(subdir_path, assay, sep ='/')

    # Create folder
    make_dir(assay_path)

    # Grab layer names
    layer_names <- SeuratObject::Layers(sobj[[assay]])
    for (layer in layer_names) {# Iterate over all layers in an assay
      # Path to current layer
      layer_path <- paste(assay_path, layer, sep ='/')

      # Create folder
      make_dir(layer_path)

      # Grab data
      data <- sobj[[assay]][[layer]]

      # Identify type of matrix
      if (is(data, "sparseMatrix")) { # Sparse matrix

        # Float elements(just checks first 100)
        if(any(data@x[1:100] %% 1 != 0)) {
          write_sparse_float(data, path = layer_path, compression, compression_level)
        }

        # Integer elements(just checks first 100)
        else if (all(data@x[1:100] %% 1 == 0)) {
          write_sparse_int(data, path = layer_path, compression, compression_level)
        }
      }

      else if (setequal(class(data), c("matrix", "array"))) { # Dense matrix
        # Float elements(just checks first 100)
        if(any(data[1:100] %% 1 != 0)) {
          write_dense_float(data, path = layer_path, compression, compression_level)
        }

        # Integer elements(just checks first 100)
        else if (all(data[1:100] %% 1 == 0)) {
          write_dense_int(data, path = layer_path, compression, compression_level)
        }

      }

      else {
        stop("Assay[", assay, "] Layer[", layer, "] does not contain a dense or sparse matrix.\nInstead data has class: ", class(data))
      }

    }
  }

}


## Load Seurat object ----
#' Load Seurat Object
#'
#' Loads a Seurat object that has been written with save_seurat
#'
#' @param dir_path file path to the saved Seurat object
#'
#' @return A Seurat object
#' @export
load_seurat <- function(dir_path) {
  # Check if main folder exists !!
  if (!file.exists(dir_path)) {
    stop("Directory does not exist! Please ensure you inputted the correct filepath.")
  }

  ### Metadata ----
  subdir_path <- paste(dir_path, "cell_metadata", sep = '/') # Path to current sub directory

  # Check if cell_metadata folder exists
  if (!file.exists(subdir_path)) {
    stop("Cell metadata does not exist!\nAll Seurat objects should have cell metadata written out with save_seurat()")
  }

  # Read cell metadata
  cell_metadata <- load_feather(subdir_path)

  # Format cell metadata
  rownames(cell_metadata) <- cell_metadata$cell_names
  cell_metadata$cell_names <- NULL

  ### Assays ----
  subdir_path <- paste(dir_path, "assays", sep='/')

  # Check if assays folder exists !!
  if (!file.exists(subdir_path)) {
    stop("Assays folder does not exist!")
  }

  # Create list to add data to
  tmp <- list()

  # Grab assay names
  assay_names <- list.files(subdir_path)
  for (assay in assay_names) { # Iterate over all assays
    # Path to current assay
    assay_path <- paste(subdir_path, assay, sep ='/')

    # Grab layer names
    layer_names <- list.files(assay_path)
    for (layer in layer_names) { # Iterate over all layers in an assay
      # Path to current layer
      layer_path <- paste(assay_path, layer, sep ='/')

      # Load data
      data <- load_feather(layer_path)

      # Append data to list
      tmp[[assay]][[layer]] <- data
    }
  }

  ### Create Seurat Object ----
  # Initialize with first assay
  sobj <- Seurat::CreateSeuratObject(create_assayv5(tmp[[1]]), assay = names(tmp)[1])
  tmp[[1]] <- NULL; gc() # Free up space

  # Add cell metadata
  sobj[[]] <- cell_metadata

  # Add other assays
  if (length(tmp) != 0) {
    for (assay in names(tmp)) {
      sobj[[assay]] <- create_assayv5(tmp[[1]])
      tmp[[1]] <- NULL; gc() # Free up space
    }
  }

  return(sobj)
}
