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
#' @param name The name of the file that will contain the data
#' @param compression The type of compression to use, default is "zstd"
#' @param compression_level If compression "zstd" is specified, the level of compression may also be specified
#'
#' @return Outputs a folder with the file path specified
#' @export
save_seurat <- function(sobj, dir_path = getwd(), name= "scdata", compression = "lz4", compression_level = NULL) {
  # Update dir_path
  dir_path <- paste(dir_path, name, sep = '/')

  # Create main folder
  make_dir(dir_path)

  ### Metadata ----
  message("Writing out metadata...")
  #### Cell Metadata ----
  subdir_path <- paste(dir_path, "cell_metadata", sep = '/') # Path to current sub directory

  # Create cell metadata directory
  make_dir(subdir_path)

  # Write out cell metadata
  write_dataframe(sobj[[]], path = subdir_path, name_of_rows = "cell_names")


  ### Assays ----
  message("Writing out assays...")
  subdir_path <- paste(dir_path, "assays", sep = '/') # Path to current sub directory

  # Create assays directory
  make_dir(subdir_path)

  #### Assay Data ----
  assay_names <- SeuratObject::Assays(sobj) # Grab assay names
  for (assay in assay_names) { # Iterate over all assays
    # Path to current assay
    assay_path <- paste(subdir_path, assay, sep ='/')

    # Create folder
    make_dir(assay_path)

    # Grab layer names
    layer_names <- SeuratObject::Layers(sobj[[assay]])
    for (layer in layer_names) { # Iterate over all layers in an assay
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
          write_sparse_float(data, path = layer_path,
                             compression, compression_level)
        }

        # Integer elements(just checks first 100)
        else if (all(data@x[1:100] %% 1 == 0)) {
          write_sparse_int(data, path = layer_path,
                           compression, compression_level)
        }
      }

      else if (setequal(class(data), c("matrix", "array"))) { # Dense matrix
        # Float elements(just checks first 100)
        if(any(data[1:100] %% 1 != 0)) {
          write_dense_float(data, path = layer_path,
                            compression, compression_level)
        }

        # Integer elements(just checks first 100)
        else if (all(data[1:100] %% 1 == 0)) {
          write_dense_int(data, path = layer_path,
                          compression, compression_level)
        }

      }

      else {
        stop("Assay[", assay, "] Layer[", layer, "] does not contain a dense or sparse matrix.\nInstead data has class: ", class(data))
      }

    }

    #### Assay Metadata ----
    # Create folder for assay metadata
    make_dir(paste(assay_path, "assay_metadata", sep ='/'))

    # Write assay metadata
    write_dataframe(df = sobj[[assay]]@meta.data, path = paste(assay_path, "assay_metadata", sep ='/'),
                    name_of_rows = "gene_names", compression, compression_level)

  }


  ### Dimension Reductions ----
  reduc_names <- Seurat::Reductions(sobj)

  if (!is.null(reduc_names)) { # Save if it exists
    message("Writing out reductions...")
    subdir_path <- paste(dir_path, "reductions", sep = '/') # Path to current sub directory

    # Create reductions directory
    make_dir(subdir_path)

    for (reduc in reduc_names) { # Iterate through all reductions
      # Path to current reduction
      reduc_path <- paste(subdir_path, reduc, sep ='/')

      # Create folder
      make_dir(reduc_path)

      # Grab reduction
      data <- sobj[[reduc]]


      # Path to cell embeddings
      embeddings_path <- paste(reduc_path, "embeddings", sep = '/')

      # Create folder
      make_dir(embeddings_path)

      # Write out cell embeddings
      write_dense_float(Seurat::Embeddings(data), embeddings_path)


      # Write out feature loadings(if they exist)
      if (length(Seurat::Loadings(data)) != 0) {
        # Path to feature loadings
        loadings_path <- paste(reduc_path, "loadings", sep = '/')

        # Create folder
        make_dir(loadings_path)

        # Write out feature loadings
        write_dense_float(Seurat::Loadings(data), loadings_path)
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
  message("Reading in metadata...")

  #### Cell Metadata ----
  subdir_path <- paste(dir_path, "cell_metadata", sep = '/') # Path to current sub directory

  # Check if cell_metadata folder exists !!
  if (!file.exists(subdir_path)) {
    stop("Cell metadata directory does not exist!")
  }

  # Read cell metadata
  cell_metadata <- load_feather(subdir_path)

  # Format cell metadata
  rownames(cell_metadata) <- cell_metadata$cell_names
  cell_metadata$cell_names <- NULL


  ### Assays ----
  message("Reading in assays...")
  subdir_path <- paste(dir_path, "assays", sep='/')

  # Check if assays folder exists !!
  if (!file.exists(subdir_path)) {
    stop("Assays directory does not exist!")
  }

  # Create list to add data to
  assay_tmp <- list()
  assay_md_tmp <- list()

  # Grab assay names
  assay_names <- list.files(subdir_path)

  for (assay in assay_names) { # Iterate over all assays
    # Path to current assay
    assay_path <- paste(subdir_path, assay, sep ='/')


    #### Assay Data ----
    # Grab layer names
    layer_names <- list.files(assay_path)
    layer_names <- layer_names[layer_names != "assay_metadata"] # Remove assay_metadata
    for (layer in layer_names) { # Iterate over all layers in an assay
      # Path to current layer
      layer_path <- paste(assay_path, layer, sep ='/')

      # Load data
      data <- load_feather(layer_path)

      # Append data to list
      assay_tmp[[assay]][[layer]] <- data
    }


    #### Assay Metadata ----
    # Load metadata
    data <- load_feather(paste(assay_path, "assay_metadata", sep = '/'))

    # Get rid of gene_names column
    rownames(data) <- data$gene_names
    data$gene_names <- NULL; gc()

    # Append metadata to list
    assay_md_tmp[[assay]] <- data
  }


  ### Dimension Reductions ----
  subdir_path <- paste(dir_path, "reductions", sep='/') # Path to current sub directory

  # Create list to add data to
  reduc_tmp <- list()

  if (file.exists(subdir_path)) {
    message("Reading in reductions...")

    # Grab reduction names
    reduc_names <- list.files(subdir_path)
    for (reduc in reduc_names) { # Iterate through all reductions
      # Path to current reduction
      reduc_path <- paste(subdir_path, reduc, sep ='/')

      # Path to cell embeddings
      embeddings_path <- paste(reduc_path, "embeddings", sep = '/')

      # Check if cell embeddings exist !!
      if (!file.exists(embeddings_path)) {
        stop("Cell embeddings do not exist for: ", reduc, " reduction!")
      }

      # Grab cell embeddings
      reduc_tmp[[reduc]]$embeddings <- load_feather(embeddings_path)

      # Path to feature loadings
      loadings_path <- paste(reduc_path, "loadings", sep = '/')

      reduc_tmp[[reduc]]$loadings <- new(Class = "matrix")
      # Check if feature loadings exists
      if (file.exists(loadings_path)) {
        reduc_tmp[[reduc]]$loadings <- load_feather(loadings_path)
      }

    }
  }


  ### Create Seurat Object ----
  message("Generating Seurat object...")
  # Initialize with first assay
  sobj <- Seurat::CreateSeuratObject(create_assayv5(assay_tmp[[1]]), assay = names(assay_tmp)[1])
  sobj[[names(assay_tmp)[1]]]@meta.data <- assay_md_tmp[[1]]

  assay_tmp[[1]] <- NULL; gc() # Free up space
  assay_md_tmp[[1]] <- NULL; gc() # Free up space

  # Add cell metadata
  sobj[[]] <- cell_metadata

  # Add other assays(if there are others)
  if (length(assay_tmp) != 0) {
    for (assay in names(assay_tmp)) {
      sobj[[assay]] <- create_assayv5(assay_tmp[[1]])
      assay_tmp[[1]] <- NULL; gc() # Free up space

      sobj[[assay]]@meta.data <- assay_md_tmp[[1]]
      assay_md_tmp[[1]] <- NULL; gc() # Free up space
    }
  }

  # Add dimension reductions(if there are any)
  if (length(reduc_tmp) != 0) {

    for (reduc in names(reduc_tmp)) {
      # Slot dimension reduced object
      sobj[[reduc]] <- Seurat::CreateDimReducObject(embeddings = reduc_tmp[[reduc]]$embeddings,
                                                    loadings = reduc_tmp[[reduc]]$loadings,
                                                    assay = Seurat::Assays(sobj)[1]) # Use first assay by default
    }


  }

  return(sobj)
}
