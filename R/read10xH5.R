
#' Read 10X hdf5 file
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#'
#' @param path Path to HDF5 file]
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param use.BPCells Whether or not to use \pkg{BPCells} to hold the matrix.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom hdf5r.Extra h5Read h5List
#' @importFrom Matrix sparseMatrix
#' @importFrom rlang is_list
#' @export
read10xH5 <- function(path, use.names = TRUE, use.BPCells = FALSE) {
  path <- file_path_as_absolute(path)
  groups <- h5List(x = path)
  read_features_func <- ifelse(
    test = "matrix" %in% groups,
    yes = .read_10x_h5_features_v3,
    no = .read_10x_h5_features_v2
  )
  output <- list()
  for (genome in groups) {
    features <- read_features_func(path, genome)
    barcodes <- h5Read(path, name = file.path(genome, "barcodes"))
    mat <- open_matrix_10x_hdf5(path = path)
    if (is_list(mat)) {
      mat <- mat[[genome]]
    }
    if (!use.BPCells) {
      mat <- as(mat, "dgCMatrix")
    }
    if (use.names) {
      rownames(mat) <- make.unique(names = features$Name)
    } else {
      rownames(mat) <- make.unique(names = features$ID)
    }
    colnames(x = mat) <- barcodes
    if (use.BPCells) {
      mat <- write_matrix_dir(mat = mat, dir = tempfile("tenx_matrix_h5"))
    }
    mat <- SingleCellExperiment(
      assays = list(counts = mat),
      rowData = as(features, "DataFrame")
    )
    if ("Type" %in% colnames(features)) {
      types <- unique(features$Type)
      if (length(types) > 1) {
        message(
          "Genome ", genome, " has multiple modalities, ",
          "returning a list of SingleCellExperiment for this genome."
        )
        mat <- sapply(
          X = types,
          FUN = function(x) mat[features$Type == x, ],
          simplify = FALSE,
          USE.NAMES = TRUE
        )
        gc(verbose = FALSE)
      }
    }
    output[[genome]] <- mat
  }
  if (length(output) == 1) {
    return(output[[genome]])
  }
  return(output)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom hdf5r.Extra h5Read
#' @importFrom dplyr select
.read_10x_h5_features_v3 <- function(path, name) {
  features <- h5Read(path, name = file.path(name, "features"))
  features[["_all_tag_keys"]] <- NULL
  features <- as.data.frame(features, stringsAsFactors = FALSE)
  query.cols <- c(
    ID = "id",
    Name = "name",
    Type = "feature_type",
    genome = "genome"
  )
  found.cols <- query.cols[query.cols %in% colnames(features)]
  features <- features %>%
    select(found.cols)
  return(features)
}

#' @importFrom hdf5r.Extra h5Read
.read_10x_h5_features_v2 <- function(path, name) {
  data.frame(
    ID = h5Read(path, name = file.path(name, "genes")),
    Name = h5Read(path, name = file.path(name, "gene_names")),
    stringsAsFactors = FALSE
  )
}
