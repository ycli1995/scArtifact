
#' Read H5AD file
#'
#' Create a \code{\link[SingleCellExperiment]{SingleCellExperiment}} from a
#' H5AD file.
#'
#' @param path A H5AD file
#' @param name Name of the HDF5 group where the AnnData is stored.
#' @param layers Which \code{layers} to read. If \code{NULL}, read all layers.
#' @param obsm Which \code{obsm} to read. If \code{NULL}, read all \code{obsm}.
#' @param obsp,varp Which \code{obsp} (\code{varp}) to read. If \code{NULL},
#' read all \code{obsp} (\code{varp}).
#' @param uns Which slots of \code{uns} data to read. If \code{NULL}, read all
#' slots of \code{uns}.
#' @param raw Whether or not to read the \code{raw} of AnnData. By default is
#' \code{TRUE}.
#' @param use.BPCells Whether or not to use \pkg{BPCells} to hold the matrix.
#'
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment assay<- rowData
#' @importFrom SingleCellExperiment altExp<- colPair<- LinearEmbeddingMatrix
#' reducedDim<- rowPair<-
#' @importFrom hdf5r.Extra h5Read h5List
#' @export
readH5AD <- function(
    path,
    name = "/",
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    varp = NULL,
    uns = NULL,
    raw = NULL,
    use.BPCells = FALSE
) {
  path <- file_path_as_absolute(path)
  all.slots <- h5List(x = path, name = name)
  obs <- h5Read(x = path, name = file.path(name, "obs"))
  sce <- .read_h5ad_raw(
    path = path,
    name = name,
    obs = obs,
    use.BPCells = use.BPCells
  )
  var <- rowData(sce)
  all.layers <- h5List(path, name = file.path(name, "layers"))
  layers <- layers %||% all.layers
  layers <- intersect(layers, all.layers)
  for (i in layers) {
    assay(sce, i = i) <- .read_h5ad_layer(
      path = path,
      name = file.path(name, "layers", i),
      obs = obs,
      var = var,
      use.BPCells = use.BPCells
    )
  }
  all.obsm <- h5List(path, name = file.path(name, "obsm"))
  obsm <- obsm %||% all.obsm
  obsm <- intersect(obsm, all.obsm)
  for (i in obsm) {
    obsm <- h5Read(path, name = file.path(name, "obsm", i)) %>%
      t()
    rownames(obsm) <- colnames(sce)
    colnames(obsm) <- paste0(i, "_", seq(ncol(obsm)))
    featureLoadings <- matrix(0, nrow = 0, ncol = ncol(obsm))
    colnames(featureLoadings) <- colnames(obsm)
    reducedDim(sce, type = i) <- LinearEmbeddingMatrix(
      sampleFactors = obsm,
      featureLoadings = featureLoadings
    )
  }
  all.obsp <- h5List(path, name = file.path(name, "obsp"))
  obsp <- obsp %||% all.obsp
  obsp <- intersect(obsp, all.obsp)
  for (i in obsp) {
    colPair(sce, type = i) <- h5Read(path, name = file.path(name, "obsp", i))
  }
  all.varp <- h5List(path, name = file.path(name, "varp"))
  varp <- varp %||% all.varp
  varp <- intersect(varp, all.varp)
  for (i in varp) {
    rowPair(sce, type = i) <- h5Read(path, name = file.path(name, "varp", i))
  }
  all.uns <- h5List(path, name = file.path(name, "uns"))
  uns <- uns %||% all.uns
  uns <- intersect(uns, all.uns)
  for (i in uns) {
    metadata(sce)[[i]] <- h5Read(path, name = file.path(name, "uns", i))
  }
  raw <- raw %||% TRUE
  if (!"raw" %in% all.slots) {
    return(sce)
  }
  raw.slots <- h5List(path, name = file.path(name, "raw"))
  if (!(length(raw.slots) > 0 && is_true(raw))) {
    return(sce)
  }
  altExp(sce, e = "raw") <- .read_h5ad_raw(
    path = path,
    name = file.path(name, "raw"),
    obs = obs,
    use.BPCells = use.BPCells
  )
  return(sce)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom HDF5Array HDF5Array HDF5ArraySeed
#' @importFrom S4Vectors new2
#' @importClassesFrom HDF5Array Dense_H5ADMatrixSeed
.read_h5ad_layer_dense <- function(path, name, obs, var) {
  dimnames <- list(rownames(var), rownames(obs))
  ans0 <- HDF5ArraySeed(filepath = path, name = name)
  if (length(dim(ans0)) == 2L) {
    return(new2("Dense_H5ADMatrixSeed", ans0, dimnames = dimnames))
  }
  warning(
    "HDF5 dataset '", name, "' in file '",
    path, "' does not have exactly 2 dimensions. ",
    "Using 'HDF5Array' to access this dataset.",
    immediate. = TRUE, call. = FALSE
  )
  mat <- HDF5Array(filepath = ans0)
  dimnames(mat) <- dimnames
  return(mat)
}

#' @importFrom HDF5Array HDF5Array
#' @importFrom hdf5r.Extra h5Read is.H5Group
.read_h5ad_layer <- function(path, name, obs, var, use.BPCells = FALSE) {
  dimnames <- list(rownames(var), rownames(x = obs))
  if (use.BPCells && !is.H5Group(file = path, name = name)) {
    warning(
      "'BPCells' doesn't support dense matrix: ", name,
      "\nUse 'H5ADMatrix' instead.",
      immediate. = TRUE, call. = FALSE
    )
    return(.read_h5ad_layer_dense(
      path = path,
      name = name,
      obs = obs,
      var = var
    ))
  }
  if (use.BPCells) {
    mat <- open_matrix_anndata_hdf5(path = path, group = name)
    dimnames(mat) <- dimnames
    return(write_matrix_dir(mat = mat, dir = tempfile("anndata_matrix_h5")))
  }
  mat <- h5Read(x = path, name = name)
  dimnames(mat) <- dimnames
  return(mat)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom hdf5r.Extra h5List h5Read
#' @importClassesFrom DelayedArray DelayedArray
.read_h5ad_raw <- function(path, name, obs, use.BPCells = FALSE) {
  var <- h5Read(path, name = file.path(name, "var"))
  mat <- .read_h5ad_layer(
    path = path,
    name = file.path(name, "X"),
    obs = obs,
    var = var,
    use.BPCells = use.BPCells
  )
  if (inherits(x = mat, what = "DelayedArray")) {
    use.BPCells <- FALSE
  }
  return(SingleCellExperiment(
    assays = list(X = mat),
    rowData = as(var, "DataFrame"),
    colData = as(obs, "DataFrame")
  ))
}
