
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## prepH5AD ####################################################################

#' @param X Which `r .doc_links("assay")` to be used as `adata.X`. Default is
#' the first one.
#' @param layers Which `r .doc_links("assays")` to be used in `adata.layers`.
#' Default will be all assays except the `X` one.
#' @param obsm Which `r .doc_links("reducedDim")` to be used in `adata.obsm`.
#' Default will be all reduction data.
#' @param obsp Which `r .doc_links("colPairs")` to be used in `adata.obsp`.
#' @param varp Which `r .doc_links("rowPairs")` to be used in `adata.varp`.
#' @param uns Which `r .doc_links("metadata")` to be used in `adata.uns`.
#' @param raw Which `r .doc_links("altExp")` to be used in `adata.raw`.
#' The `nrows()` of `raw` must be larger than the `mainExp`.
#'
#' @importFrom SummarizedExperiment assayNames colData
#' @importFrom SingleCellExperiment altExp altExpNames colPair colPairNames
#' reducedDim rowPair rowPairNames
#'
#' @rdname exportH5AD
#' @export
#' @method prepH5AD SingleCellExperiment
prepH5AD.SingleCellExperiment <- function(
    object,
    X = NULL,
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    varp = NULL,
    uns = NULL,
    raw = NULL,
    ...
) {
  X <- X %||% assayNames(object)[1]
  out <- .prep_sce2h5ad_raw(object = object, X = X)

  layers <- layers %||% assayNames(object)
  layers <- layers %>%
    setdiff(X)
  out$layers <- .prep_sce2h5ad_layers(object = object, layers = layers)

  out$obs <- colData(object) %>%
    as.data.frame()

  obsm <- obsm %||% reducedDimNames(object)
  obsm <- obsm %>%
    intersect(reducedDimNames(object))
  out$obsm <- list()
  for (i in obsm) {
    out$obsm[[i]] <- reducedDim(object, type = i, withDimnames = TRUE) %>%
      as.matrix() %>%
      t()
  }

  obsp <- obsp %||% colPairNames(object)
  obsp <- obsp %>%
    intersect(colPairNames(object))
  out$obsp <- list()
  for (i in obsp) {
    out$obsp[[i]] <- colPair(object, type = i, asSparse = FALSE) %>%
      as("dgCMatrix")
  }
  varp <- varp %||% rowPairNames(object)
  varp <- varp %>%
    intersect(rowPairNames(object))
  out$varp <- list()
  for (i in varp) {
    out$varp[[i]] <- rowPair(object, type = i, asSparse = FALSE) %>%
      as("dgCMatrix")
  }
  gc(verbose = FALSE)
  uns <- uns %||% names(metadata(object))
  uns <- metadata(object) %>%
    names() %>%
    intersect(x = uns)
  out$uns <- metadata(object)[uns]

  raw <- raw[1] %||% altExpNames(object)[1]
  raw <- raw %>%
    intersect(altExpNames(object))
  if (length(raw) == 0) {
    return(out)
  }
  raw_exp <- altExp(object, e = raw)
  if (any(rownames(object) %in% rownames(raw_exp))) {
    out$raw <- .prep_sce2h5ad_raw(object = raw_exp, X = raw_X)
  }
  gc(verbose = FALSE)
  return(out)
}

#' @param assayt Which Seurat `assay` to use. If `NULL`, use
#' `DefaultAssay(object)`.
#'
#' @importFrom SeuratObject DefaultAssay Embeddings Graphs JoinLayers Layers
#' Misc Reductions
#' @importClassesFrom SeuratObject Assay5
#'
#' @rdname exportH5AD
#' @export
#' @method prepH5AD Seurat
prepH5AD.Seurat <- function(
    object,
    assay = NULL,
    X = NULL,
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    uns = NULL,
    raw = NULL
) {
  assay <- assay %||% DefaultAssay(object)
  X <- X %||% "data"
  assay.obj <- object@assays[[assay]]
  if (!inherits(assay.obj, "Assay5")) {
    assay.obj <- as(assay.obj, "Assay5")
  }
  assay <- JoinLayers(assay.obj)
  out <- .prep_seurat2h5ad_raw(assay = assay.obj, X = X)

  all.cells <- colnames(out$X)
  all.features <- rownames(out$X)
  all.layers <- Layers(assay.obj)

  raw <- raw[1] %||% all.layers[1]
  raw <- raw %>%
    intersect(all.layers)
  all.layers <- all.layers %>%
    setdiff(X) %>%
    setdiff(raw)
  if (length(raw) > 0 & !identical(rownames(assay.obj), rownames(out$X))) {
    out$raw <- .prep_seurat2h5ad_raw(assay = assay.obj, X = raw)
  }

  layers <- layers %||% all.layers
  out$layers <- .prep_seurat2h5ad_layers(assay = assay.obj, layers = layers)
  for (i in seq_along(out$layers)) {
    chk.cells <- identical(colnames(out$layers[[i]]), all.cells)
    chk.features <- identical(rownames(out$layers[[i]]), all.features)
    if (!all(chk.cells, chk.features)) {
      out$layers[[i]] <- out$layers[[i]][all.features, all.cells, drop = FALSE]
    }
  }
  out$obs <- obj[[]][all.cells, , drop = FALSE] %>%
    as.data.frame()

  all.reducs <- Reductions(object)
  obsm <- obsm %||% all.reducs
  obsm <- obsm %>%
    intersect(all.reducs)
  out$obsm <- list()
  for (i in obsm) {
    out$obsm[[i]] <- Embeddings(object, reduction = i)[all.cells, ] %>%
      as.matrix() %>%
      t()
  }

  all.graphs <- Graphs(object)
  obsp <- obsp %||% all.graphs
  obsp <- obsp %>%
    intersect(all.graphs)
  out$obsp <- list()
  for (i in obsp) {
    out$obsp[[i]] <- object[[i]] %>%
      as("dgCMatrix")
    if (!identical(colnames(out$obsp[[i]]), all.cells)) {
      out$obsp[[i]] <- out$obsp[[i]][all.cells, all.cells]
    }
    gc(verbose = FALSE)
  }

  uns <- uns %||% names(Misc(object))
  uns <- uns %>%
    intersect(names(Misc(object)))
  out$uns <- Misc(object)[uns]

  gc(verbose = FALSE)
  return(out)
}

## exportH5AD ##################################################################

#' @param file Path to the H5AD file.
#' @param name Name of the HDF5 group where to save \code{object}.
#' @param overwrite Whether or not to overwrite the existing HDF5 link.
#' @param gzip_level Enable zipping at the level given here.
#' @param verbose `r .vb_param`
#'
#' @importFrom hdf5r.Extra h5AbsLinkName
#'
#' @rdname exportH5AD
#' @export
#' @method exportH5AD SingleCellExperiment
exportH5AD.SingleCellExperiment <- function(
    object,
    file,
    name = "/",
    overwrite = FALSE,
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  name <- h5AbsLinkName(name = name)
  out.file <- normalizePath(path = file, mustWork = FALSE)
  file <- .h5_overwrite_before(
    file = file,
    name = name,
    overwrite = overwrite
  )
  if (!identical(x = file, y = out.file)) {
    on.exit(expr = unlink(x = file, force = TRUE))
  }
  object <- prepH5AD(
    object = object,
    X = X,
    layers = layers,
    obsm = obsm,
    obsp = obsp,
    varp = varp,
    uns = uns,
    raw = raw
  )
  .write_h5ad(
    object = object,
    file = file,
    name = name,
    gzip_level = gzip_level,
    verbose = verbose
  )
  file <- file_path_as_absolute(x = file)
  out.file <- out.file %>%
    file_path_as_absolute() %>%
    .h5_overwrite_after(file = file, name = name)
  return(invisible(x = NULL))
}

#' @param X Which \code{\link[SummarizedExperiment]{assay}} to be used as
#' `adata.X`. Default is the first one.
#' @param layers Which \code{\link[SummarizedExperiment]{assays}} to be used in
#' `adata.layers`. Default will be all assays except the `X` one.
#' @param obsm Which \code{\link[SingleCellExperiment]{reducedDim}} to be used
#' in `adata.obsm`. Default will be all reduction data.
#' @param obsp Which \code{\link[SingleCellExperiment]{colPair}} to be used in
#' `adata.obsp`.
#' @param varp Which \code{\link[SingleCellExperiment]{rowPair}} to be used in
#' `adata.varp`.
#' @param uns Which \code{\link[S4Vectors]{metadata}} to be used in `adata.uns`.
#' @param raw Which \code{\link[SingleCellExperiment]{altExp}} to be used in
#' `adata.raw`. The `nrows()` of `raw` must be larger than the `mainExp`.
#'
#' @importClassesFrom SeuratObject Seurat
#'
#' @export
#' @rdname exportH5AD
#' @method exportH5AD Seurat
exportH5AD.Seurat <- function(
    object,
    file,
    name = "/",
    overwrite = FALSE,
    gzip_level = 0L,
    verbose = TRUE,
    assay = NULL,
    X = NULL,
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    uns = NULL,
    ...
) {
  name <- h5AbsLinkName(name = name)
  out.file <- normalizePath(path = file, mustWork = FALSE)
  file <- .h5_overwrite_before(
    file = file,
    name = name,
    overwrite = overwrite
  )
  if (!identical(x = file, y = out.file)) {
    on.exit(expr = unlink(x = file, force = TRUE))
  }
  object <- prepObj(
    object = object,
    assay = assay,
    X = X,
    layers = layers,
    obsm = obsm,
    obsp = obsp,
    uns = uns
  )
  .write_h5ad(
    object = object,
    file = file,
    name = name,
    gzip_level = gzip_level,
    verbose = verbose
  )
  file <- file_path_as_absolute(x = file)
  out.file <- out.file %>%
    file_path_as_absolute() %>%
    .h5_overwrite_after(file = file, name = name)
  return(invisible(x = NULL))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Handle overwrite ############################################################

#' Check 'overwrite' before writing SCE into HDF5 file
#'
#' @return A character specifying the target file where to be actually written.
#'
#' @details
#' \itemize{
#' \item First check whether `file` already exists. If not, return the original
#' `file` directly.
#' \item When `file` exists, check whether `name` exists. If not, return the
#' orignial `file`.
#' \item When both `file` and `name` exist, check `overwrite`. If `FALSE`, raise
#' an error. When actually do `overwrite`, generate a temporary HDF5 file name
#' (without creating file) to hold the data.
#' }
#'
#' @noRd
#' @importFrom hdf5r.Extra h5AbsLinkName h5Backup h5CreateFile h5Exists
.h5_overwrite_before <- function(file, name, overwrite) {
  name <- h5AbsLinkName(name)
  if (!file.exists(file)) {
    return(normalizePath(file, mustWork = FALSE))
  }
  file <- normalizePath(file)
  tmp.file <- tempfile(tmpdir = dirname(file), fileext = ".h5")
  if (!h5Exists(file, name)) {
    message(
      "HDF5 file already exists. Use temporary file to write:",
      "\n  File: ", file,
      "\n  Object: ", name,
      "\n  Temporary file: ", tmp.file
    )
    return(tmp.file)
  }
  if (!overwrite) {
    stop(
      "\nFound object that already exists: ",
      "\n  File: ", file,
      "\n  Object: ", name,
      "\nSet 'overwrite = TRUE' to remove it."
    )
  }
  message(
    "HDF5 object already exists. Use temporary file to overwrite:",
    "\n  File: ", file,
    "\n  Object: ", name,
    "\n  Temporary file: ", tmp.file
  )
  return(normalizePath(tmp.file, mustWork = FALSE))
}

#' @importFrom hdf5r.Extra h5Copy h5Overwrite
.h5_overwrite_after <- function(file, out.file, name) {
  if (identical(file, out.file)) {
    return(out.file)
  }
  out.file <- h5Overwrite(file = out.file, name = name, overwrite = TRUE)
  h5Copy(
    from.file = file,
    from.name = name,
    to.file = out.file,
    to.name = name,
    overwrite = TRUE
  )
  return(out.file)
}

## Write H5AD from a list ######################################################

.write_h5ad <- function(
    object,
    file,
    name = "/",
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  .write_raw_h5ad(
    object = object,
    file = file,
    name = name,
    gzip_level = gzip_level,
    verbose = verbose
  )
  for (i in names(x = object$layers)) {
    verboseMsg("Writing layer '", i, "'")
    .write_mat_h5ad(
      mat = object$layers[[i]],
      file = file,
      name = file.path(name, "layers", i),
      gzip_level = gzip_level
    )
  }
  if ("raw" %in% names(x = object)) {
    verboseMsg("Writing 'raw'")
    .write_raw_h5ad(
      object = object$raw,
      file = file,
      name = file.path(name, "raw"),
      gzip_level = gzip_level,
      verbose = verbose
    )
  }
  verboseMsg("Writing 'obs'")
  h5Write(
    x = object$obs,
    file = file,
    name = file.path(name, "obs"),
    overwrite = TRUE,
    gzip_level = gzip_level
  )
  for (i in c("obsm", "obsp", "varm", "varp", "uns")) {
    verboseMsg("Writing ", i)
    h5Write(
      x = object[[i]],
      file = file,
      name = file.path(name, i),
      gzip_level = gzip_level
    )
  }
  return(invisible(x = NULL))
}

#' @importFrom S4Arrays is_sparse
#' @importFrom HDF5Array writeHDF5Array writeTENxMatrix
#' @importFrom BPCells write_matrix_anndata_hdf5
.write_mat_h5ad <- function(
    mat,
    file,
    name = "X",
    gzip_level = 0L,
    verbose = FALSE
) {
  name <- h5AbsLinkName(name = name)
  if (inherits(x = mat, what = "IterableMatrix")) {
    mat <- write_matrix_anndata_hdf5(
      mat = mat,
      path = file,
      group = name,
      gzip_level = gzip_level
    )
    # h5Delete(x = file, name = file.path(dirname(path = name), "obs"))
    # h5Delete(x = file, name = file.path(dirname(path = name), "var"))
    return(invisible(x = NULL))
  }
  if (is_sparse(x = mat)) {
    mat <- writeTENxMatrix(
      x = mat,
      filepath = file,
      group = name,
      level = gzip_level,
      verbose = verbose
    )
    # h5Delete(x = file, name = file.path(name, "shape"))
    # h5Delete(x = file, name = file.path(name, "barcodes"))
    # h5Delete(x = file, name = file.path(name, "genes"))
    return(invisible(x = NULL))
  }
  mat <- writeHDF5Array(
    x = mat,
    filepath = file,
    name = name,
    level = gzip_level,
    with.dimnames = FALSE,
    verbose = verbose
  )
  return(invisible(x = NULL))
}

#' @importFrom hdf5r.Extra h5Write
.write_raw_h5ad <- function(
    object,
    file,
    name = "raw",
    gzip_level = 0L,
    verbose = FALSE
) {
  name <- h5AbsLinkName(name)
  if ("X" %in% names(object)) {
    verboseMsg("Writing 'X'")
    .write_mat_h5ad(
      mat = object$X,
      file = file,
      name = file.path(name, "X"),
      gzip_level = gzip_level,
      verbose = verbose
    )
  }
  verboseMsg("Writing 'var'")
  h5Write(
    x = object$var,
    file = file,
    name = file.path(name, "var"),
    overwrite = TRUE,
    gzip_level = gzip_level
  )
  return(invisible(x = NULL))
}

## SingleCellExperiment ########################################################

#' @importFrom SummarizedExperiment assay assayNames
.prep_sce2h5ad_layers <- function(object, layers = NULL) {
  out <- list()
  layers <- assayNames(object) %>%
    intersect(layers)
  if (length(layers) == 0) {
    return(out)
  }
  for (i in layers) {
    out[[i]] <- assay(object, i = i, withDimnames = TRUE)
    if (inherits(out[[i]], "dgCMatrix")) {
      out[[i]] <- as(out[[i]], "IterableMatrix")
    }
  }
  return(out)
}

#' @importFrom SummarizedExperiment assayNames rowData
.prep_sce2h5ad_raw <- function(object, X = NULL) {
  out <- list(varm = list())
  out$var <- rowData(object) %>%
    as.data.frame()
  X <- X %||% assayNames(object)[1]
  if (length(X) == 0) {
    return(out)
  }
  out$X <- .prep_sce2h5ad_layers(object = object, layers = X)[[1]]
  return(out)
}

## Seurat ######################################################################

.prep_seurat2h5ad_raw <- function(assay, X = NULL) {
  X <- X %||% "data"
  out <- list(varm = list())
  out$var <- assay[[]] %>%
    as.data.frame()
  X <- X %||% Layers(assay)[1]
  if (length(X) == 0) {
    return(out)
  }
  out$X <- .prep_seurat2h5ad_layers(assay = assay, layers = X)[[1]]
  out$var <- out$var[rownames(out$X), , drop = FALSE]
  return(out)
}

.prep_seurat2h5ad_layers <- function(assay, layers = NULL) {
  out <- list()
  layers <- Layers(assay) %>%
    intersect(layers)
  if (length(layers) == 0) {
    return(out)
  }
  for (i in layers) {
    out[[i]] <- GetAssayData(assay, i)
    if (inherits(out[[i]], "dgCMatrix")) {
      out[[i]] <- as(out[[i]], "IterableMatrix")
    }
  }
  return(out)
}
