
#' @rdname Obj-IO
#' @export
#' @method dumpObj default
dumpObj.default <- function(object, path, verbose = TRUE, ...) {
  if (is.null(object)) {
    return(invisible(NULL))
  }
  stop("Saving ", class(object), " to disk is not supported.")
}

#' @rdname Obj-IO
#' @export
#' @method dumpObj list
dumpObj.list <- function(object, path, verbose = TRUE, class = "List", ...) {
  verboseMsg(class, ": ", path)
  nms <- names(object)
  for (i in seq_along(object)) {
    nm <- gsub("\\s|\\/", "_", nms[i])
    filepath <- file.path(path, nm)
    verboseMsg(class, ": Writing element '", nms[i], "' to ", filepath)
    dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
    filepath <- file_path_as_absolute(filepath)
    dumpObj(object[[i]], path = filepath, verbose = verbose, ...)
  }
  extra <- prepInfo(object)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom S4Vectors SimpleList
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj SimpleList
dumpObj.SimpleList <- function(
    object,
    path,
    verbose = TRUE,
    class = "List",
    ...
) {
  old_func <- getS3method(f = "dumpObj", class = "list")
  old_func(object, path = path, verbose = verbose, class = class, ...)
}

#' @importFrom data.table fwrite
#' @importFrom tibble rownames_to_column
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj data.frame
dumpObj.data.frame <- function(
    object,
    path,
    verbose = TRUE,
    class = "DataFrame",
    file.name = "DataFrame.csv",
    idx.name = "_index",
    ...
) {
  verboseMsg(class, ": ", path)
  filepath <- file.path(path, file.name)
  verboseMsg(class, ": Writing the data.frame to ", filepath)
  fwrite(
    object %>% rownames_to_column(idx.name),
    file = filepath,
    row.names = FALSE,
    ...
  )
  extra <- prepInfo(object, file.name = file.name, idx.name = idx.name)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom S4Vectors DataFrame
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj DataFrame
dumpObj.DataFrame <- function(
    object,
    path,
    verbose = TRUE,
    class = "DataFrame",
    file.name = "DataFrame.csv",
    idx.name = "_index",
    ...
) {
  dumpObj(
    as.data.frame(object),
    path = path,
    verbose = verbose,
    class = class,
    file.name = file.name,
    idx.name = idx.name,
    ...
  )
}

#' @importFrom HDF5Array writeHDF5Array
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj matrix
dumpObj.matrix <- function(
    object,
    path,
    verbose = TRUE,
    class = "HDF5Matrix",
    file.name = "matrix.h5",
    name = "matrix",
    with.dimnames = TRUE,
    ...
) {
  verboseMsg(class, ": ", path)
  filepath <- file.path(path, file.name)
  verboseMsg(class, ": Writing the ", class(object)[1], " to ", filepath)
  object <- writeHDF5Array(
    x = object,
    filepath = filepath,
    name = name,
    verbose = FALSE,
    with.dimnames = with.dimnames,
    ...
  )
  extra <- prepInfo(
    object,
    file.name = file.name,
    name = name,
    with.dimnames = with.dimnames
  )
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom DelayedArray DelayedArray
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj DelayedArray
dumpObj.DelayedArray <- function(
    object,
    path,
    verbose = TRUE,
    class = "HDF5Matrix",
    file.name = "matrix.h5",
    name = "matrix",
    with.dimnames = TRUE,
    ...
) {
  old_func <- getS3method(f = "dumpObj", class = "matrix")
  old_func(
    object,
    path = path,
    verbose = verbose,
    file.name = file.name,
    name = name,
    class = class,
    ...
  )
}

#' @rdname Obj-IO
#' @export
#' @method dumpObj IterableMatrix
dumpObj.IterableMatrix <- function(
    object,
    path,
    verbose = TRUE,
    dir.name = "matrix",
    class = "MatrixDir",
    ...
) {
  verboseMsg(class, ": ", path)
  filepath <- file.path(path, dir.name)
  verboseMsg(class, ": Writing ", class(object)[1], " to ", filepath)
  object <- write_matrix_dir(mat = object, dir = filepath, overwrite = TRUE)
  extra <- prepInfo(object, dir.name = dir.name)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom Matrix CsparseMatrix generalMatrix
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj CsparseMatrix
dumpObj.CsparseMatrix <- function(
    object,
    path,
    verbose = TRUE,
    dir.name = "matrix",
    class = "dgCMatrix",
    ...
) {
  verboseMsg(class, ": ", path)
  filepath <- file.path(path, dir.name)
  verboseMsg(class, ": Writing ", class(object)[1], " to ", filepath)
  if (!inherits(object, "generalMatrix")) {
    object <- as(object, Class = "generalMatrix")
  }
  mat <- as(object, "IterableMatrix")
  mat <- write_matrix_dir(mat = mat, dir = filepath, overwrite = TRUE)
  extra <- prepInfo(mat, dir.name = dir.name)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom Matrix TsparseMatrix generalMatrix
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj TsparseMatrix
dumpObj.TsparseMatrix <- function(
    object,
    path,
    verbose = TRUE,
    dir.name = "matrix",
    class = "dgTMatrix",
    ...
) {
  dumpObj.CsparseMatrix(
    object = object,
    path = path,
    verbose = verbose,
    dir.name = dir.name,
    class = class,
    ...
  )
}

#' @rdname Obj-IO
#' @export
#' @method dumpObj IterableFragments
dumpObj.IterableFragments <- function(
    object,
    path,
    verbose = TRUE,
    dir.name = "fragments",
    class = "FragmentsDir",
    ...
) {
  verboseMsg(class, ": ", path)
  filepath <- file.path(path, dir.name)
  verboseMsg(class, ": Writing ", class(object)[1], " to ", filepath)
  out <- write_fragments_dir(
    fragments = object,
    dir = filepath,
    overwrite = TRUE
  )
  extra <- prepInfo(object)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @rdname Obj-IO
#' @export
#' @method dumpObj SelfHits
dumpObj.SelfHits <- function(
    object,
    path,
    verbose = TRUE,
    dir.name = "matrix",
    class = "SelfHits",
    ...
) {
  object %>%
    as("TsparseMatrix") %>%
    dumpObj(
      path = path,
      verbose = verbose,
      dir.name = dir.name,
      class = class,
      ...
    )
}

#' @importClassesFrom GenomeInfoDb Seqinfo
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj Seqinfo
dumpObj.Seqinfo <- function(
    object,
    path,
    verbose = TRUE,
    class = "Seqinfo",
    ...
) {
  .S4_dumpObj(
    object = object,
    path = path,
    verbose = verbose,
    class = class,
    file.name = "Seqinfo.csv",
    ...
  )
}

#' @importFrom rtracklayer export
#' @importClassesFrom GenomicRanges GRanges
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj GRanges
dumpObj.GRanges <- function(
    object,
    path,
    verbose = TRUE,
    file.name = "GRanges.gtf",
    class = "GRanges",
    ...
) {
  verboseMsg(class, ": ", path)
  filepath <- file.path(path, file.name)
  verboseMsg(class, ": Writing ", class(object)[1], " to ", filepath)
  export(object, filepath, format = "gtf")
  extra <- prepInfo(object, file.name = file.name)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importFrom HDF5Array writeHDF5Array
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj LinearEmbeddingMatrix
dumpObj.LinearEmbeddingMatrix <- function(
    object,
    path,
    verbose = TRUE,
    class = "LinearEmbeddingMatrix",
    ...
) {
  if (any(dim(object) == 0)) {
    warning("Skip writing an empty ", class(object)[1], immediate. = TRUE)
    return(invisible(NULL))
  }
  out <- prepObj(object)
  verboseMsg(class, ": ", path)
  for (i in c("sampleFactors", "featureLoadings")) {
    if (!is_empty(out[[i]])) {
      filepath <- file.path(path, paste0(i, ".h5"))
      verboseMsg(class, ": Writing '", i, "' to ", filepath)
      out[[i]] <- writeHDF5Array(
        out[[i]],
        filepath = filepath,
        name = "matrix",
        with.dimnames = FALSE,
        level = 0L,
        verbose = FALSE
      )
    }
  }
  for (i in c("samples", "features", "factors")) {
    if (length(out[[i]]) > 0) {
      filepath <- file.path(path, i)
      verboseMsg(class, ": Writing '", i, "' to ", filepath)
      write(out[[i]], filepath)
    }
  }
  filepath <- file.path(path, .md_file)
  verboseMsg(class, ": Writing 'metadata' to ", filepath)
  saveRDS(out$metadata, file = filepath)
  verboseMsg(
    class, ": Writing 'factorData' to ",
    file.path(path, "factorData.csv")
  )
  dumpObj(
    object = out$factorData,
    path = path,
    verbose = verbose,
    file.name = "factorData.csv",
    class = class
  )
  return(invisible(NULL))
}

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj RangedSummarizedExperiment
dumpObj.RangedSummarizedExperiment <- function(
    object,
    path,
    verbose = TRUE,
    class = "RangedSummarizedExperiment",
    ...
) {
  .S4_dumpObj(
    object = object,
    path = path,
    class = class,
    verbose = verbose,
    ...
  )
}

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj SingleCellExperiment
dumpObj.SingleCellExperiment <- function(
    object,
    path,
    verbose = TRUE,
    class = "SingleCellExperiment",
    ...
) {
  .S4_dumpObj(
    object = object,
    path = path,
    class = class,
    verbose = verbose,
    ...
  )
}

#' @rdname Obj-IO
#' @export
#' @method dumpObj ChromExperiment
dumpObj.ChromExperiment <- function(
    object,
    path,
    verbose = TRUE,
    class = "ChromExperiment",
    ...
) {
  .S4_dumpObj(
    object = object,
    path = path,
    class = class,
    verbose = verbose,
    ...
  )
}

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj MultiAssayExperiment
dumpObj.MultiAssayExperiment <- function(
    object,
    path,
    verbose = TRUE,
    class = "MultiAssayExperiment",
    ...
) {
  .S4_dumpObj(
    object = object,
    path = path,
    class = class,
    verbose = verbose,
    ...
  )
}

#' @rdname Obj-IO
#' @export
#' @method dumpObj SingleCellMultiExperiment
dumpObj.SingleCellMultiExperiment <- function(
    object,
    path,
    verbose = TRUE,
    class = "SingleCellMultiExperiment",
    ...
) {
  .S4_dumpObj(
    object = object,
    path = path,
    class = class,
    verbose = verbose,
    ...
  )
}

#' @importClassesFrom SeuratObject Assay5
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj Assay5
dumpObj.Assay5 <- function(
    object,
    path,
    verbose = TRUE,
    class = "Assay5",
    ...
) {
  verboseMsg(class, ": ", path)
  extra <- prepInfo(object)
  object <- prepObj(object)
  tmp.slots <- c("cells", "features", "misc")
  tmp <- object[tmp.slots]
  object <- object[setdiff(names(object), tmp.slots)]
  dumpObj(
    object = object,
    path = path,
    verbose = verbose,
    class = class,
    ...
  )
  for (i in c("cells", "features")) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Writing '", i, "' to ", filepath)
    write(tmp[[i]], file = filepath)
  }
  filepath <- file.path(path, .misc_file)
  verboseMsg(class, ": Writing 'misc' to ", filepath)
  saveRDS(tmp[['misc']], file = filepath)
  writeObjFile(path = path, class = class, extra = extra)
}

#' @importClassesFrom SeuratObject DimReduc
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj DimReduc
dumpObj.DimReduc <- function(
    object,
    path,
    verbose = TRUE,
    class = "DimReduc",
    ...
) {
  out <- prepObj(object)
  if (is_empty(out$embeddings)) {
    warning("Skip writing an empty ", class(object)[1], immediate. = TRUE)
    return(invisible(NULL))
  }
  verboseMsg(class, ": ", path)
  for (i in c("embeddings", "loadings", "projected_loadings")) {
    if (!is_empty(out[[i]])) {
      filepath <- file.path(path, paste0(i, ".h5"))
      verboseMsg(class, ": Writing '", i, "' to ", filepath)
      out[[i]] <- writeHDF5Array(
        out[[i]],
        filepath = filepath,
        name = "matrix",
        with.dimnames = FALSE,
        level = 0L,
        verbose = FALSE
      )
    }
  }
  for (i in c("cells", "features", "column_names")) {
    if (length(out[[i]]) > 0) {
      filepath <- file.path(path, i)
      verboseMsg(class, ": Writing '", i, "' to ", filepath)
      write(out[[i]], filepath)
    }
  }
  for (i in c("jackstraw", "misc")) {
    filepath <- file.path(path, paste0(i, ".rds"))
    verboseMsg(class, ": Writing '", i, "' to ", filepath)
    saveRDS(out[[i]], file = filepath)
  }
  extra <- prepInfo(object)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom SeuratObject Graph
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj Graph
dumpObj.Graph <- function(
    object,
    path,
    verbose = TRUE,
    dir.name = "matrix",
    class = "Graph",
    ...
) {
  dumpObj.CsparseMatrix(
    as(object, "dgCMatrix"),
    path = path,
    verbose = verbose,
    dir.name = dir.name,
    class = class,
    ...
  )
  extra <- prepInfo(object)
  extra$dir_name <- dir.name
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom SeuratObject Neighbor
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj Neighbor
dumpObj.Neighbor <- function(
    object,
    path,
    verbose = TRUE,
    class = "Neighbor",
    ...
) {
  verboseMsg(class, ": ", path)
  out <- prepObj(object)
  for (i in c("idx", "dist")) {
    if (!is_empty(out[[i]])) {
      filepath <- file.path(path, paste0(i, ".h5"))
      verboseMsg(class, ": Writing '", i, "' to ", filepath)
      out[[i]] <- writeHDF5Array(
        out[[i]],
        filepath = filepath,
        name = "matrix",
        with.dimnames = FALSE,
        level = 0L,
        verbose = FALSE
      )
    }
  }
  if (length(out[["cells"]]) > 0) {
    filepath <- file.path(path, "cells")
    verboseMsg(class, ": Writing '", i, "' to ", filepath)
    write(out[["cells"]], filepath)
  }
  for (i in "alg") {
    filepath <- file.path(path, paste0(i, ".rds"))
    verboseMsg(class, ": Writing '", i, "' to ", filepath)
    saveRDS(out[[i]], file = filepath)
  }
  extra <- prepInfo(object)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

#' @importClassesFrom SeuratObject Seurat
#'
#' @rdname Obj-IO
#' @export
#' @method dumpObj Seurat
dumpObj.Seurat <- function(
    object,
    path,
    verbose = TRUE,
    class = "Seurat",
    ...
) {
  verboseMsg(class, ": ", path)
  extra <- prepInfo(object)
  object <- prepObj(object)
  tmp.slots <- c("misc", "commands", "tools")
  tmp <- object[tmp.slots]
  object <- object[setdiff(names(object), tmp.slots)]
  dumpObj(
    object = object,
    path = path,
    verbose = verbose,
    class = class,
    ...
  )
  for (i in tmp.slots) {
    filepath <- file.path(path, paste0(i, ".rds"))
    verboseMsg(class, ": Writing '", i, "' to ", filepath)
    saveRDS(tmp[[i]], file = filepath)
  }
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

## prepInfo ####################################################################

#' @rdname ObjFile-IO
#' @export
#' @method prepInfo default
prepInfo.default <- function(object, ...) {
  list(version = .writeobj_ver)
}

#' @rdname ObjFile-IO
#' @export
#' @method prepInfo list
prepInfo.list <- function(object, ...) {
  list(version = .writeobj_ver, names = names(object))
}

#' @importClassesFrom S4Vectors SimpleList
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo SimpleList
prepInfo.SimpleList <- function(object, ...) {
  list(version = .writeobj_ver, names = names(object))
}

#' @rdname ObjFile-IO
#' @export
#' @method prepInfo data.frame
prepInfo.data.frame <- function(
    object,
    file.name = "DataFrame.csv",
    idx.name = "_index",
    ...
) {
  factor_cols <- object %>%
    Filter(f = is.factor) %>%
    names()
  categories <- sapply(
    factor_cols,
    function(i) levels(object[[i]]),
    simplify = FALSE
  )
  list(
    version = .writeobj_ver,
    file_name = file.name,
    colnames = c(idx.name, colnames(object)),
    index_column = idx.name,
    categories = categories
  )
}

#' @importClassesFrom HDF5Array HDF5Array
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo HDF5Array
prepInfo.HDF5Array <- function(
    object,
    file.name = "matrix.h5",
    name = "matrix",
    with.dimnames = TRUE,
    ...
) {
  list(
    version = .writeobj_ver,
    file_name = file.name,
    data_name = name,
    with_dimnames = with.dimnames
  )
}

#' @rdname ObjFile-IO
#' @export
#' @method prepInfo IterableMatrix
prepInfo.IterableMatrix <- function(object, dir.name = "matrix", ...) {
  list(version = .writeobj_ver, dir_name = dir.name)
}

#' @rdname ObjFile-IO
#' @export
#' @method prepInfo IterableFragments
prepInfo.IterableFragments <- function(object, dir.name = "fragments", ...) {
  list(version = .writeobj_ver, dir_name = dir.name)
}

#' @importClassesFrom GenomicRanges GRanges
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo GRanges
prepInfo.GRanges <- function(object, file.name = "GRanges.gtf", ...) {
  list(version = .writeobj_ver, file_name = file.name)
}

#' @rdname ObjFile-IO
#' @export
#' @method prepInfo SingleCellMultiExperiment
prepInfo.SingleCellMultiExperiment <- function(object, ...) {
  defaultExp <- object@defaultExp
  if (is.na(defaultExp)) {
    defaultExp <- ""
  }
  list(version = .writeobj_ver, defaultExp = defaultExp)
}

#' @importFrom SeuratObject DefaultLayer Key
#' @importClassesFrom SeuratObject StdAssay
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo StdAssay
prepInfo.StdAssay <- function(object, ...) {
  list(
    key = Key(object),
    assay_orig = object@assay.orig,
    default = DefaultLayer(object)
  )
}

#' @importClassesFrom SeuratObject Graph
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo Graph
prepInfo.Graph <- function(object, ...) {
  list(assay_used = object@assay.used)
}

#' @importFrom SeuratObject IsGlobal Key Stdev
#' @importClassesFrom SeuratObject DimReduc
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo DimReduc
prepInfo.DimReduc <- function(object, ...) {
  out <- list(
    assay_used = object@assay.used,
    global = IsGlobal(object),
    key = Key(object)
  )
  if (length(Stdev(object)) > 0) {
    out$stdev <- Stdev(object)
  }
  return(out)
}

#' @importFrom SeuratObject DefaultAssay Project Version
#' @importClassesFrom SeuratObject DimReduc
#'
#' @rdname ObjFile-IO
#' @export
#' @method prepInfo Seurat
prepInfo.Seurat <- function(object, ...) {
  list(
    active_assay = DefaultAssay(object),
    project_name = Project(object),
    seuratobject_version = as.character(Version(object))
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.md_file <- "metadata.rds"

.misc_file <- "misc.rds"

.S4_dumpObj <- function(
    object,
    path,
    class,
    verbose = TRUE,
    ...
) {
  verboseMsg(class, ": ", path)
  extra <- prepInfo(object)
  object <- prepObj(object)
  md <- object$metadata
  object$metadata <- NULL
  dumpObj(
    object = object,
    path = path,
    verbose = verbose,
    class = class,
    ...
  )
  filepath <- file.path(path, .md_file)
  verboseMsg(class, ": Writing 'metadata' to ", filepath)
  saveRDS(md, file = filepath)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}

.Seurat_dumpObj <- function(
    object,
    path,
    class,
    verbose = TRUE,
    ...
) {
  verboseMsg(class, ": ", path)
  extra <- prepInfo(object)
  object <- prepObj(object)
  misc <- object$misc
  object$misc <- NULL
  dumpObj(
    object = object,
    path = path,
    verbose = verbose,
    class = class,
    ...
  )
  filepath <- file.path(path, .misc_file)
  verboseMsg(class, ": Writing 'misc' to ", misc.path)
  saveRDS(misc, file = filepath)
  writeObjFile(path = path, class = class, extra = extra)
  return(invisible(NULL))
}
