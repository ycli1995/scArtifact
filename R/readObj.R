
#' @importFrom HDF5Array HDF5Array
#' @importFrom rtracklayer import
#'
#' @rdname Obj-IO
#' @export
readObj <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- switch(
    EXPR = class,
    "List" = .readObj_List(path, verbose = verbose, ...),
    "DataFrame" = .readObj_DF(path, verbose = verbose, ...),
    "HDF5Matrix" = HDF5Array(
      filepath = file.path(path, extra$file_name),
      name = extra$data_name,
      ...
    ),
    "MatrixDir" = open_matrix_dir(
      dir = file.path(path, extra$dir_name),
      ...
    ),
    "dgCMatrix" = open_matrix_dir(
      dir = file.path(path, extra$dir_name),
      ...
    ) %>%
      as("dgCMatrix"),
    "dgTMatrix" = open_matrix_dir(
      dir = file.path(path, extra$dir_name),
      ...
    ) %>%
      as("dgTMatrix"),
    "SelfHits" = open_matrix_dir(
      dir = file.path(path, extra$dir_name),
      ...
    ) %>%
      as("dgTMatrix") %>%
      as("SelfHits"),
    "GRanges" = import(file.path(path, extra$file_name), ...),
    "Seqinfo" = .readObj_DF(path, verbose = verbose, ...) %>% .df_to_seqinfo(),
    "FragmentsDir" = open_fragments_dir(
      dir = file.path(path, extra$dir_name),
      ...
    ),
    "LinearEmbeddingMatrix" = .readObj_LEM(path, verbose = verbose, ...),
    "RangedSummarizedExperiment" = .readObj_RSE(path, verbose = verbose, ...),
    "SingleCellExperiment" = .readObj_SCE(path, verbose = verbose, ...),
    "ChromExperiment" = .readObj_SCCE(path, verbose = verbose, ...),
    "MultiAssayExperiment" = .readObj_MAE(path, verbose = verbose, ...),
    "SingleCellMultiExperiment" = .readObj_SCME(path, verbose = verbose, ...),

    "Assay5" = .readObj_Assay5(path, verbose = verbose, ...),
    "Seurat" = .readObj_Seurat(path, verbose = verbose, ...),
    "DimReduc" = .readObj_DimReduc(path, verbose = verbose, ...),
    "Neighbor" = .readObj_Neighbor(path, verbose = verbose, ...),
    "Graph" = open_matrix_dir(
      dir = file.path(path, extra$dir_name),
      ...
    ) %>%
      as("dgCMatrix") %>%
      new(Class = "Graph", assay.used = extra$assay_used),

    stop("Unsupported class '", class, "' for readObj().")
  )
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## readObj for each class ######################################################

.readObj_List <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  all.dirs <- list.files(path = path, full.names = FALSE) %>%
    setdiff(.objfile)
  out <- list()
  for (i in all.dirs) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose, ...)
  }
  names(out) <- extra$names
  out
}

#' @importFrom data.table fread
.readObj_DF <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  filepath <- file.path(path, extra$file_name)
  verboseMsg(class, ": Reading ", filepath)
  df <- fread(
    filepath,
    header = TRUE,
    stringsAsFactors = FALSE,
    data.table = FALSE,
    ...
  )[, extra$colnames, drop = FALSE]
  rownames(df) <- df[, extra$index_column]
  df[[extra$index_column]] <- NULL
  for (i in names(extra$categories)) {
    df[[i]] <- factor(df[[i]], levels = extra$categories[[i]])
  }
  df
}

#' @importFrom HDF5Array HDF5Array
#' @importFrom SingleCellExperiment LinearEmbeddingMatrix
.readObj_LEM <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- list()
  sf.file <- file.path(path, "sampleFactors.h5")
  if (!file.exists(sf.file)) {
    stop("Cannot read the ", class, ", missing sampleFactors.h5")
  }
  for (i in c("samples", "features", "factors")) {
    name.file <- file.path(path, i)
    if (file.exists(name.file)) {
      verboseMsg(class, ": Reading '", i, "': ", name.file)
      out[[i]] <- readLines(name.file)
    }
  }
  verboseMsg(class, ": Reading 'sampleFactors': ", sf.file)
  out$sampleFactors <- HDF5Array(sf.file, name = "matrix") %>%
    as.matrix()
  colnames(out$sampleFactors) <- out$factors
  rownames(out$sampleFactors) <- out$samples
  fl.file <- file.path(path, "featureLoadings.h5")
  if (file.exists(fl.file)) {
    verboseMsg(class, ": Reading 'featureLoadings' ", fl.file)
    out$featureLoadings <- HDF5Array(fl.file, name = "matrix") %>%
      as.matrix()
  } else {
    out$featureLoadings <- matrix(0, nrow = 0, ncol = ncol(out$sampleFactors))
  }
  colnames(out$featureLoadings) <- out$factors
  rownames(out$featureLoadings) <- out$features
  out$factorData <- .readObj_DF(path = path, verbose = verbose, ...)
  out$metadata <- .read_metadata(path = path, class = class, verbose = verbose)
  LinearEmbeddingMatrix(
    sampleFactors = out$sampleFactors,
    featureLoadings = out$featureLoadings,
    factorData = DataFrame(out$factorData),
    metadata = out$metadata
  )
}

.readObj_RSE <- function(path, verbose = TRUE, ...) {
  out <- .readObj_RSE_list(path = path, verbose = verbose, ...)
  .list2RSE(out)
}

.readObj_SCE <- function(path, verbose = TRUE, ...) {
  out <- .readObj_SCE_list(path = path, verbose = verbose, ...)
  out <- .list2SCE(out)
  int_metadata(out)[[.path_key]] <- file_path_as_absolute(path)
  out
}

.readObj_SCCE <- function(path, verbose = TRUE, ...) {
  out <- .readObj_SCCE_list(path = path, verbose = verbose, ...)
  out <- .list2SCCE(out)
  int_metadata(out)[[.path_key]] <- file_path_as_absolute(path)
  out
}

.readObj_MAE <- function(path, verbose = TRUE, ...) {
  out <- .readObj_MAE_list(path = path, verbose = verbose, ...)
  .list2MAE(out)
}

.readObj_SCME <- function(path, verbose = TRUE, ...) {
  out <- .readObj_SCME_list(path = path, verbose = verbose, ...)
  out <- .list2SCME(out)
  int_metadata(out)[[.path_key]] <- file_path_as_absolute(path)
  out
}

.readObj_DimReduc <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- list()
  sf.file <- file.path(path, "embeddings.h5")
  if (!file.exists(sf.file)) {
    stop("Cannot read the ", class, ", missing embeddings.h5")
  }
  for (i in c("cells", "features", "column_names")) {
    name.file <- file.path(path, i)
    if (file.exists(name.file)) {
      verboseMsg(class, ": Reading '", i, "': ", name.file)
      out[[i]] <- readLines(name.file)
    }
  }
  verboseMsg(class, ": Reading 'embeddings.h5': ", sf.file)
  out$embeddings <- HDF5Array(sf.file, name = "matrix") %>%
    as.matrix()
  colnames(out$embeddings) <- out$column_names
  rownames(out$embeddings) <- out$cells

  for (i in c("loadings", "projected_loadings")) {
    fl.file <- file.path(path, paste0(i, ".h5"))
    if (file.exists(fl.file)) {
      verboseMsg(class, ": Reading '", i, "' ", fl.file)
      out[[i]] <- HDF5Array(fl.file, name = "matrix") %>%
        as.matrix()
      colnames(out[[i]]) <- out$column_names
      rownames(out[[i]]) <- out$features
    } else {
      out[[i]] <- matrix(0, 0, 0)
    }
  }
  for (i in c("misc", "jackstraw")) {
    filepath <- file.path(path, paste0(i, ".rds"))
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readRDS(filepath)
  }
  CreateDimReducObject(
    embeddings = out$embeddings,
    loadings = out$loadings,
    projected = out$projected_loadings,
    assay = extra$assay_used,
    stdev = extra$stdev %||% numeric(0L),
    key = extra$key,
    global = extra$global,
    jackstraw = out$jackstraw,
    misc = out$misc
  )
}

.readObj_Neighbor <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- list()
  name.file <- file.path(path, "cells")
  verboseMsg(class, ": Reading '", i, "': ", name.file)
  out[["cells"]] <- readLines(name.file)
  for (i in c("idx", "dist")) {
    filepath <- file.path(path, paste0(i, ".h5"))
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- HDF5Array(filepath, name = "matrix") %>%
      as.matrix()
    colnames(out[[i]]) <- out$cells
  }
  filepath <- file.path(path, "alg.rds")
  verboseMsg(class, ": Reading 'alg.idx' and 'alg.info' ", filepath)
  out[["alg"]] <- readRDS(filepath)
  new(
    Class = "Neighbor",
    nn.idx = out$idx,
    nn.dist = out$dist,
    alg.idx = out$alg$alg.idx,
    alg.info = out$alg$alg.info,
    cell.names = out$cells
  )
}

.readObj_Assay5 <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  assay.orig <- extra$assay_orig
  if (length(assay.orig) == 0) {
    assay.orig <- character(0)
  }
  out <- list()

  for (i in c("layers", "meta.data")) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose, ...)
  }
  for (i in c("cells", "features")) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readLines(filepath)
  }
  filepath <- file.path(path, .misc_file)
  verboseMsg(class, ": Reading 'misc' ", filepath)
  out[['misc']] <- readRDS(filepath)
  object <- .CreateStdAssay(
    counts = out$layers,
    cells = out$cells,
    features = out$features,
    type = "Assay5"
  )
  object@meta.data <- out$meta.data[Features(object), ]
  object@misc <- out$misc
  Key(object) <- extra$key
  DefaultAssay(object) <- assay.orig
  DefaultLayer(object) <- extra$default
  return(object)
}

.readObj_Seurat <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  slot.names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  out <- list()
  for (i in slot.names) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose)
  }
  cells <- rownames(out$meta.data)
  for (i in c("misc", "tools", "commands")) {
    filepath <- file.path(path, paste0(i, ".rds"))
    if (file.exists(filepath)) {
      out[[i]] <- readRDS(filepath)
    } else {
      out[[i]] <- list()
    }
  }
  suppressWarnings(new(
    Class = "Seurat",
    assays = out$assays,
    meta.data = out$meta.data,
    active.assay = extra$active_assay,
    active.ident = out$meta.data$active.ident,
    graphs = out$graphs %||% list(),
    neighbors = out$neighbors %||% list(),
    reductions = out$reductions %||% list(),
    images = out$images %||% list(),
    project.name = extra$project_name,
    misc = out$misc,
    version = package_version(extra$seuratobject_version),
    commands = out$commands,
    tools = out$tools
  ))
}



## read S4 object as list ######################################################

.readObj_RSE_list <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  slot.names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  to.read <- c("rowData", "colData", "rowRanges") %>%
    intersect(slot.names)
  out <- list()
  for (i in to.read) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose)
  }
  filepath <- file.path(path, "assays")
  verboseMsg(class, ": Reading 'assays' ", filepath)
  out$assays <- .read_Exps(
    filepath,
    rownames = rownames(rowData),
    colnames = rownames(colData),
    verbose = verbose,
    ...
  )
  out$metadata <- .read_metadata(path = path, class = class, verbose = verbose)
  return(out)
}

.readObj_SCE_list <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- .readObj_RSE_list(path = path, verbose = verbose)
  filepath <- file.path(path, "altExps")
  verboseMsg(class, ": Reading 'altExps' ", filepath)
  out$altExps <- .read_Exps(
    filepath,
    key = "altExp",
    colnames = rownames(out$colData),
    verbose = verbose,
    ...
  )
  slot.names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  to.read <- c(.red_key, .colp_key, .rowp_key) %>%
    intersect(slot.names)
  for (i in to.read) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose)
  }
  return(out)
}

.readObj_SCCE_list <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- .readObj_SCE_list(path = path, verbose = verbose)
  slot.names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  to.read <- c(.frag_key, .sinfo_key, .annot_key) %>%
    intersect(slot.names)
  for (i in to.read) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose)
  }
  return(out)
}

.readObj_MAE_list <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  slot.names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  to.read <- c("colData", "sampleMap") %>%
    intersect(slot.names)
  out <- list()
  for (i in to.read) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose)
  }
  filepath <- file.path(path, "experiments")
  verboseMsg(class, ": Reading 'experiments' ", filepath)
  out$experiments <- .read_Exps(
    filepath,
    key = "experiment",
    verbose = verbose,
    ...
  )
  out$metadata <- .read_metadata(path = path, class = class, verbose = verbose)
  return(out)
}

.readObj_SCME_list <- function(path, verbose = TRUE, ...) {
  extra <- readObjFile(path)
  class <- extra$class
  out <- .readObj_MAE_list(path = path, verbose = verbose)

  slot.names <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  to.read <- c(.red_key, .colp_key) %>%
    intersect(slot.names)
  for (i in to.read) {
    filepath <- file.path(path, i)
    verboseMsg(class, ": Reading '", i, "' ", filepath)
    out[[i]] <- readObj(filepath, verbose = verbose)
  }
  out$defaultExp <- extra$defaultExp
  return(out)
}

## reading helpers #############################################################

.read_Exps <- function(
    path,
    key = "assay",
    rownames = NULL,
    colnames = NULL,
    verbose = TRUE,
    ...
) {
  assays <- .readObj_List(path = path, verbose = verbose)
  nms <- names(assays)
  for (i in seq_along(assays)) {
    verboseMsg("Checking dimension names for ", key, " '",  nms[i], "' ")
    assays[[i]] <- .check_set_rownames(assays[[i]], rownames)
    assays[[i]] <- .check_set_colnames(assays[[i]], colnames)
  }
  return(assays)
}

.read_metadata <- function(
    path,
    class,
    name = NULL,
    key = "metadata",
    verbose = TRUE
) {
  name <- name %||% .md_file
  out <- list()
  md.file <- file.path(path, name)
  if (!file.exists(md.file)) {
    return(out)
  }
  verboseMsg(class, ": Reading '", key, "': ", md.file)
  return(readRDS(md.file))
}

## Convert S3 to S4 ############################################################

#' @importFrom GenomeInfoDb Seqinfo
.df_to_seqinfo <- function(df) {
  seqlengths <- df[['seqlengths']] %||% NA
  isCircular <- df[['isCircular']] %||% NA
  genome <- df[['genome']] %||% NA
  return(Seqinfo(
    seqnames = rownames(df),
    seqlengths = seqlengths,
    isCircular = isCircular,
    genome = genome
  ))
}

#' @importFrom SummarizedExperiment rowData<- rowRanges<- SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
.list2RSE <- function(rlist) {
  out <- SummarizedExperiment(
    assays = rlist$assays,
    colData = as(rlist$colData, "DataFrame"),
    metadata = rlist$metadata
  )
  if (length(rlist$rowRanges) > 0) {
    rowRanges(out) <- rlist$rowRanges
  }
  rowData(out) <- as(rlist$rowData, "DataFrame")
  return(out)
}

#' @importFrom SummarizedExperiment rowData<- rowRanges<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
.list2SCE <- function(rlist) {
  out <- SingleCellExperiment(
    assays = rlist$assays,
    colData = as(rlist$colData, "DataFrame"),
    reducedDims = rlist[[.red_key]],
    colPairs = rlist[[.colp_key]],
    rowPairs = rlist[[.rowp_key]],
    altExps = rlist[[.alt_key]],
    metadata = rlist$metadata
  )
  if (length(rlist$rowRanges) > 0) {
    rowRanges(out) <- rlist$rowRanges
  }
  rowData(out) <- as(rlist$rowData, "DataFrame")
  return(out)
}

.list2SCCE <- function(rlist) {
  out <- SingleCellExperiment(
    assays = rlist$assays,
    colData = as(rlist$colData, "DataFrame"),
    reducedDims = rlist[[.red_key]],
    colPairs = rlist[[.colp_key]],
    rowPairs = rlist[[.rowp_key]],
    altExps = rlist[[.alt_key]],
    metadata = rlist$metadata
  )
  .sce_to_csce(
    out,
    ranges = rlist$rowRanges,
    genome = rlist[[.sinfo_key]],
    fragments = rlist[[.frag_key]],
    annotations = rlist[[.annot_key]]
  )
}

#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
.list2MAE <- function(rlist) {
  MultiAssayExperiment(
    experiments = ExperimentList(rlist$experiments),
    colData = as(rlist$colData, "DataFrame"),
    sampleMap = as(rlist$sampleMap, "DataFrame"),
    metadata = rlist$metadata
  )
}

#' @importFrom MultiAssayExperiment ExperimentList
.list2SCME <- function(rlist) {
  if (rlist$defaultExp == "") {
    rlist$defaultExp <- NA
  }
  if (length(rlist$defaultExp) == 0) {
    rlist$defaultExp <- NULL
  }
  SingleCellMultiExperiment(
    experiments = ExperimentList(rlist$experiments),
    sampleMap = as(rlist$sampleMap, "DataFrame"),
    colData = as(rlist$colData, "DataFrame"),
    defaultExp = rlist$defaultExp,
    reducedDims = rlist[[.red_key]],
    colPairs = rlist[[.colp_key]],
    metadata = rlist$metadata
  )
}







