
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## toSCE #######################################################################

#' @importFrom SeuratObject GetAssayData VariableFeatures
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SeuratObject Assay
toSCE.Assay <- function(object, ...) {
  rdata <- object[[]] %>%
    as("DataFrame")
  sce_assays <- .SeuratAssay2SCEassays(object)
  sce <- SingleCellExperiment(assays = sce_assays, rowData = rdata)
  variableFeatures(sce) <- VariableFeatures(object)
  return(sce)
}

#' @importClassesFrom SeuratObject StdAssay
toSCE.StdAssay <- function(object, ...) {
  old_func <- getS3method(f = "toSCE", class = "Assay")
  old_func(object)
}

#' @importFrom Signac Fragments Annotation
#' @importClassesFrom Signac ChromatinAssay
toSCE.ChromatinAssay <- function(object, ...) {
  old_func <- getS3method(f = "toSCE", class = "Assay")
  sce <- old_func(object)
  .sce_to_csce(
    sce = sce,
    ranges = granges(object),
    fragments = Fragments(object),
    genome = genome(object),
    annotations = Annotation(object),
    ...
  )
}

#' @importFrom SeuratObject Assays
toSCE.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay[1] %||% DefaultAssay(object)
  if (!assay %in% Assays(object)) {
    stop("The assay you are trying to convert is not in the Seurat object")
  }

  curr_assay <- object@assays[[a]]
  sce <- toSCE(curr_assay)
  cells <- colnames(sce)
  colData(sce) <- object[[]][cells, , drop = FALSE] %>%
    as("DataFrame")
  colData(sce)$ident <- Idents(object)[cells]
  reducedDims(sce) <- .reductions2reducedDims(object, cells = cells)
  colPairs(sce) <- .graphs2colPairs(object, cells = cells)
  metadata(sce) <- Misc(object)
  return(sce)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SeuratObject Cells Key Reductions
.reductions2reducedDims <- function(
    object,
    reductions = NULL,
    cells = NULL,
    ...
) {
  all.cells <- Cells(object)
  cells <- cells %||% all.cells
  cells <- intersect(cells, all.cells)

  all.reducs <- Reductions(object)
  reductions <- reductions %||% all.reducs
  reductions <- intersect(reductions, all.reducs)
  out <- list()
  for (i in reductions) {
    reduc <- subset(object@reductions[[i]], cells = cells)
    i2 <- gsub("_*$", "", Key(reduc))
    out[[i2]] <- as(reduc, "LinearEmbeddingMatrix")
  }
  return(out)
}

#' @importFrom SeuratObject Cells Graphs
.graphs2colPairs <- function(
    object,
    graphs = NULL,
    cells = NULL,
    SelfHits = TRUE,
    ...
) {
  all.cells <- Cells(object)
  cells <- cells %||% all.cells
  cells <- intersect(cells, all.cells)
  if (identical(cells, all.cells)) {
    cells <- NULL
  }

  all.graphs <- Graphs(object)
  graphs <- graphs %||% all.graphs
  graphs <- intersect(graphs, all.graphs)
  out <- list()
  for (i in graphs) {
    out[[i]] <- as(object@graphs[[i]], "dgCMatrix")
    if (!is.null(cells)) {
      out[[i]] <- out[[i]][cells, cells, drop = FALSE]
    }
    if (SelfHits) {
      out[[i]] <- as(out[[i]], "SelfHits")
    }
  }
  gc(verbose = FALSE)
  return(out)
}

.SeuratAssay2SCEassays <- function(object, ...) {
  sce_assays <- list(
    counts = GetAssayData(object, "counts"),
    logcounts = GetAssayData(object, "data")
  )
  scale.data <- GetAssayData(object, "scale.data")
  if (identical(dimnames(sce_assays$counts), dimnames(scale.data))) {
    sce_assays[['scaledata']] <- scale.data
  }
  for (i in names(sce_assays)) {
    if (any(dim(sce_assays[[i]]) == 0)) {
      sce_assays[[i]] <- NULL
    }
  }
  return(sce_assays)
}
